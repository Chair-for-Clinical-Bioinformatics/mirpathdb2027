#!/usr/bin/env python3
"""Compute the authoritative machine-readable release census (REL-001/REL-002).

``release_census.json`` is derived DIRECTLY from a finalized release payload —
never from manifests, database state, or hand-maintained constants — and is the
single source every census consumer (api/stats, MCP data status, landing tiles,
download registry, release manifest, manuscript/supplement macros, release
notes) must be generated or populated from::

    final release files -> release_census.json -> API / MCP / website / LaTeX / manifests

Definitions (fixed by contract):

- *retained*: every row of ``enrichment_significant/<org>/<FAMILY>_<TIER>.tsv``
  (the shipped store, BH-FDR <= 0.25).
- *significant*: retained rows with ``p_adjust`` <= 0.05; *stringent*: <= 0.005.
  ``p_adjust`` is always the PRIMARY column: the full-family BH value in a
  contract-v2 payload, the overlap-family value in a legacy (v1) payload. A v2
  payload additionally carries the legacy ``p_adjust_overlap`` column; its
  ``<= 0.05`` count is reported ONLY as a within-retained-set figure
  (``significant_overlap_family_within_retained_total``), because the retained
  set itself is selected on ``p_adjust`` and is not the legacy retained set.
- ``mirnas_retained``: distinct miRNA identifiers present in the retained
  association files (REL-002: this is 6,430 for the audited build — the miRNA
  entity table of the webserver is a DIFFERENT population and must be exposed
  under its own name, never as ``n_mirnas_retained``).
- *pathway accessions*: distinct pathway identifiers (union across organisms);
  *organism-resolved pathway records*: the same identifier counted once per
  organism that serves it.

Run::

    python scripts/build_release_census.py --payload <release-dir> \
        --out <release-dir>/release_census.json

The output is deterministic for a given payload (keys sorted, no timestamps in
the counted body; generation metadata lives under ``_meta``).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ALPHA = 0.05
STRINGENT_ALPHA = 0.005
TIERS = ("HIGH", "MEDIUM", "LOW")
LEGACY_P_ADJUST_COLUMN = "p_adjust_overlap"
LEGACY_ADJUSTMENT_FAMILY = "overlap_family_bh (legacy)"
PRIMARY_P_ADJUST_FULL = "p_adjust (full-family BH over m_family eligible terms)"
PRIMARY_P_ADJUST_LEGACY = "p_adjust (overlap-family BH, legacy)"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_census(payload: Path) -> dict:
    """Stream the payload's enrichment store and label table into the census."""
    enr = payload / "enrichment_significant"
    if not enr.is_dir():
        raise SystemExit(f"not a release payload (no enrichment_significant/): {payload}")

    cells_per_org: dict[str, int] = defaultdict(int)
    rows_ret = defaultdict(int)  # org -> retained rows
    rows_sig = defaultdict(int)
    rows_str = defaultdict(int)
    fam_ret = defaultdict(int)  # family -> rows
    fam_sig = defaultdict(int)
    org_fam_ret = defaultdict(int)  # (org, family) -> rows
    org_fam_sig = defaultdict(int)
    tier_ret = defaultdict(int)  # tier -> rows
    tier_sig = defaultdict(int)
    mirnas_ret: dict[str, set] = defaultdict(set)  # org -> mirnas
    mirnas_sig: dict[str, set] = defaultdict(set)
    paths_ret: dict[str, set] = defaultdict(set)  # org -> pathway ids
    paths_sig: dict[str, set] = defaultdict(set)
    org_fam_paths: dict[tuple, set] = defaultdict(set)
    # Contract v2 only: legacy overlap-family significance WITHIN the retained set,
    # and the family size M per cell (must be one value per cell and per family).
    rows_sig_legacy = defaultdict(int)
    m_family_cell: dict[tuple, set] = defaultdict(set)  # (org, family, tier) -> {M}
    legacy_cells = 0

    for cell in sorted(enr.glob("*/*.tsv")):
        org = cell.parent.name
        stem = cell.stem  # e.g. GO_BP_HIGH / WIKIPATHWAYS_MEDIUM
        tier = next(t for t in TIERS if stem.endswith(f"_{t}"))
        family = stem[: -(len(tier) + 1)]
        cells_per_org[org] += 1
        with cell.open(encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = next(reader)
            i_mirna, i_path = header.index("mirna"), header.index("pathway_id")
            i_padj = header.index("p_adjust")
            i_legacy = i_mfam = None
            if LEGACY_P_ADJUST_COLUMN in header:
                i_legacy = header.index(LEGACY_P_ADJUST_COLUMN)
                i_mfam = header.index("m_family")
                legacy_cells += 1
            for row in reader:
                mirna, pid, padj = row[i_mirna], row[i_path], float(row[i_padj])
                if i_legacy is not None:
                    if float(row[i_legacy]) <= ALPHA:
                        rows_sig_legacy[org] += 1
                    m_family_cell[(org, family, tier)].add(int(row[i_mfam]))
                rows_ret[org] += 1
                fam_ret[family] += 1
                org_fam_ret[(org, family)] += 1
                tier_ret[tier] += 1
                mirnas_ret[org].add(mirna)
                paths_ret[org].add(pid)
                org_fam_paths[(org, family)].add(pid)
                if padj <= ALPHA:
                    rows_sig[org] += 1
                    fam_sig[family] += 1
                    org_fam_sig[(org, family)] += 1
                    tier_sig[tier] += 1
                    mirnas_sig[org].add(mirna)
                    paths_sig[org].add(pid)
                    if padj <= STRINGENT_ALPHA:
                        rows_str[org] += 1

    labels_path = payload / "pathway_labels.tsv"
    n_labels = n_inferred = 0
    with labels_path.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_labels += 1
            if r.get("inferred") == "true":
                n_inferred += 1

    orgs = sorted(cells_per_org)
    union_mirnas_ret = set().union(*mirnas_ret.values())
    union_mirnas_sig = set().union(*mirnas_sig.values())
    union_paths = set().union(*paths_ret.values())
    payload_files = sorted(p for p in payload.rglob("*") if p.is_file())

    if legacy_cells not in (0, sum(cells_per_org.values())):
        raise SystemExit(
            f"mixed contract layouts: {legacy_cells} of {sum(cells_per_org.values())} "
            f"cells carry {LEGACY_P_ADJUST_COLUMN}"
        )
    has_legacy = legacy_cells > 0
    m_family = _m_family_per_organism_per_family(m_family_cell) if has_legacy else None
    info = _release_info(payload)

    associations = {
        "retained_total": sum(rows_ret.values()),
        "retained_per_organism": {o: rows_ret[o] for o in orgs},
        "significant_total": sum(rows_sig.values()),
        "significant_per_organism": {o: rows_sig[o] for o in orgs},
        "stringent_total": sum(rows_str.values()),
        "stringent_per_organism": {o: rows_str[o] for o in orgs},
    }
    if has_legacy:
        # Rows with the LEGACY overlap-family value <= alpha, counted only among the
        # rows retained on the primary p_adjust -- not the legacy release's count.
        associations["significant_overlap_family_within_retained_total"] = sum(
            rows_sig_legacy.values()
        )
        associations["significant_overlap_family_within_retained_per_organism"] = {
            o: rows_sig_legacy[o] for o in orgs
        }

    census = {
        "_meta": {
            "generator": "scripts/build_release_census.py",
            "computed_from": payload.name,  # package basename only, no host path
            "release_id": info.get("release_id"),
            "contract_version": _contract_version(payload),
            "adjustment_family": info.get("adjustment_family") or LEGACY_ADJUSTMENT_FAMILY,
            "primary_p_adjust": info.get("primary_p_adjust")
            or (PRIMARY_P_ADJUST_FULL if has_legacy else PRIMARY_P_ADJUST_LEGACY),
            "sha256sums_digest": _sha256(payload / "SHA256SUMS")
            if (payload / "SHA256SUMS").exists()
            else None,
            "alpha": ALPHA,
            "stringent_alpha": STRINGENT_ALPHA,
            "retained_definition": "every shipped enrichment row (BH-FDR <= 0.25 store)",
        },
        "organisms": orgs,
        "cells_per_organism": dict(cells_per_org),
        "cells_total": sum(cells_per_org.values()),
        "associations": associations,
        "mirnas": {
            "retained_distinct_union": len(union_mirnas_ret),
            "retained_distinct_per_organism": {o: len(mirnas_ret[o]) for o in orgs},
            "significant_distinct_union": len(union_mirnas_sig),
            "significant_distinct_per_organism": {o: len(mirnas_sig[o]) for o in orgs},
        },
        "pathways": {
            "distinct_accessions_union": len(union_paths),
            "organism_resolved_records": sum(len(paths_ret[o]) for o in orgs),
            "organism_resolved_records_per_organism": {o: len(paths_ret[o]) for o in orgs},
            "records_with_significant": sum(len(paths_sig[o]) for o in orgs),
            "records_with_significant_per_organism": {o: len(paths_sig[o]) for o in orgs},
            "per_organism_per_family": {
                o: {
                    fam: len(org_fam_paths[(o, fam)])
                    for fam in sorted({f for (oo, f) in org_fam_paths if oo == o})
                }
                for o in orgs
            },
        },
        "families": {
            "retained_rows": dict(sorted(fam_ret.items())),
            "significant_rows": dict(sorted(fam_sig.items())),
            "per_organism_retained_rows": {
                o: {
                    fam: org_fam_ret[(o, fam)]
                    for fam in sorted({f for (oo, f) in org_fam_ret if oo == o})
                }
                for o in orgs
            },
            "per_organism_significant_rows": {
                o: {
                    fam: org_fam_sig[(o, fam)]
                    for fam in sorted({f for (oo, f) in org_fam_sig if oo == o})
                }
                for o in orgs
            },
        },
        "evidence": {
            "retained_rows_per_tier": {t: tier_ret[t] for t in TIERS},
            "significant_rows_per_tier": {t: tier_sig[t] for t in TIERS},
        },
        "labels": {
            "pathway_label_rows": n_labels,
            "inferred_pathway_rows": n_inferred,
        },
        "package": {
            "data_row_total": sum(rows_ret.values()) + n_labels,
            "payload_file_count": len(payload_files),
            "enrichment_file_count": sum(cells_per_org.values()),
        },
    }
    if m_family is not None:
        census["families"]["m_family_per_organism_per_family"] = m_family
    return census


def _m_family_per_organism_per_family(
    m_family_cell: dict[tuple, set],
) -> dict[str, dict[str, int]]:
    """Collapse per-cell ``m_family`` readings to ``{org: {family: M}}``.

    A cell must carry exactly one value, and the tiers of one (organism, family)
    share the family (same term2gene + universe), so they must agree too. Cells
    with no data rows contribute nothing.
    """
    per: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for (org, family, tier), values in sorted(m_family_cell.items()):
        if len(values) != 1:
            raise SystemExit(
                f"{org}/{family}_{tier}: m_family is not one value per cell: {sorted(values)}"
            )
        per[org][family] |= values
    out: dict[str, dict[str, int]] = {}
    for org in sorted(per):
        out[org] = {}
        for family in sorted(per[org]):
            values = per[org][family]
            if len(values) != 1:
                raise SystemExit(f"{org}/{family}: m_family differs across tiers: {sorted(values)}")
            out[org][family] = next(iter(values))
    return out


def _release_info(payload: Path) -> dict:
    info = payload / "RELEASE_INFO.json"
    if info.exists():
        try:
            doc = json.loads(info.read_text())
        except (OSError, ValueError):
            return {}
        return doc if isinstance(doc, dict) else {}
    return {}


def _contract_version(payload: Path) -> str | None:
    stamp = payload / "contract_version"
    if stamp.exists():
        return stamp.read_text(encoding="utf-8").strip() or None
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--payload", required=True, help="finalized release payload directory")
    ap.add_argument("--out", required=True, help="where to write release_census.json")
    args = ap.parse_args()
    census = build_census(Path(args.payload))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(census, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(f"census written -> {out}")
    for line in (
        f"organisms={census['organisms']}",
        f"cells={census['cells_total']}",
        f"retained={census['associations']['retained_total']:,}",
        f"significant={census['associations']['significant_total']:,}",
        f"mirnas_retained={census['mirnas']['retained_distinct_union']:,}",
        f"accessions={census['pathways']['distinct_accessions_union']:,}",
        f"records={census['pathways']['organism_resolved_records']:,}",
        f"labels={census['labels']['pathway_label_rows']:,} "
        f"(inferred {census['labels']['inferred_pathway_rows']:,})",
    ):
        print("  " + line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
