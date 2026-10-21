"""Production miRNA -> target-gene sets with the locked HIGH/MEDIUM/LOW rollup.

This is the **production** target-set builder for rat (rno, 10116) and macaque
(mml, 9544) -- the species that have NO miRPathDB2 ground truth, so we build the
real enrichment inputs straight from the four official sources (no v2
reconstruction). It turns raw source rows into per-confidence-band maps

    ``{tier: {miRNA_display_name: set(target_gene_SYMBOL)}}``    tier in HIGH/MEDIUM/LOW

following the LOCKED evidence rollup (see ``config/target_sources.yaml``):

* **Rat (rno):**
    - HIGH   = miRTarBase strong (``Functional MTI``) + DIANA verified-strong.
    - MEDIUM = miRTarBase weak (everything else) + DIANA CLIP/AGO-supported.
    - LOW    = TargetScan (conserved) + miRDB (predicted).
* **Macaque (mml):**  (no miRTarBase / miRDB coverage)
    - HIGH   = DIANA verified-strong.
    - MEDIUM = DIANA CLIP/AGO-supported.
    - LOW    = TargetScan (predicted).

ID space: every band is emitted in **gene SYMBOL** space, matching the GO
TERM2GENE built by :mod:`mirpathdb_prod.pathways.go_genesets`:

* DIANA ``gene_name``, miRTarBase ``Target Gene``, TargetScan ``Gene Symbol`` are
  ALREADY symbols -- used verbatim.
* miRDB targets are RefSeq mRNA accessions -- mapped to symbol **through the same
  ``org.*.eg.db``** the GO sets came from (:func:`refseq_to_symbol`) so the two
  spaces cannot drift. Accessions that do not resolve are dropped (recorded in the
  per-source census).

The bands are cumulative-by-evidence but emitted as DISJOINT-by-construction maps
the caller can union as it likes; HIGH ⊆ HIGH∪MEDIUM ⊆ HIGH∪MEDIUM∪LOW is the
intended evidence hierarchy. The DIANA HIGH/MEDIUM split is delegated to
:func:`mirpathdb_prod.sources.tarbase.rollup_tier`; the miRTarBase strong/weak
split reuses :data:`mirpathdb_prod.sources.mirtarbase.STRONG_SUPPORT_TYPE`.

Pure transform over already-loaded rows (no I/O here except the optional R-backed
RefSeq->symbol mapping, which shells out exactly like the GO materialiser).
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
import warnings
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from mirpathdb_prod.sources import mirtarbase, tarbase
from mirpathdb_prod.sources.gene_identity import is_date_artifact

# The user-facing confidence bands (the locked rollup).
TIERS: tuple[str, ...] = ("HIGH", "MEDIUM", "LOW")

# Default build-cache directory for the cached raw downloads (overridable).
# Tests import this; production callers thread the config value
# (paths.source_cache_dir) through ``cache_dir`` instead of relying on it.
DEFAULT_CACHE_DIR = "/tmp/mpd-data"

# Per-tier minimum miRNA floors, seeded from the production run manifests
# (provenance/production_enrichment_manifest.json). These guard against a silent
# zero-row / degenerate build: when a builder is run with ``allow_missing=False``
# a tier below its floor raises. The macaque MEDIUM floor is intentionally absent
# (macaque MEDIUM = 0 is CORRECT — its only validated MTI lands in HIGH) and the
# macaque HIGH floor is 1 (the lone marked validated MTI). These are the
# OWNER-CONFIRMABLE defaults; pass ``min_counts`` to override per run.
RAT_MIN_COUNTS: dict[str, int] = {
    # rat DIANA contributes >=1125 rows (26 miRNAs) and the four-source rollup
    # yields populous HIGH/MEDIUM/LOW tiers; floor each tier at >=1 miRNA here and
    # leave the heavier DIANA-row floor to the eager loader. (owner-confirmable)
    "HIGH": 1,
    "MEDIUM": 1,
    "LOW": 1,
}
MACAQUE_MIN_COUNTS: dict[str, int] = {
    # macaque HIGH = the single marked validated MTI (mml-miR-125b-5p->ANGPT2);
    # MEDIUM is exempt (legitimately 0); LOW = the TargetScan-predicted band
    # (293 miRNAs in the pinned run — floor a little below to allow drift).
    "HIGH": 1,
    # MEDIUM omitted on purpose -> exempt from the floor (0 is correct).
    "LOW": 293,
}

# Minimum DIANA rows expected for rat in the Other_species bundle (the pinned run
# carried ~1125 rat MTIs). Passed as ``min_diana_rows`` by the rat production
# runners into :func:`build_rat_target_sets`, which raises (under
# ``allow_missing=False``) if the rat DIANA slice is thinner than this -- so a
# silently-shrunk rat DIANA source cannot pass the tier floors (RAT_MIN_COUNTS
# only requires >=1 miRNA/tier, which TargetScan alone satisfies). (owner-confirmable)
RAT_DIANA_MIN_ROWS = 1125

# miRTarBase column accessors (already validated by its adapter).
_MTB_MIRNA = "miRNA"
_MTB_SUPPORT = "Support Type"
_MTB_SYMBOL = "Target Gene"

# DIANA column accessors.
_DIANA_MIRNA = "mirna_name"
_DIANA_SYMBOL = "gene_name"

# TargetScan column accessors.
_TS_MIRNA = "miRNA"
_TS_SYMBOL = "Gene Symbol"

# miRDB column accessors (the adapter's positional names).
_MIRDB_MIRNA = "mirna"
_MIRDB_TARGET = "target"

RSCRIPT_ENV = "MPD_RSCRIPT"


@dataclass
class TargetSetCensus:
    """Per-source provenance for one organism's production target-set build.

    Attributes:
        per_source: ``{source_id: {n_rows, n_mirnas, n_genes, ...}}``.
        per_tier: ``{tier: {n_mirnas, n_genes, n_pairs}}`` after rollup.
        notes: Free-text per-source notes (e.g. unmapped RefSeq count).
    """

    per_source: dict[str, dict[str, Any]] = field(default_factory=dict)
    per_tier: dict[str, dict[str, int]] = field(default_factory=dict)
    notes: dict[str, str] = field(default_factory=dict)


def _clean(value: Any) -> str | None:
    """Trim a cell to a non-empty string, else None.

    Excel-date guard (2026-08-24 audit, defect (b)): date/datetime-typed cells
    and stringified date/timestamp renderings (``"2022-09-02 00:00:00"``,
    ``"1-Mar"``, ...) are autocorrect artifacts, never identifiers -- they are
    dropped (None) instead of being stringified into the tier maps.
    :func:`mirpathdb_prod.sources.mirtarbase.load_mti` recovers such cells via
    the workbook's own Entrez column before they get here; this is the
    defense-in-depth backstop for every band-builder input path.
    """
    if value is None or is_date_artifact(value):
        return None
    text = str(value).strip()
    if not text or is_date_artifact(text):
        return None
    return text


def _add(target_map: dict[str, set[str]], mirna: str | None, gene: str | None) -> None:
    """Add one (miRNA, gene) pair to a tier map if both are present."""
    if mirna and gene:
        target_map.setdefault(mirna, set()).add(gene)


# --------------------------------------------------------------------------- #
# Per-source contributions (each returns {miRNA -> set(symbol)} for ONE band).
# --------------------------------------------------------------------------- #
def mirtarbase_band(
    mti_rows: Iterable[dict[str, Any]], *, strong_only: bool
) -> dict[str, set[str]]:
    """miRTarBase contribution to a band (symbol space).

    ``strong_only=True`` keeps only ``Functional MTI`` rows (HIGH); ``False``
    keeps the weak/other rows that are NOT strong (the MEDIUM-only remainder).
    """
    out: dict[str, set[str]] = {}
    for row in mti_rows:
        # Strip the Support Type cell before the strong/weak compare so a
        # whitespace/case variant (e.g. "Functional MTI ") is not silently demoted
        # HIGH->MEDIUM. mirtarbase.validate_vocabulary strips before comparing; this
        # path must agree with it.
        is_strong = _clean(row.get(_MTB_SUPPORT)) == mirtarbase.STRONG_SUPPORT_TYPE
        if strong_only and not is_strong:
            continue
        if not strong_only and is_strong:
            continue
        _add(out, _clean(row.get(_MTB_MIRNA)), _clean(row.get(_MTB_SYMBOL)))
    return out


def diana_bands(mti_rows: Iterable[dict[str, Any]]) -> dict[str, dict[str, set[str]]]:
    """DIANA contribution split into HIGH / MEDIUM via :func:`tarbase.rollup_tier`."""
    high: dict[str, set[str]] = {}
    medium: dict[str, set[str]] = {}
    for row in mti_rows:
        mirna = _clean(row.get(_DIANA_MIRNA))
        gene = _clean(row.get(_DIANA_SYMBOL))
        if tarbase.rollup_tier(row) == "HIGH":
            _add(high, mirna, gene)
        else:
            _add(medium, mirna, gene)
    return {"HIGH": high, "MEDIUM": medium}


def targetscan_band(pred_rows: Iterable[dict[str, Any]]) -> dict[str, set[str]]:
    """TargetScan contribution to LOW (symbol space, verbatim Gene Symbol)."""
    out: dict[str, set[str]] = {}
    for row in pred_rows:
        _add(out, _clean(row.get(_TS_MIRNA)), _clean(row.get(_TS_SYMBOL)))
    return out


# Soft ceiling on the miRDB RefSeq->SYMBOL unmapped fraction. Above this, the
# RefSeq->symbol join likely drifted (wrong org.db / stale accessions) and we
# WARN (warn-only -- the band still builds). Overridable per call.
MIRDB_UNMAPPED_WARN_FRACTION = 0.10


def mirdb_band(
    pred_rows: Iterable[dict[str, Any]],
    refseq_symbol: dict[str, str],
    *,
    unmapped_warn_fraction: float = MIRDB_UNMAPPED_WARN_FRACTION,
) -> tuple[dict[str, set[str]], int, int]:
    """miRDB contribution to LOW, mapping RefSeq targets to symbol.

    Returns ``(band, n_mapped_pairs, n_unmapped_pairs)``. Targets whose RefSeq
    accession has no symbol in ``refseq_symbol`` are dropped and counted.

    L4 (warn-only): if the unmapped fraction
    ``unmapped / (mapped + unmapped)`` exceeds ``unmapped_warn_fraction``
    (default :data:`MIRDB_UNMAPPED_WARN_FRACTION` = 0.10), ``warnings.warn`` so a
    drifted RefSeq->SYMBOL join (wrong org.db / stale accessions silently dropping
    targets) surfaces. Reuses the mapped/unmapped counters already tallied; does
    NOT change the band.
    """
    out: dict[str, set[str]] = {}
    mapped = 0
    unmapped = 0
    for row in pred_rows:
        mirna = _clean(row.get(_MIRDB_MIRNA))
        acc = _clean(row.get(_MIRDB_TARGET))
        if not mirna or not acc:
            continue
        # miRDB accessions are versionless (NM_xxxxx); strip any version suffix.
        base = acc.split(".", 1)[0]
        symbol = refseq_symbol.get(base) or refseq_symbol.get(acc)
        if symbol:
            _add(out, mirna, symbol)
            mapped += 1
        else:
            unmapped += 1
    total = mapped + unmapped
    if total and (unmapped / total) > unmapped_warn_fraction:
        warnings.warn(
            f"miRDB RefSeq->SYMBOL unmapped fraction "
            f"{unmapped}/{total} = {unmapped / total:.1%} exceeds "
            f"{unmapped_warn_fraction:.0%}: the org.db RefSeq join may have drifted "
            f"(wrong org.*.eg.db / stale accessions). Targets dropped silently "
            f"otherwise -- confirm the annotation package matches the organism.",
            stacklevel=2,
        )
    return out, mapped, unmapped


# --------------------------------------------------------------------------- #
# RefSeq -> SYMBOL mapping via org.*.eg.db (same package as the GO sets).
# --------------------------------------------------------------------------- #
_R_REFSEQ_DRIVER = r"""
suppressMessages(library(AnnotationDbi))
args     <- commandArgs(trailingOnly = TRUE)
pkg      <- args[[1]]
acc_file <- args[[2]]
out_file <- args[[3]]
suppressMessages(library(pkg, character.only = TRUE))
db <- get(pkg)
acc <- readLines(acc_file)
acc <- unique(trimws(acc))
acc <- acc[nzchar(acc)]
res <- suppressMessages(AnnotationDbi::select(
  db, keys = acc, columns = "SYMBOL", keytype = "REFSEQ"))
res <- res[!is.na(res$SYMBOL), ]
res <- res[!duplicated(res$REFSEQ), ]
write.table(res[, c("REFSEQ", "SYMBOL")], file = out_file, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
"""


def _resolve_rscript(rscript: str | None = None) -> str:
    if rscript:
        return rscript
    override = os.environ.get(RSCRIPT_ENV)
    if override:
        return override
    found = shutil.which("Rscript")
    if found:
        return found
    raise RuntimeError(
        "Rscript not found: set $MPD_RSCRIPT to the clusterProfiler conda env's Rscript"
    )


def refseq_to_symbol(
    accessions: Iterable[str], *, org_db_package: str, rscript: str | None = None
) -> dict[str, str]:
    """Map RefSeq mRNA accessions to gene SYMBOL via an ``org.*.eg.db`` package.

    Shells out to R once (the same annotation package the GO sets come from), so
    the miRDB symbol space cannot drift from the GO TERM2GENE space. Accessions
    are stripped of any version suffix before lookup. Returns a versionless-base
    accession -> symbol dict (accessions with no symbol are simply absent).

    Args:
        accessions: RefSeq mRNA accessions (versioned or not).
        org_db_package: The ``org.*.eg.db`` package name (e.g. ``org.Rn.eg.db``).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        ``{base_accession: symbol}``.

    Raises:
        RuntimeError: if the R subprocess fails.
    """
    bases = sorted({a.split(".", 1)[0].strip() for a in accessions if a and a.strip()})
    if not bases:
        return {}
    rscript_bin = _resolve_rscript(rscript)
    with tempfile.TemporaryDirectory(prefix="mpd_refseq_") as tmp:
        tmp_path = Path(tmp)
        acc_file = tmp_path / "accessions.txt"
        out_file = tmp_path / "refseq_symbol.tsv"
        driver_file = tmp_path / "driver_refseq.R"
        acc_file.write_text("\n".join(bases) + "\n", encoding="utf-8")
        driver_file.write_text(_R_REFSEQ_DRIVER, encoding="utf-8")
        proc = subprocess.run(
            [rscript_bin, str(driver_file), org_db_package, str(acc_file), str(out_file)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"RefSeq->symbol Rscript failed for {org_db_package} "
                f"(exit {proc.returncode}).\n--- stderr ---\n{proc.stderr}"
            )
        mapping: dict[str, str] = {}
        if out_file.exists():
            with out_file.open("r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    acc = (row.get("REFSEQ") or "").strip()
                    sym = (row.get("SYMBOL") or "").strip()
                    if acc and sym:
                        mapping[acc] = sym
        return mapping


# --------------------------------------------------------------------------- #
# Rollup: combine the per-source bands into HIGH/MEDIUM/LOW for one organism.
# --------------------------------------------------------------------------- #
def _merge_into(dst: dict[str, set[str]], src: dict[str, set[str]]) -> None:
    for mirna, genes in src.items():
        if genes:
            dst.setdefault(mirna, set()).update(genes)


def _summ(band: dict[str, set[str]]) -> dict[str, int]:
    genes: set[str] = set()
    pairs = 0
    for g in band.values():
        genes |= g
        pairs += len(g)
    return {"n_mirnas": len(band), "n_genes": len(genes), "n_pairs": pairs}


def _enforce_floors(
    organism: str,
    per_tier: dict[str, dict[str, int]],
    min_counts: dict[str, int],
    *,
    allow_missing: bool,
) -> None:
    """Raise unless every tier with a floor meets it (H10 fail-loud guard).

    A tier whose key is absent from ``min_counts`` is EXEMPT (e.g. macaque MEDIUM,
    which is legitimately 0). When ``allow_missing`` is True the check is skipped
    entirely (the documented degraded-source path). Reuses the per-tier miRNA
    counts already computed in :func:`_summ`.
    """
    if allow_missing:
        return
    for tier, floor in min_counts.items():
        got = per_tier.get(tier, {}).get("n_mirnas", 0)
        if got < floor:
            raise ValueError(
                f"[{organism}] tier {tier}: {got} miRNAs < floor {floor} "
                f"(expected non-empty production tier). Pass allow_missing=True or "
                f"override min_counts only if this thinness is real + documented."
            )


def build_rat_target_sets(
    *,
    mirtarbase_rows: Iterable[dict[str, Any]] | None = None,
    diana_rows: Iterable[dict[str, Any]] | None = None,
    targetscan_rows: Iterable[dict[str, Any]] | None = None,
    mirdb_rows: Iterable[dict[str, Any]] | None = None,
    refseq_symbol: dict[str, str] | None = None,
    min_counts: dict[str, int] | None = None,
    min_diana_rows: int | None = None,
    allow_missing: bool = False,
) -> tuple[dict[str, dict[str, set[str]]], TargetSetCensus]:
    """Build rat (rno) HIGH/MEDIUM/LOW target sets from the four sources.

    HIGH = miRTarBase strong + DIANA strong; MEDIUM = miRTarBase weak + DIANA
    CLIP; LOW = TargetScan + miRDB (RefSeq mapped to symbol via ``refseq_symbol``).

    Fail-loud (H10/M2): with ``allow_missing=False`` (the default), an expected
    source that is PASSED but contributes 0 rows raises ``ValueError`` rather than
    silently producing a thin/degenerate build. A source omitted entirely
    (``None``) is treated as "not requested" and skipped. Pass ``min_counts``
    (e.g. :data:`RAT_MIN_COUNTS`) to ALSO floor each tier's miRNA count -- the
    production runners do this; when ``min_counts`` is ``None`` no tier floor is
    enforced (so a caller supplying a single source still degrades gracefully).
    Set ``allow_missing=True`` to skip every fail-loud check (used by the
    source-degradation test).

    Returns ``({tier: {mirna: set(symbol)}}, census)``.
    """
    floors = min_counts
    high: dict[str, set[str]] = {}
    medium: dict[str, set[str]] = {}
    low: dict[str, set[str]] = {}
    census = TargetSetCensus()

    mtb_rows = list(mirtarbase_rows) if mirtarbase_rows is not None else None
    if mtb_rows is not None and not mtb_rows and not allow_missing:
        raise ValueError(
            "[rno] miRTarBase: 0 rows for a requested source "
            "(expected non-empty). Pass allow_missing=True for a documented gap."
        )
    if mtb_rows:
        mtb_strong = mirtarbase_band(mtb_rows, strong_only=True)
        mtb_weak = mirtarbase_band(mtb_rows, strong_only=False)
        _merge_into(high, mtb_strong)
        _merge_into(medium, mtb_weak)
        census.per_source["mirtarbase"] = {
            "n_rows": len(mtb_rows),
            "n_mirnas_strong": len(mtb_strong),
            "n_mirnas_weak": len(mtb_weak),
        }

    diana_list = list(diana_rows) if diana_rows is not None else None
    if diana_list is not None and not diana_list and not allow_missing:
        raise ValueError(
            "[rno] DIANA-TarBase: 0 rows for a requested source "
            "(expected non-empty). Pass allow_missing=True for a documented gap."
        )
    if (
        min_diana_rows is not None
        and diana_list is not None
        and not allow_missing
        and len(diana_list) < min_diana_rows
    ):
        raise ValueError(
            f"[rno] DIANA-TarBase: {len(diana_list)} rows < expected floor "
            f"{min_diana_rows} (RAT_DIANA_MIN_ROWS). A thinned DIANA slice would "
            f"silently under-populate the rat HIGH/MEDIUM tiers while TargetScan "
            f"alone still satisfies RAT_MIN_COUNTS. Pass allow_missing=True or lower "
            f"min_diana_rows only if this thinness is real + documented."
        )
    if diana_list:
        bands = diana_bands(diana_list)
        _merge_into(high, bands["HIGH"])
        _merge_into(medium, bands["MEDIUM"])
        census.per_source["diana_tarbase"] = {
            "n_rows": len(diana_list),
            "n_mirnas_high": len(bands["HIGH"]),
            "n_mirnas_medium": len(bands["MEDIUM"]),
            "tier_rule": tarbase.DIANA_TIER_RULE,
        }

    ts_rows = list(targetscan_rows) if targetscan_rows is not None else None
    if ts_rows is not None and not ts_rows and not allow_missing:
        raise ValueError(
            "[rno] TargetScan: 0 rows for a requested source "
            "(expected non-empty). Pass allow_missing=True for a documented gap."
        )
    if ts_rows:
        ts_band = targetscan_band(ts_rows)
        _merge_into(low, ts_band)
        census.per_source["targetscan"] = {"n_rows": len(ts_rows), "n_mirnas": len(ts_band)}

    mdb_rows = list(mirdb_rows) if mirdb_rows is not None else None
    if mdb_rows is not None and not mdb_rows and not allow_missing:
        raise ValueError(
            "[rno] miRDB: 0 rows for a requested source "
            "(expected non-empty). Pass allow_missing=True for a documented gap."
        )
    if mdb_rows:
        mdb_band, mapped, unmapped = mirdb_band(mdb_rows, refseq_symbol or {})
        _merge_into(low, mdb_band)
        census.per_source["mirdb"] = {
            "n_rows": len(mdb_rows),
            "n_mirnas": len(mdb_band),
            "n_mapped_pairs": mapped,
            "n_unmapped_pairs": unmapped,
        }
        if unmapped:
            census.notes["mirdb"] = (
                f"{unmapped} RefSeq targets dropped (no SYMBOL in org.db); {mapped} mapped"
            )

    tiers = {"HIGH": high, "MEDIUM": medium, "LOW": low}
    census.per_tier = {t: _summ(b) for t, b in tiers.items()}
    if floors is not None:
        _enforce_floors("rno", census.per_tier, floors, allow_missing=allow_missing)
    return tiers, census


def build_macaque_target_sets(
    *,
    diana_rows: Iterable[dict[str, Any]] | None = None,
    targetscan_rows: Iterable[dict[str, Any]] | None = None,
    min_counts: dict[str, int] | None = None,
    allow_missing: bool = False,
) -> tuple[dict[str, dict[str, set[str]]], TargetSetCensus]:
    """Build macaque (mml) HIGH/MEDIUM/LOW target sets.

    Macaque has NO miRTarBase / miRDB coverage: HIGH/MEDIUM come ONLY from DIANA
    (verified -> HIGH, CLIP -> MEDIUM) and LOW from TargetScan. In DIANA v9 the
    macaque validated set is a single MTI, so HIGH is intentionally tiny -- that is
    the real extent of macaque's experimentally-validated evidence; the bulk of
    macaque enrichment input is the TargetScan-predicted LOW band.

    Fail-loud (H10/M2): with ``allow_missing=False`` (the default), a source that
    is PASSED but contributes 0 rows raises ``ValueError``. The macaque-specific
    truth -- MEDIUM = 0 and HIGH = 1 -- is CORRECT, so :data:`MACAQUE_MIN_COUNTS`
    floors HIGH at 1 and LOW at 293 but EXEMPTS MEDIUM (no MEDIUM key); pass it as
    ``min_counts`` (the runners do). When ``min_counts`` is ``None`` no tier floor
    is enforced. ``allow_missing=True`` skips every fail-loud check.

    Returns ``({tier: {mirna: set(symbol)}}, census)``.
    """
    floors = min_counts
    high: dict[str, set[str]] = {}
    medium: dict[str, set[str]] = {}
    low: dict[str, set[str]] = {}
    census = TargetSetCensus()

    diana_list = list(diana_rows) if diana_rows is not None else None
    if diana_list is not None and not diana_list and not allow_missing:
        raise ValueError(
            "[mml] DIANA-TarBase: 0 rows for a requested source. Macaque has a "
            "single validated MTI in the Other_species bundle -- 0 rows means the "
            "slice failed (wrong species string / reshaped file), not a real gap. "
            "Pass allow_missing=True only if this is genuinely documented."
        )
    if diana_list:
        bands = diana_bands(diana_list)
        _merge_into(high, bands["HIGH"])
        _merge_into(medium, bands["MEDIUM"])
        census.per_source["diana_tarbase"] = {
            "n_rows": len(diana_list),
            "n_mirnas_high": len(bands["HIGH"]),
            "n_mirnas_medium": len(bands["MEDIUM"]),
            "tier_rule": tarbase.DIANA_TIER_RULE,
        }

    ts_rows = list(targetscan_rows) if targetscan_rows is not None else None
    if ts_rows is not None and not ts_rows and not allow_missing:
        raise ValueError(
            "[mml] TargetScan: 0 rows for a requested source (macaque's "
            "substantive LOW band -- 293 predicted miRNAs -- cannot be empty). "
            "Pass allow_missing=True for a documented gap."
        )
    if ts_rows:
        ts_band = targetscan_band(ts_rows)
        _merge_into(low, ts_band)
        census.per_source["targetscan"] = {"n_rows": len(ts_rows), "n_mirnas": len(ts_band)}

    tiers = {"HIGH": high, "MEDIUM": medium, "LOW": low}
    census.per_tier = {t: _summ(b) for t, b in tiers.items()}
    if floors is not None:
        _enforce_floors("mml", census.per_tier, floors, allow_missing=allow_missing)
    return tiers, census


def build_symbol_casemap(go_symbols: Iterable[str]) -> dict[str, str]:
    """Build a casefold -> canonical-symbol map from the GO TERM2GENE symbols.

    The org.*.eg.db GO space uses the organism's native symbol casing
    (rat: ``Aadat``, ``Actb``; macaque/human: ``AADAT``). TargetScan ships
    human-anchored UPPERCASE symbols for EVERY organism, so for rat its symbols
    (``AADAT``) would miss the GO space (``Aadat``) entirely -- the classic
    symbol-case silent failure. This map lets :func:`normalize_to_go_symbols`
    re-case a band's symbols onto the GO space. If two GO symbols collide under
    casefold (rare), the first in iteration order wins; pass a ``sorted`` iterable
    for determinism.
    """
    out: dict[str, str] = {}
    for sym in go_symbols:
        s = sym.strip()
        if s:
            out.setdefault(s.casefold(), s)
    return out


def normalize_to_go_symbols(
    band: dict[str, set[str]], casemap: dict[str, str]
) -> tuple[dict[str, set[str]], int, int]:
    """Re-case a band's gene symbols onto the GO symbol space via ``casemap``.

    A symbol already in the GO space (its casefold maps back to itself) is kept;
    a symbol that differs only by case (TargetScan ``AADAT`` vs GO ``Aadat``) is
    rewritten to the canonical GO casing; a symbol with no casefold match in the
    GO space is DROPPED (it cannot be enriched against the GO universe anyway).

    Returns ``(normalized_band, n_kept_pairs, n_dropped_pairs)``.
    """
    out: dict[str, set[str]] = {}
    kept = 0
    dropped = 0
    for mirna, genes in band.items():
        for gene in genes:
            canon = casemap.get(gene.casefold())
            if canon:
                out.setdefault(mirna, set()).add(canon)
                kept += 1
            else:
                dropped += 1
    return out, kept, dropped


def normalize_tiers_to_go_symbols(
    tiers: dict[str, dict[str, set[str]]], go_symbols: Iterable[str]
) -> tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, int]]]:
    """Re-case ALL tiers' symbols onto the GO space (see :func:`normalize_to_go_symbols`).

    Returns ``(normalized_tiers, {tier: {n_kept, n_dropped}})``. Apply this after
    building tiers and before running the ORA so every band lives in the GO
    TERM2GENE symbol space (the case fix for TargetScan-derived LOW symbols).
    """
    casemap = build_symbol_casemap(sorted(go_symbols))
    out: dict[str, dict[str, set[str]]] = {}
    stats: dict[str, dict[str, int]] = {}
    for tier, band in tiers.items():
        norm, kept, dropped = normalize_to_go_symbols(band, casemap)
        out[tier] = norm
        stats[tier] = {"n_kept_pairs": kept, "n_dropped_pairs": dropped}
    return out, stats


# --------------------------------------------------------------------------- #
# Production materialisation for the all-four-source organisms (hsa, mmu, rno).
# --------------------------------------------------------------------------- #
# Canonical taxid -> cached miRTarBase per-species workbook filename.
_MIRTARBASE_FILE: dict[int, str] = {
    9606: "mirtarbase_hsa_MTI.xlsx",
    10090: "mirtarbase_mmu_MTI.xlsx",
    10116: "mirtarbase_rno_MTI.xlsx",
}

# Canonical taxid -> cached DIANA-TarBase v9 per-species TSV filename. (rat +
# macaque live in the Other_species bundle and are handled by the rno/mml runner;
# this materialiser covers the per-species-file organisms hsa + mmu.)
_DIANA_FILE: dict[int, str] = {
    9606: "Homo_sapiens_TarBase-v9.tsv.gz",
    10090: "Mus_musculus_TarBase-v9.tsv.gz",
}

# Cached TargetScan + miRDB are single multi-species downloads sliced by taxid.
_TARGETSCAN_FILE = "Predicted_Targets_Context_Scores.default_predictions.txt.zip"
_MIRDB_FILE = "miRDB_v6.0_prediction_result.txt.gz"


def load_and_build_target_sets(
    taxid: int,
    *,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    org_db_package: str | None = None,
    rscript: str | None = None,
) -> tuple[dict[str, dict[str, set[str]]], TargetSetCensus]:
    """Load the four cached sources for an all-four-source organism and roll up.

    This is the production target-set materialiser for the organisms that ship a
    miRTarBase per-species workbook AND a DIANA per-species TSV -- **human (9606)
    and mouse (10090)** (rat 10116 is also wired here but its production run uses
    the Other_species DIANA bundle via the rno/mml runner). It locates each cached
    download under ``cache_dir``, loads it with the matching adapter, maps miRDB
    RefSeq accessions to SYMBOL through ``org_db_package`` (the same annotation the
    GO/Reactome/WP sets use, so the symbol space cannot drift), and rolls the four
    sources up into HIGH/MEDIUM/LOW via :func:`build_rat_target_sets` (the locked
    all-four-source rollup: miRTarBase strong + DIANA strong -> HIGH; miRTarBase
    weak + DIANA CLIP -> MEDIUM; TargetScan + miRDB -> LOW).

    Args:
        taxid: Canonical NCBI taxid (one of 9606 / 10090 / 10116).
        cache_dir: Build-cache dir holding the cached raw downloads.
        org_db_package: ``org.*.eg.db`` for the miRDB RefSeq->SYMBOL join; if
            ``None`` it is looked up from ``taxid`` via the GO source mapping.
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        ``({tier: {mirna: set(symbol)}}, census)`` -- raw (non-recased) symbol
        space; apply :func:`normalize_tiers_to_go_symbols` before ORA.

    Raises:
        ValueError: if ``taxid`` is not an all-four-source organism.
        FileNotFoundError: if a required cached download is absent.
    """
    # Local import keeps the module import cheap (mirdb/targetscan pull stdlib zip
    # + gzip only, but keeping the import here also avoids a cycle with the GO
    # source when org_db_package must be resolved).
    from mirpathdb_prod.sources import mirdb, targetscan

    if taxid not in _MIRTARBASE_FILE or taxid not in _DIANA_FILE:
        raise ValueError(
            f"load_and_build_target_sets is for the all-four-source organisms "
            f"{sorted(set(_MIRTARBASE_FILE) & set(_DIANA_FILE))}; "
            f"got taxid {taxid}. rat/macaque use the Other_species runner."
        )

    cache = Path(cache_dir)

    def _cached(name: str) -> Path:
        p = cache / name
        if not p.exists():
            raise FileNotFoundError(f"cached source not found: {p}")
        return p

    if org_db_package is None:
        from mirpathdb_prod.pathways.go_genesets import org_db_for_taxid

        org_db_package = org_db_for_taxid(taxid)

    diana_rows = tarbase.load_mti(_cached(_DIANA_FILE[taxid]), taxid=taxid)
    mtb_rows = mirtarbase.load_mti(_cached(_MIRTARBASE_FILE[taxid]), taxid=taxid)
    ts_rows = targetscan.load_predictions(_cached(_TARGETSCAN_FILE), taxid=taxid)
    mdb_rows = mirdb.load_predictions(_cached(_MIRDB_FILE), taxid=taxid)

    accs = {r["target"] for r in mdb_rows}
    refseq_symbol = refseq_to_symbol(accs, org_db_package=org_db_package, rscript=rscript)

    return build_rat_target_sets(
        mirtarbase_rows=mtb_rows,
        diana_rows=diana_rows,
        targetscan_rows=ts_rows,
        mirdb_rows=mdb_rows,
        refseq_symbol=refseq_symbol,
    )


def cumulative_tier(
    tiers: dict[str, dict[str, set[str]]], up_to: str
) -> dict[str, set[str]]:
    """Union the bands HIGH..``up_to`` into one ``{miRNA: set(symbol)}`` map.

    ``cumulative_tier(t, "HIGH")`` is HIGH only; ``"MEDIUM"`` is HIGH∪MEDIUM;
    ``"LOW"`` is HIGH∪MEDIUM∪LOW -- the evidence hierarchy used as ORA input.
    """
    if up_to not in TIERS:
        raise ValueError(f"unknown tier {up_to!r}; expected one of {TIERS}")
    order = TIERS[: TIERS.index(up_to) + 1]
    out: dict[str, set[str]] = {}
    for tier in order:
        _merge_into(out, tiers.get(tier, {}))
    return out
