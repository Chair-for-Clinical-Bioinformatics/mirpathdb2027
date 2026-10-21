#!/usr/bin/env python3
"""Derive the label-corrected immutable release package (DATA-001 / WS7).

The enrichment contract of the source package is UNCHANGED (byte-asserted); this
derivation replaces only the display-label crosswalk (regenerated through the fixed
normalization pipeline), adds the release census + provenance records, and
regenerates every manifest under a NEW package identity — corrected bytes are never
published under an existing checksum.

Fail-closed: everything is assembled under <out>.building and renamed at the end.

    python scripts/derive_label_corrected_package.py \
        --source .../mirpathdb-data-2026-08-30 \
        --labels .../rebuild-labels-2026-08-31/pathway_labels.tsv \
        --out    .../mirpathdb-data-2026-08-31 \
        --release-id mirpathdb-contract-v1.2-2026-08-31 \
        --authoritative-build docs/AUTHORITATIVE_BUILD.md
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mirpathdb_prod.pathways.label_normalization import (  # noqa: E402
    assert_labels_clean,
    label_violations,
    load_name_repairs,
    normalize_label,
)

NEW_ADDONS = ("LABEL_CORRECTIONS.tsv", "release_census.json", "AUTHORITATIVE_BUILD.md")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fail(msg: str) -> None:
    print(f"FATAL: {msg}", file=sys.stderr)
    raise SystemExit(1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--labels", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--release-id", required=True)
    ap.add_argument("--authoritative-build", type=Path, required=True)
    args = ap.parse_args()

    src, out = args.source, args.out
    building = out.with_suffix(".building")
    if out.exists():
        fail(f"{out} already exists — a package is immutable, never rebuilt in place")
    if building.exists():
        shutil.rmtree(building)

    # ---- 1. copy source, preserving relative symlinks --------------------------------------
    shutil.copytree(src, building, symlinks=True)

    # ---- 2. replace the label crosswalk (gate + row-count + diff record) --------------------
    old_labels = {}
    with (src / "pathway_labels.tsv").open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            old_labels[r["pathway_id"]] = r["display_name"]
    new_rows = []
    with args.labels.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            new_rows.append(r)
    if len(new_rows) != len(old_labels):
        fail(f"label row count changed: {len(old_labels)} -> {len(new_rows)}")
    assert_labels_clean(((r["pathway_id"], r["display_name"]) for r in new_rows),
                        context="derived package labels")
    repairs = load_name_repairs(release="20260610")
    corrections = []
    for r in new_rows:
        old = old_labels[r["pathway_id"]]
        if old != r["display_name"]:
            reason = ("upstream_gmt_entity_mangling_repaired_from_gpml"
                      if r["pathway_id"] in repairs else "whitespace_normalization")
            corrections.append((r["pathway_id"], old, r["display_name"], reason))
    if not corrections:
        fail("no label corrections — nothing to derive")
    for pid, old, new, _ in corrections:
        if normalize_label(old) != new and pid not in repairs:
            fail(f"{pid}: correction is not the pipeline normalization of the old label")
        if label_violations(new):
            fail(f"{pid}: corrected label still flagged")
    shutil.copy2(args.labels, building / "pathway_labels.tsv")
    with (building / "LABEL_CORRECTIONS.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pathway_id", "previous_display_name", "display_name", "reason"])
        w.writerows(corrections)
    print(f"labels: {len(new_rows)} rows, {len(corrections)} corrected "
          f"({sum(1 for c in corrections if c[3].startswith('upstream'))} repairs, "
          f"{sum(1 for c in corrections if c[3] == 'whitespace_normalization')} whitespace)")

    # ---- 3. byte-assert the untouched contract ----------------------------------------------
    cells = sorted((src / "enrichment_significant").glob("*/*.tsv"))
    for cell in cells:
        rel = cell.relative_to(src)
        if sha256(cell) != sha256(building / rel):
            fail(f"cell changed during derivation: {rel}")
    if sha256(src / "contract_version") != sha256(building / "contract_version"):
        fail("contract_version changed")
    print(f"cells byte-identical: {len(cells)}/57 expected 57" if len(cells) == 57
          else fail(f"cell count {len(cells)} != 57"))

    # ---- 4. provenance records ---------------------------------------------------------------
    shutil.copy2(args.authoritative_build, building / "AUTHORITATIVE_BUILD.md")
    info = json.loads((src / "RELEASE_INFO.json").read_text(encoding="utf-8"))
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    git_head = subprocess.run(
        ["git", "-C", str(Path(__file__).resolve().parents[1]), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=True).stdout.strip()
    info.update({
        "release_id": args.release_id,
        "generated_utc": now,
        "prior_release": info["release_id"],
        "derived_from": {
            "package": src.name,
            "sha256sums_digest": sha256(src / "SHA256SUMS"),
            "note": ("label-only derivation: the 57 enrichment cells and contract_version "
                     "are byte-identical to the prior package (asserted at build time); "
                     "only pathway_labels.tsv changed, plus added provenance records."),
        },
        "label_corrections": {
            "rows_changed": len(corrections),
            "apostrophe_repairs_from_gpml": sum(
                1 for c in corrections if c[3].startswith("upstream")),
            "whitespace_normalizations": sum(
                1 for c in corrections if c[3] == "whitespace_normalization"),
            "root_cause": ("upstream WikiPathways 20260610 GMT export mangles HTML entities "
                           "in pathway titles; see AUTHORITATIVE_BUILD.md and the backend's "
                           "docs/DATA-001_wikipathways_label_root_cause.md"),
            "details_file": "LABEL_CORRECTIONS.tsv",
            "fix_commit": git_head,
        },
    })
    (building / "RELEASE_INFO.json").write_text(
        json.dumps(info, indent=2) + "\n", encoding="utf-8")
    readme = building / "README.md"
    readme.write_text(readme.read_text(encoding="utf-8") + f"""
## {args.release_id} — display-label correction (DATA-001)

Derived from `{info['prior_release']}`: the 57 enrichment files and every
association are byte-identical; only `pathway_labels.tsv` changed
({len(corrections)} display names — 7 upstream-mangled WikiPathways apostrophe
titles repaired from the release's own GPML, 42 upstream Reactome double-space
names collapsed; full diff in `LABEL_CORRECTIONS.tsv`). `release_census.json`
is the machine-readable census every served statistic must equal
(`AUTHORITATIVE_BUILD.md` records the census decision).
""", encoding="utf-8")

    # ---- 5+6 ORDER (fix): the census records the digest of the FINAL contract manifest, so
    # the 59-file RELEASE_MANIFEST/SHA256SUMS MUST be regenerated BEFORE the census runs, and
    # the census MUST exist before ADDON_SHA256SUMS/UPLOAD_MANIFEST hash it. Previously the
    # census ran here (step 5) against the copied-in PRIOR SHA256SUMS, so _meta.sha256sums_digest
    # named the superseded package. Topology is non-circular: SHA256SUMS covers only the 59
    # scientific files (never the census); the census references SHA256SUMS' digest; the add-on
    # and upload manifests cover the census.
    # ---- 6a. contract manifest first (RELEASE_MANIFEST.tsv + SHA256SUMS: the 59-file public contract).
    old_man = {}
    with (src / "RELEASE_MANIFEST.tsv").open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            old_man[parts[0]] = parts
    man_lines = [header]
    sums_lines = []
    payload_files = ["contract_version", "pathway_labels.tsv"] + sorted(
        p.relative_to(building).as_posix()
        for p in (building / "enrichment_significant").glob("*/*.tsv"))
    for rel in payload_files:
        p = building / rel
        h, b = sha256(p), p.stat().st_size
        old = old_man[rel]
        if rel == "pathway_labels.tsv":
            dr = str(len(new_rows))
        else:
            dr = old[3]
            if (h, str(b)) != (old[2], old[1]):
                fail(f"unexpected byte change in contract file {rel}")
        man_lines.append(f"{rel}\t{b}\t{h}\t{dr}\t{old[4]}\t{old[5]}\t{old[6]}")
        sums_lines.append(f"{h}  ./{rel}")
    if len(sums_lines) != 59:
        fail(f"SHA256SUMS has {len(sums_lines)} entries, expected 59")
    (building / "RELEASE_MANIFEST.tsv").write_text("\n".join(man_lines) + "\n")
    (building / "SHA256SUMS").write_text("\n".join(sums_lines) + "\n")

    # ---- 6b. census AFTER the final SHA256SUMS exists, so it binds to the deployed manifest.
    subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent / "build_release_census.py"),
         "--payload", str(building), "--out", str(building / "release_census.json")],
        check=True)
    census_doc = json.loads((building / "release_census.json").read_text(encoding="utf-8"))
    if census_doc["_meta"]["sha256sums_digest"] != sha256(building / "SHA256SUMS"):
        fail("census sha256sums_digest does not bind to the final SHA256SUMS (ordering defect)")

    # ---- 6c. ADDON_SHA256SUMS: the source addon set (paths preserved) + the new records.
    addon_paths = []
    with (src / "ADDON_SHA256SUMS").open(encoding="utf-8") as fh:
        for line in fh:
            addon_paths.append(line.split("  ", 1)[1].strip().lstrip("./"))
    addon_paths += list(NEW_ADDONS)
    addon_paths = sorted(set(addon_paths))
    (building / "ADDON_SHA256SUMS").write_text("".join(
        f"{sha256(building / rel)}  ./{rel}\n" for rel in addon_paths), encoding="utf-8")

    # UPLOAD_MANIFEST.sha256: every real file except itself.
    upload = sorted(p for p in building.rglob("*")
                    if p.is_file() and not p.is_symlink()
                    and p.name != "UPLOAD_MANIFEST.sha256")
    (building / "UPLOAD_MANIFEST.sha256").write_text("".join(
        f"{sha256(p)}  ./{p.relative_to(building).as_posix()}\n" for p in upload),
        encoding="utf-8")

    # ---- 7. verify + seal --------------------------------------------------------------------
    for manifest in ("SHA256SUMS", "ADDON_SHA256SUMS", "UPLOAD_MANIFEST.sha256"):
        res = subprocess.run(["sha256sum", "-c", "--quiet", manifest], cwd=building,
                             capture_output=True, text=True)
        if res.returncode != 0:
            fail(f"{manifest} verification failed:\n{res.stderr[:500]}{res.stdout[:500]}")
        print(f"{manifest}: verified")
    building.rename(out)
    print(f"\nSEALED: {out}")
    print(f"  release_id: {args.release_id}")
    print(f"  SHA256SUMS digest:       {sha256(out / 'SHA256SUMS')}")
    print(f"  RELEASE_MANIFEST digest: {sha256(out / 'RELEASE_MANIFEST.tsv')}")
    print(f"  upload files: {len(upload)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
