#!/usr/bin/env python3
"""Regenerate ``data/reference/wikipathways_name_repairs.tsv`` from GPML titles.

The upstream WikiPathways GMT export mangles HTML entities in pathway titles
(DATA-001; see ``docs/DATA-001_wikipathways_label_root_cause.md``). This script
rebuilds the checked-in repair map deterministically for the PINNED release:

1. scan the cached (or freshly downloaded) per-organism GMTs for display names
   that trip the residue gates in
   :mod:`mirpathdb_prod.pathways.label_normalization`;
2. download the SAME release's per-organism GPML archives from
   ``data.wikipathways.org/<release>/gpml/`` (cached beside the GMTs);
3. take each flagged WPID's canonical title from its GPML ``Name`` attribute
   (XML parsing decodes entities correctly);
4. write the TSV with full provenance (release, corrupt name, canonical name,
   GPML member filename, GPML zip sha256).

Run inside the production environment::

    python scripts/build_wp_name_repairs.py --cache-dir /tmp/mpd-data \
        --out data/reference/wikipathways_name_repairs.tsv

The map is release-pinned reference data and is committed; a rebuild of the
release does not need the network as long as the map covers every flagged id
(the gates fail the build otherwise).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mirpathdb_prod.pathways.label_normalization import (  # noqa: E402
    label_violations,
    normalize_label,
)
from mirpathdb_prod.pathways.wikipathways_genesets import (  # noqa: E402
    WP_DATA_BASE,
    WP_RELEASE,
    WP_SPECIES_FILE,
    fetch_gmt,
    parse_wp_gmt,
)

_WPID_IN_MEMBER = re.compile(r"_(WP\d+)_")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch_gpml_zip(species: str, release: str, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"wikipathways-{release}-gpml-{species}.zip"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = f"{WP_DATA_BASE}/{release}/gpml/{dest.name}"
    with urllib.request.urlopen(url, timeout=600) as resp:  # noqa: S310 (https only)
        dest.write_bytes(resp.read())
    return dest


def _gpml_titles(zip_path: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    """Return ``{wpid: (canonical title, member filename)}`` for ``wanted``."""
    out: dict[str, tuple[str, str]] = {}
    with zipfile.ZipFile(zip_path) as zf:
        for member in zf.namelist():
            m = _WPID_IN_MEMBER.search(member)
            if not m or m.group(1) not in wanted:
                continue
            root = ET.fromstring(zf.read(member))
            name = root.get("Name")
            if name:
                out[m.group(1)] = (normalize_label(name), member)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", required=True, help="source cache dir (GMTs + GPML zips)")
    ap.add_argument("--release", default=WP_RELEASE)
    ap.add_argument("--out", default="data/reference/wikipathways_name_repairs.tsv")
    args = ap.parse_args()

    rows: list[dict[str, str]] = []
    for taxid in sorted(WP_SPECIES_FILE):
        species = WP_SPECIES_FILE[taxid]
        gmt = fetch_gmt(taxid, cache_dir=args.cache_dir, release=args.release)
        _term2entrez, term2name = parse_wp_gmt(gmt)
        flagged = {
            wpid: normalize_label(name)
            for wpid, name in term2name.items()
            if label_violations(normalize_label(name))
        }
        if not flagged:
            continue
        gpml_zip = _fetch_gpml_zip(species, args.release, Path(args.cache_dir) / "wikipathways")
        zip_sha = _sha256(gpml_zip)
        titles = _gpml_titles(gpml_zip, set(flagged))
        for wpid, corrupt in sorted(flagged.items()):
            if wpid not in titles:
                print(f"ERROR: no GPML title for flagged {wpid} ({species})", file=sys.stderr)
                return 1
            canonical, member = titles[wpid]
            if label_violations(canonical):
                print(f"ERROR: GPML title for {wpid} still flagged: {canonical!r}", file=sys.stderr)
                return 1
            rows.append(
                {
                    "release": args.release,
                    "pathway_id": wpid,
                    "corrupt_name": corrupt,
                    "display_name": canonical,
                    "source": f"gpml:{member}",
                    "source_sha256": zip_sha,
                }
            )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["release", "pathway_id", "corrupt_name", "display_name", "source", "source_sha256"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} repair row(s) -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
