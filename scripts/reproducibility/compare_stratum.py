#!/usr/bin/env python3
"""Compare a recomputed enrichment stratum with the shipped one, keyed on miRNA and pathway.

A stratum is one ``enrichment_significant/<org>/<FAMILY>_<TIER>.tsv`` of a release package
(contract version 2: ``mirna, pathway_id, p_value, p_adjust, overlap, p_adjust_overlap,
m_overlap, m_family``). The recomputed file reproduces the shipped one when

* both carry the contract header, parse cleanly and hold no duplicate (mirna, pathway_id) key;
* the key sets are identical (and so are the row counts);
* the integer columns (``overlap``, ``m_overlap``, ``m_family``) are exactly equal;
* every P column (``p_value``, ``p_adjust``, ``p_adjust_overlap``) agrees to a relative 1e-9,
  ``|recomputed - shipped| <= 1e-9 * |shipped|`` (a shipped zero needs an exact zero).

Byte identity and row order are reported, not required (exact equality of the text failed in
the 2026-09-15 migration check; relative 1e-9 passed with a largest difference of 5.6e-12).

With ``--package`` and ``--relpath`` the shipped file is read from the package and first checked
against the package's ``SHA256SUMS``; with ``--sums-copy`` that ``SHA256SUMS`` is itself checked
against an independent copy (for 2026-09-15 the one tracked in this repository under
``docs/reproducibility/checksums/``), so a swapped or edited file is never parsed.

Usage::

    python scripts/reproducibility/compare_stratum.py \\
        --rerun data/processed/enrichment_significant/rno/GO_BP_HIGH.tsv \\
        --package <mirpathdb-data-2026-09-15> --relpath enrichment_significant/rno/GO_BP_HIGH.tsv \\
        --sums-copy docs/reproducibility/checksums/mirpathdb-data-2026-09-15/SHA256SUMS \\
        --expected-rows 13560 --out comparison.json

Standard library only. Exit status: 0 the strata agree, 1 they do not, 2 an input is missing
or fails its checksum.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

REL_TOL = 1e-9
HEADER = (
    "mirna",
    "pathway_id",
    "p_value",
    "p_adjust",
    "overlap",
    "p_adjust_overlap",
    "m_overlap",
    "m_family",
)
KEY_COLUMNS = ("mirna", "pathway_id")
INTEGER_COLUMNS = ("overlap", "m_overlap", "m_family")
P_COLUMNS = ("p_value", "p_adjust", "p_adjust_overlap")
MAX_EXAMPLES = 5

_SUMS_LINE = re.compile(r"^([0-9a-f]{64}) [ *](.+)$")


class IntegrityError(RuntimeError):
    """An input is missing or does not match its recorded checksum."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_sha256sums(path: Path) -> dict[str, str]:
    """Parse ``sha256sum`` output into ``{relative path without './': digest}``."""
    sums: dict[str, str] = {}
    for n, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        m = _SUMS_LINE.match(line)
        if not m:
            raise IntegrityError(f"{path}:{n}: not a sha256sum line: {line!r}")
        sums[m.group(2).removeprefix("./")] = m.group(1)
    return sums


def read_release_manifest_rows(path: Path) -> dict[str, int]:
    """``RELEASE_MANIFEST.tsv`` (path, bytes, sha256, data_rows, ...) -> ``{path: data_rows}``."""
    rows: dict[str, int] = {}
    with Path(path).open(encoding="utf-8", newline="") as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if rec.get("data_rows", "").strip().isdigit():
                rows[rec["path"].removeprefix("./")] = int(rec["data_rows"])
    return rows


@dataclass
class Cell:
    """One stratum file, parsed: header, keys in file order, the row of each key."""

    path: Path
    header: tuple[str, ...]
    keys: list[tuple[str, str]] = field(default_factory=list)
    rows: dict[tuple[str, str], list[str]] = field(default_factory=dict)
    duplicate_keys: list[tuple[str, str]] = field(default_factory=list)
    malformed_lines: list[int] = field(default_factory=list)

    @property
    def n_rows(self) -> int:
        return len(self.keys)


def load_cell(path: Path) -> Cell:
    """Read a stratum TSV (tab-separated, no quoting). A repeated key keeps its first row."""
    path = Path(path)
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        header = tuple(next(reader, ()))
        cell = Cell(path=path, header=header)
        if not all(c in header for c in KEY_COLUMNS):
            return cell
        k0, k1 = (header.index(c) for c in KEY_COLUMNS)
        for lineno, rec in enumerate(reader, start=2):
            if len(rec) != len(header):
                cell.malformed_lines.append(lineno)
                continue
            key = (rec[k0], rec[k1])
            if key in cell.rows:
                cell.duplicate_keys.append(key)
            else:
                cell.rows[key] = rec
            cell.keys.append(key)
    return cell


def relative_difference(recomputed: float, shipped: float) -> float:
    """``|recomputed - shipped| / |shipped|``; 0 when equal, infinite when shipped is 0 or NaN."""
    if recomputed == shipped:
        return 0.0
    if math.isnan(recomputed) or math.isnan(shipped) or shipped == 0.0:
        return math.inf
    return abs(recomputed - shipped) / abs(shipped)


def _compare_integers(rerun: Cell, shipped: Cell, common, col: str) -> dict:
    i_r, i_s = rerun.header.index(col), shipped.header.index(col)
    bad = []
    for key in common:
        a, b = rerun.rows[key][i_r], shipped.rows[key][i_s]
        if a == b:
            continue
        try:
            same = int(a) == int(b)
        except ValueError:
            same = False
        if not same:
            bad.append((key, a, b))
    return {
        "n_mismatch": len(bad),
        "examples": [
            {"key": list(k), "rerun": a, "shipped": b} for k, a, b in bad[:MAX_EXAMPLES]
        ],
    }


def _compare_p(rerun: Cell, shipped: Cell, common, col: str, rel_tol: float) -> dict:
    i_r, i_s = rerun.header.index(col), shipped.header.index(col)
    n_text_identical = 0
    max_rel = 0.0
    over = []
    for key in common:
        a, b = rerun.rows[key][i_r], shipped.rows[key][i_s]
        if a == b:
            n_text_identical += 1
            continue
        try:
            d = relative_difference(float(a), float(b))
        except ValueError:
            d = math.inf
        max_rel = max(max_rel, d)
        if not d <= rel_tol:
            over.append((key, a, b, d))
    return {
        "n_text_identical": n_text_identical,
        "max_relative_difference": max_rel if math.isfinite(max_rel) else "inf",
        "n_over_tolerance": len(over),
        "examples": [
            {
                "key": list(k),
                "rerun": a,
                "shipped": b,
                "relative_difference": d if math.isfinite(d) else "inf",
            }
            for k, a, b, d in over[:MAX_EXAMPLES]
        ],
    }


def compare_cells(
    rerun_path: Path,
    shipped_path: Path,
    *,
    rel_tol: float = REL_TOL,
    expected_rows: int | None = None,
) -> dict:
    """Compare two stratum files; the result's ``pass`` is the verdict."""
    rerun, shipped = load_cell(rerun_path), load_cell(shipped_path)
    rerun_sha, shipped_sha = sha256_file(rerun_path), sha256_file(shipped_path)
    res: dict = {
        "rerun_path": str(rerun_path),
        "shipped_path": str(shipped_path),
        "rerun_sha256": rerun_sha,
        "shipped_sha256": shipped_sha,
        "byte_identical": rerun_sha == shipped_sha,
        "header_rerun": list(rerun.header),
        "header_shipped": list(shipped.header),
        "header_matches": rerun.header == shipped.header == HEADER,
        "n_rows_rerun": rerun.n_rows,
        "n_rows_shipped": shipped.n_rows,
        "malformed_lines_rerun": rerun.malformed_lines[:MAX_EXAMPLES],
        "malformed_lines_shipped": shipped.malformed_lines[:MAX_EXAMPLES],
        "n_malformed_rerun": len(rerun.malformed_lines),
        "n_malformed_shipped": len(shipped.malformed_lines),
        "duplicate_keys_rerun": len(rerun.duplicate_keys),
        "duplicate_keys_shipped": len(shipped.duplicate_keys),
        "row_order_identical": rerun.keys == shipped.keys,
        "expected_rows": expected_rows,
    }
    only_r = sorted(set(rerun.rows) - set(shipped.rows))
    only_s = sorted(set(shipped.rows) - set(rerun.rows))
    res["keys_only_in_rerun"] = len(only_r)
    res["keys_only_in_shipped"] = len(only_s)
    res["keys_only_in_rerun_examples"] = [list(k) for k in only_r[:MAX_EXAMPLES]]
    res["keys_only_in_shipped_examples"] = [list(k) for k in only_s[:MAX_EXAMPLES]]

    columns_present = all(
        c in rerun.header and c in shipped.header for c in INTEGER_COLUMNS + P_COLUMNS
    )
    common = [k for k in shipped.rows if k in rerun.rows]
    res["n_keys_compared"] = len(common) if columns_present else 0
    res["integers"] = (
        {c: _compare_integers(rerun, shipped, common, c) for c in INTEGER_COLUMNS}
        if columns_present
        else {}
    )
    res["p_values"] = (
        {c: _compare_p(rerun, shipped, common, c, rel_tol) for c in P_COLUMNS}
        if columns_present
        else {}
    )
    res["rel_tol"] = rel_tol
    res["expected_rows_met"] = expected_rows is None or (
        rerun.n_rows == expected_rows and shipped.n_rows == expected_rows
    )
    res["pass"] = bool(
        res["header_matches"]
        and columns_present
        and not rerun.malformed_lines
        and not shipped.malformed_lines
        and not rerun.duplicate_keys
        and not shipped.duplicate_keys
        and rerun.n_rows == shipped.n_rows
        and not only_r
        and not only_s
        and all(v["n_mismatch"] == 0 for v in res["integers"].values())
        and all(v["n_over_tolerance"] == 0 for v in res["p_values"].values())
        and res["expected_rows_met"]
    )
    return res


def check_shipped_file(package: Path, relpath: str, sums_copy: Path | None = None) -> dict:
    """Verify a package file against the package's SHA256SUMS (and that against a copy).

    Raises :class:`IntegrityError` on any mismatch; returns what was checked.
    """
    package = Path(package)
    sums_path = package / "SHA256SUMS"
    if not sums_path.is_file():
        raise IntegrityError(f"{sums_path} not found")
    out: dict = {"package": str(package), "relpath": relpath}
    if sums_copy is not None:
        copy_sha, pkg_sha = sha256_file(sums_copy), sha256_file(sums_path)
        out["sha256sums_digest"] = pkg_sha
        out["sha256sums_copy"] = str(sums_copy)
        if copy_sha != pkg_sha:
            raise IntegrityError(
                f"{sums_path} (sha256 {pkg_sha}) differs from the copy {sums_copy} "
                f"(sha256 {copy_sha})"
            )
    sums = read_sha256sums(sums_path)
    if relpath not in sums:
        raise IntegrityError(f"{relpath} is not listed in {sums_path}")
    shipped = package / relpath
    if not shipped.is_file():
        raise IntegrityError(f"{shipped} not found")
    digest = sha256_file(shipped)
    if digest != sums[relpath]:
        raise IntegrityError(
            f"{shipped}: sha256 {digest} does not match SHA256SUMS ({sums[relpath]})"
        )
    out["shipped_sha256"] = digest
    out["shipped_sha256_matches_SHA256SUMS"] = True
    manifest = package / "RELEASE_MANIFEST.tsv"
    if manifest.is_file():
        out["release_manifest_rows"] = read_release_manifest_rows(manifest).get(relpath)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--rerun", type=Path, required=True, help="recomputed stratum TSV")
    ap.add_argument("--shipped", type=Path, help="shipped stratum TSV (or --package/--relpath)")
    ap.add_argument("--package", type=Path, help="release package directory")
    ap.add_argument("--relpath", help="the stratum's path inside the package")
    ap.add_argument("--sums-copy", type=Path, help="independent copy of the package SHA256SUMS")
    ap.add_argument("--expected-rows", type=int, help="row count both files must have")
    ap.add_argument("--rel-tol", type=float, default=REL_TOL, help="relative P tolerance")
    ap.add_argument("--out", type=Path, help="write the comparison as JSON here")
    args = ap.parse_args(argv)

    integrity: dict = {}
    try:
        if args.package is not None:
            if not args.relpath:
                ap.error("--package needs --relpath")
            integrity = check_shipped_file(args.package, args.relpath, args.sums_copy)
            shipped = args.package / args.relpath
        elif args.shipped is not None:
            shipped = args.shipped
        else:
            ap.error("give --shipped or --package with --relpath")
        if not args.rerun.is_file():
            raise IntegrityError(f"{args.rerun} not found")
        if not shipped.is_file():
            raise IntegrityError(f"{shipped} not found")
    except IntegrityError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    res = compare_cells(
        args.rerun, shipped, rel_tol=args.rel_tol, expected_rows=args.expected_rows
    )
    if integrity:
        res["integrity"] = integrity
        rows_manifest = integrity.get("release_manifest_rows")
        if rows_manifest is not None and rows_manifest != res["n_rows_shipped"]:
            res["pass"] = False
            res["release_manifest_rows_mismatch"] = rows_manifest
    payload = json.dumps(res, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload, encoding="utf-8")
    p_max = {c: v["max_relative_difference"] for c, v in res["p_values"].items()}
    print(
        f"rows {res['n_rows_rerun']}/{res['n_rows_shipped']}  "
        f"keys only in rerun/shipped {res['keys_only_in_rerun']}/{res['keys_only_in_shipped']}  "
        f"integer mismatches {sum(v['n_mismatch'] for v in res['integers'].values())}  "
        f"max relative P difference {p_max}  byte-identical {res['byte_identical']}  "
        f"{'PASS' if res['pass'] else 'FAIL'}"
    )
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
