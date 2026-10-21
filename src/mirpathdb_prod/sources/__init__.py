"""Source adapters.

Adapters in this subpackage ONLY fetch and load raw data from upstream or
reference locations. They do no normalization and produce no pass/fail
judgments — that belongs to :mod:`mirpathdb_prod.normalize` and
:mod:`mirpathdb_prod.validation` respectively. No prototype code is imported.
"""

from __future__ import annotations

import csv
import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import Any

# The canonical organism-table schema shared by every organism TSV in this repo
# (validation/organisms_old.tsv, validation/organisms_mirpathdb2.tsv, and the
# generated organisms_new.tsv). Adapters validate the header against this so a
# silently-reshaped file fails loudly instead of producing garbage records.
CANONICAL_ORGANISM_COLUMNS: tuple[str, ...] = (
    "canonical_tax_id",
    "scientific_name",
    "source_name",
    "common_name",
    "source",
    "source_version",
    "source_record_id",
    "normalization_notes",
)


def read_organism_tsv(path: Path) -> list[dict[str, Any]]:
    """Read a canonical-schema organism TSV into a list of raw row dicts.

    Pure load: returns the rows verbatim (no normalization, no judgment). The
    header is checked against :data:`CANONICAL_ORGANISM_COLUMNS`; a mismatch
    raises ``ValueError`` so a malformed reference file is caught at load time.

    Args:
        path: Path to a tab-separated organism table with the canonical header.

    Returns:
        One ``dict`` per data row, keyed by column name (values are strings).

    Raises:
        FileNotFoundError: if ``path`` does not exist.
        ValueError: if the header does not match the canonical schema.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"organism TSV not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header != CANONICAL_ORGANISM_COLUMNS:
            raise ValueError(
                f"{path}: unexpected organism-table header\n"
                f"  expected: {CANONICAL_ORGANISM_COLUMNS}\n"
                f"  got:      {header}"
            )
        return [dict(row) for row in reader]


def validate_species_column(
    rows: Iterable[dict[str, Any]],
    *,
    column: str,
    taxid: int,
    path: Path | str,
    source_label: str,
    mismatch_error_fraction: float = 0.5,
) -> None:
    """Check rows' own species annotations against the requested ``taxid``.

    Inputs keyed to an organism must agree with the data's own species column,
    not just the filename convention: a per-species file passed with the wrong
    ``taxid`` would otherwise load silently. The expected scientific name comes
    from the curated taxonomy (:mod:`mirpathdb_prod.normalize.taxonomy`); names
    are compared after :func:`~mirpathdb_prod.normalize.taxonomy.normalize_name`
    so casing/whitespace variants do not false-alarm. Rows with a blank/absent
    species value are skipped (nothing to compare).

    Behavior:
        * ``taxid`` outside the curated taxonomy -> no check (nothing to compare
          against).
        * mismatch fraction > ``mismatch_error_fraction`` (default 0.5) ->
          ``ValueError`` naming the file, the expected species, and the observed
          values (wrong-file detection).
        * a minority of rows mismatching -> ``warnings.warn`` with a count; rows
          are NOT filtered (cross-species entries stay visible downstream).
    """
    from mirpathdb_prod.normalize.taxonomy import canonical_scientific_name, normalize_name

    expected = canonical_scientific_name(taxid)
    if expected is None:
        return
    expected_key = normalize_name(expected)

    n_checked = 0
    n_mismatch = 0
    observed: set[str] = set()
    for row in rows:
        raw = row.get(column)
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        n_checked += 1
        if normalize_name(text) != expected_key:
            n_mismatch += 1
            observed.add(text)

    if not n_checked or not n_mismatch:
        return
    if n_mismatch / n_checked > mismatch_error_fraction:
        raise ValueError(
            f"{path}: {source_label} {column!r} column contradicts taxid {taxid}: "
            f"expected {expected!r} but {n_mismatch}/{n_checked} rows carry "
            f"{sorted(observed)} -- this looks like the wrong species file for "
            f"the requested organism."
        )
    warnings.warn(
        f"{path}: {source_label}: {n_mismatch}/{n_checked} rows' {column!r} value "
        f"differs from the expected {expected!r} (taxid {taxid}); observed "
        f"{sorted(observed)}. Rows are kept (not filtered).",
        stacklevel=3,
    )
