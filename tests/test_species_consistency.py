"""Row-level species-consistency validation for the per-species MTI loaders.

Inputs keyed to an organism must agree with the data's own species annotations,
not just the filename/sheet convention: a wrong file passed with a plausible name
must fail loudly. Covered here for DIANA-TarBase (``species`` column) and
miRTarBase (``Species (miRNA)`` column) with synthetic fixtures:

* the correct per-species file loads without raising or warning;
* a swapped-species file (majority mismatch) raises with the file, the expected
  species, and the observed values in the message;
* a minority mismatch warns with a count and keeps every row (no filtering).
"""

from __future__ import annotations

import csv
import gzip
import warnings
from pathlib import Path

import pytest

from mirpathdb_prod.sources import mirtarbase, tarbase

try:
    import openpyxl

    _HAS_OPENPYXL = True
except ImportError:  # pragma: no cover - openpyxl is a validation extra
    _HAS_OPENPYXL = False


# --------------------------------------------------------------------------- #
# DIANA-TarBase per-species files (species column).
# --------------------------------------------------------------------------- #
def _tarbase_row(species: str, mirna: str, gene: str) -> dict[str, str]:
    return {
        "species": species,
        "mirna_name": mirna,
        "gene_name": gene,
        "experimental_method": "Luciferase Reporter Assay",
        "confidence": "1",
    }


def _write_tarbase(path: Path, rows: list[dict[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tarbase.TARBASE_COLUMNS), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in tarbase.TARBASE_COLUMNS})


def test_tarbase_correct_species_file_passes(tmp_path: Path) -> None:
    p = tmp_path / "Homo_sapiens_TarBase-v9.tsv.gz"
    _write_tarbase(
        p,
        [
            _tarbase_row("Homo sapiens", "hsa-let-7a-5p", "HMGA2"),
            _tarbase_row("Homo sapiens", "hsa-miR-21-5p", "PTEN"),
        ],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any species/vocab warning would fail here
        rows = tarbase.load_mti(p, taxid=9606)
    assert len(rows) == 2


def test_tarbase_swapped_species_file_raises(tmp_path: Path) -> None:
    # A mouse file loaded as human (majority mismatch) is wrong-file detection.
    p = tmp_path / "Homo_sapiens_TarBase-v9.tsv.gz"
    _write_tarbase(
        p,
        [
            _tarbase_row("Mus musculus", "mmu-let-7a-5p", "Hmga2"),
            _tarbase_row("Mus musculus", "mmu-miR-21a-5p", "Pten"),
        ],
    )
    with pytest.raises(ValueError) as excinfo:
        tarbase.load_mti(p, taxid=9606)
    msg = str(excinfo.value)
    assert str(p) in msg
    assert "Homo sapiens" in msg  # expected species named
    assert "Mus musculus" in msg  # observed values named


def test_tarbase_minority_mismatch_warns_and_keeps_rows(tmp_path: Path) -> None:
    p = tmp_path / "Homo_sapiens_TarBase-v9.tsv.gz"
    _write_tarbase(
        p,
        [
            _tarbase_row("Homo sapiens", "hsa-let-7a-5p", "HMGA2"),
            _tarbase_row("Homo sapiens", "hsa-miR-21-5p", "PTEN"),
            _tarbase_row("Homo sapiens", "hsa-miR-155-5p", "SOCS1"),
            _tarbase_row("Mus musculus", "mmu-miR-1a-3p", "Hand2"),  # minority
        ],
    )
    with pytest.warns(UserWarning, match=r"1/4 rows"):
        rows = tarbase.load_mti(p, taxid=9606)
    assert len(rows) == 4  # mismatching rows are kept, not filtered


def test_tarbase_case_and_whitespace_variants_do_not_false_alarm(tmp_path: Path) -> None:
    # Comparison runs through normalize_name, so casing/underscores are tolerated.
    p = tmp_path / "Homo_sapiens_TarBase-v9.tsv.gz"
    _write_tarbase(p, [_tarbase_row("homo_sapiens", "hsa-let-7a-5p", "HMGA2")])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rows = tarbase.load_mti(p, taxid=9606)
    assert len(rows) == 1


# --------------------------------------------------------------------------- #
# miRTarBase per-species workbooks (Species (miRNA) column).
# --------------------------------------------------------------------------- #
def _mtb_row(species: str, mirna: str, gene: str) -> tuple:
    return (
        "MIRT000000",
        mirna,
        species,
        gene,
        1234,
        species,
        "Luciferase reporter assay",
        "Functional MTI",
        12345678,
    )


def _write_mti_xlsx(path: Path, sheet: str, rows: list[tuple]) -> None:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = sheet
    ws.append(list(mirtarbase.MTI_COLUMNS))
    for r in rows:
        ws.append(list(r))
    wb.save(path)


@pytest.mark.skipif(not _HAS_OPENPYXL, reason="openpyxl (validation extra) not installed")
def test_mirtarbase_correct_species_workbook_passes(tmp_path: Path) -> None:
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        "Homo sapiens",
        [
            _mtb_row("Homo sapiens", "hsa-let-7a-5p", "HMGA2"),
            _mtb_row("Homo sapiens", "hsa-miR-21-5p", "PTEN"),
        ],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rows = mirtarbase.load_mti(p, taxid=9606)
    assert len(rows) == 2


@pytest.mark.skipif(not _HAS_OPENPYXL, reason="openpyxl (validation extra) not installed")
def test_mirtarbase_swapped_species_workbook_raises(tmp_path: Path) -> None:
    # A mouse workbook loaded as human: the single-sheet fallback would accept the
    # sheet regardless of name, so the row-level check is what catches it.
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        "Mus musculus",
        [
            _mtb_row("Mus musculus", "mmu-let-7a-5p", "Hmga2"),
            _mtb_row("Mus musculus", "mmu-miR-21a-5p", "Pten"),
        ],
    )
    with pytest.raises(ValueError) as excinfo:
        mirtarbase.load_mti(p, taxid=9606)
    msg = str(excinfo.value)
    assert str(p) in msg
    assert "Homo sapiens" in msg
    assert "Mus musculus" in msg


@pytest.mark.skipif(not _HAS_OPENPYXL, reason="openpyxl (validation extra) not installed")
def test_mirtarbase_minority_mismatch_warns_and_keeps_rows(tmp_path: Path) -> None:
    # e.g. a viral-prefix MTI inside the human workbook: warn + keep, never filter.
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        "Homo sapiens",
        [
            _mtb_row("Homo sapiens", "hsa-let-7a-5p", "HMGA2"),
            _mtb_row("Homo sapiens", "hsa-miR-21-5p", "PTEN"),
            _mtb_row("Homo sapiens", "hsa-miR-155-5p", "SOCS1"),
            _mtb_row("Epstein-Barr virus", "ebv-miR-BART1-5p", "CASP3"),  # minority
        ],
    )
    with pytest.warns(UserWarning, match=r"1/4 rows"):
        rows = mirtarbase.load_mti(p, taxid=9606)
    assert len(rows) == 4


@pytest.mark.skipif(not _HAS_OPENPYXL, reason="openpyxl (validation extra) not installed")
def test_mirtarbase_blank_species_cells_are_skipped(tmp_path: Path) -> None:
    # Rows without a species annotation carry nothing to compare; no false alarm.
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        "Homo sapiens",
        [
            _mtb_row("Homo sapiens", "hsa-let-7a-5p", "HMGA2"),
            _mtb_row("", "hsa-miR-21-5p", "PTEN"),  # blank species cell
        ],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        rows = mirtarbase.load_mti(p, taxid=9606)
    assert len(rows) == 2
