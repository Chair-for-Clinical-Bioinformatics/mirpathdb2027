"""Regression tests for the Excel-date Target Gene defect (2026-08-24 audit, (b)).

The miRTarBase MTI workbooks carry ``Target Gene`` cells Excel autocorrected
into date objects (``SEPT2`` -> ``datetime(2022, 9, 2)``, ``MARCH1`` ->
``datetime(2022, 3, 1)``); loading them verbatim shipped stringified timestamps
as gene symbols (750 rows in the target-pairs-2026-08-16 deliverable). These
tests pin the corrected semantics, mirroring the webserver importer:

* date/datetime-typed cells and stringified date/timestamp forms never survive
  :func:`mirpathdb_prod.sources.mirtarbase.load_mti`;
* recovery uses the workbook's own ``Target Gene (Entrez ID)`` column, resolved
  to a symbol only when that Entrez ID maps to exactly ONE valid symbol
  elsewhere in the same workbook;
* otherwise the ``ENTREZ:<id>`` sentinel is stored;
* rows with neither a valid symbol nor a usable Entrez ID are skipped and
  counted (per-species summary printed);
* no emitted symbol ever matches a date / timestamp / numeric pattern;
* ``production_target_sets._clean`` is a defense-in-depth backstop dropping
  date artifacts on every band-builder input path.
"""

from __future__ import annotations

import datetime
import re
from pathlib import Path

import pytest

from mirpathdb_prod.sources import mirtarbase
from mirpathdb_prod.sources import production_target_sets as pts
from mirpathdb_prod.sources.gene_identity import is_valid_gene_symbol

openpyxl = pytest.importorskip(
    "openpyxl", reason="openpyxl (validation extra) not installed"
)


def _write_mti_xlsx(path: Path, rows: list[tuple]) -> None:
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Homo sapiens"
    ws.append(list(mirtarbase.MTI_COLUMNS))
    for r in rows:
        ws.append(list(r))
    wb.save(path)


def _row(
    mirt: str,
    mirna: object,
    gene: object,
    entrez: object,
    support: str = "Functional MTI",
    pmid: int = 1,
) -> tuple:
    return (
        mirt, mirna, "Homo sapiens", gene, entrez, "Homo sapiens",
        "qPCR", support, pmid,
    )


def test_date_cell_recovered_to_workbook_symbol(tmp_path: Path, capsys) -> None:
    """datetime(2022, 9, 2) + Entrez 4735 -> the symbol the workbook uses (SEPTIN2)."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", datetime.datetime(2022, 9, 2), 4735),
            _row("MIRT2", "hsa-miR-2-3p", "SEPTIN2", 4735,
                 support="Functional MTI (Weak)"),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert len(rows) == 2
    assert [r["Target Gene"] for r in rows] == ["SEPTIN2", "SEPTIN2"]
    out = capsys.readouterr().out
    assert "rows read 2" in out
    assert "date cells 1 (recovered via Entrez 1)" in out
    assert "recovered to symbol 1" in out
    assert "skipped 0" in out


def test_date_cell_without_workbook_symbol_gets_sentinel(tmp_path: Path) -> None:
    """No valid symbol for the Entrez ID anywhere in the workbook -> ENTREZ:<id>."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", datetime.datetime(2022, 3, 1), 55016),
            _row("MIRT2", "hsa-miR-2-3p", "ACTB", 60),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert {r["Target Gene"] for r in rows} == {"ENTREZ:55016", "ACTB"}


def test_ambiguous_entrez_resolution_gets_sentinel(tmp_path: Path) -> None:
    """An Entrez ID used with TWO distinct valid symbols is ambiguous -> sentinel."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", datetime.datetime(2022, 9, 2), 9999),
            _row("MIRT2", "hsa-miR-2-3p", "AAA1", 9999,
                 support="Functional MTI (Weak)"),
            _row("MIRT3", "hsa-miR-3-3p", "BBB2", 9999,
                 support="Functional MTI (Weak)"),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert sorted(r["Target Gene"] for r in rows) == ["AAA1", "BBB2", "ENTREZ:9999"]


def test_stringified_date_and_timestamp_forms_recovered(tmp_path: Path) -> None:
    """String renderings of the artifact are as forbidden as the typed cells."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", "2022-09-02 00:00:00", 4735),
            _row("MIRT2", "hsa-miR-2-3p", "1-Mar", 51765,
                 support="Functional MTI (Weak)"),
            _row("MIRT3", "hsa-miR-3-3p", "SEPTIN2", 4735,
                 support="Functional MTI (Weak)"),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert [r["Target Gene"] for r in rows] == ["SEPTIN2", "ENTREZ:51765", "SEPTIN2"]


def test_unrecoverable_rows_are_skipped_and_counted(tmp_path: Path, capsys) -> None:
    """Neither a valid symbol nor a usable Entrez ID -> counted skip, never a timestamp."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", datetime.datetime(2022, 12, 1), None),
            _row("MIRT2", "hsa-miR-2-3p", "3/1/2022", "-"),
            _row("MIRT3", None, "TP53", 7157),  # empty miRNA -> skipped too
            _row("MIRT4", "hsa-miR-4-3p", "TP53", 7157),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert [r["Target Gene"] for r in rows] == ["TP53"]
    out = capsys.readouterr().out
    assert "skipped 3" in out
    assert "invalid_gene_no_entrez" in out
    assert "empty_mirna" in out


def test_month_like_symbols_survive_verbatim(tmp_path: Path) -> None:
    """Plain-string MARCH1/SEPT2/DEC1/JUN/MAY are real symbols, never touched."""
    genes = ["MARCH1", "SEPT2", "DEC1", "JUN", "MAY"]
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [_row(f"MIRT{i}", f"hsa-miR-{i}-3p", g, 1000 + i) for i, g in enumerate(genes)],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    assert [r["Target Gene"] for r in rows] == genes


def test_no_output_symbol_matches_date_timestamp_or_numeric_patterns(
    tmp_path: Path,
) -> None:
    """Invariant: every emitted Target Gene is a valid symbol-shaped string."""
    p = tmp_path / "hsa_MTI.xlsx"
    _write_mti_xlsx(
        p,
        [
            _row("MIRT1", "hsa-miR-1-3p", datetime.datetime(2022, 9, 2), 4735),
            _row("MIRT2", "hsa-miR-2-3p", "SEPTIN2", 4735,
                 support="Functional MTI (Weak)"),
            _row("MIRT3", "hsa-miR-3-3p", datetime.datetime(2022, 3, 1), 55016),
            _row("MIRT4", "hsa-miR-4-3p", "MARCH1", 55016,
                 support="Functional MTI (Weak)"),
            _row("MIRT5", "hsa-miR-5-3p", "2022-12-01 00:00:00", 50514,
                 support="Functional MTI (Weak)"),
            _row("MIRT6", "hsa-miR-6-3p", "4735", 4735,
                 support="Functional MTI (Weak)"),  # numeric rendering of the ID
            _row("MIRT7", "hsa-miR-7-3p", datetime.datetime(2022, 12, 1), None),
        ],
    )
    rows = mirtarbase.load_mti(p, taxid=9606)
    date_or_time = re.compile(r"\d{4}-\d{1,2}-\d{1,2}|\d{1,2}:\d{2}")
    for r in rows:
        g = r["Target Gene"]
        assert isinstance(g, str), g
        assert not isinstance(g, (datetime.date, datetime.datetime))
        assert is_valid_gene_symbol(g), g
        assert not date_or_time.search(g), g
        assert not re.fullmatch(r"\d+(\.\d+)?", g), g
    # MIRT1/MIRT6 recover to SEPTIN2 (4735), MIRT3 to MARCH1 (55016),
    # MIRT5 to the sentinel (50514 has no symbol), MIRT7 is skipped.
    assert [r["Target Gene"] for r in rows] == [
        "SEPTIN2", "SEPTIN2", "MARCH1", "MARCH1", "ENTREZ:50514", "SEPTIN2",
    ]


def test_clean_rejects_date_artifacts() -> None:
    """production_target_sets._clean drops typed AND stringified date artifacts."""
    assert pts._clean(datetime.datetime(2022, 9, 2)) is None
    assert pts._clean(datetime.date(2022, 3, 1)) is None
    assert pts._clean("2022-09-02 00:00:00") is None
    assert pts._clean("2022-09-02") is None
    assert pts._clean("1-Mar") is None
    assert pts._clean("Mar-01") is None
    assert pts._clean(None) is None
    assert pts._clean("   ") is None
    assert pts._clean(" SEPTIN2 ") == "SEPTIN2"
    assert pts._clean("MARCH1") == "MARCH1"
    assert pts._clean("ENTREZ:4735") == "ENTREZ:4735"
    assert pts._clean("Functional MTI") == "Functional MTI"


def test_mirtarbase_band_drops_raw_date_cells() -> None:
    """Rows that bypass load_mti still cannot put a timestamp into a band."""
    rows = [
        {"miRNA": "hsa-miR-1", "Target Gene": datetime.datetime(2022, 9, 2),
         "Support Type": "Functional MTI"},
        {"miRNA": "hsa-miR-1", "Target Gene": "ACTB",
         "Support Type": "Functional MTI"},
        {"miRNA": "hsa-miR-2", "Target Gene": "2022-03-01 00:00:00",
         "Support Type": "Functional MTI (Weak)"},
    ]
    strong = pts.mirtarbase_band(rows, strong_only=True)
    weak = pts.mirtarbase_band(rows, strong_only=False)
    assert strong == {"hsa-miR-1": {"ACTB"}}
    assert weak == {}
