"""Tests for :mod:`mirpathdb_prod.sources.gene_identity` (W1 Excel-date guard).

Pins the identity rules shared with the webserver importer: date/datetime cells
and every stringified date/timestamp rendering are rejected; real month-like
symbols (MARCH1, SEPT2, DEC1) are accepted; Entrez cells normalize across the
openpyxl renderings; the ``ENTREZ:<id>`` sentinel is itself a valid symbol.
"""

from __future__ import annotations

import datetime

import pytest

from mirpathdb_prod.sources.gene_identity import (
    ENTREZ_SENTINEL_PREFIX,
    entrez_sentinel,
    is_date_artifact,
    is_valid_gene_symbol,
    normalize_entrez_id,
)


@pytest.mark.parametrize(
    "value",
    [
        "HIF1A", "SEPTIN2", "MARCH1", "SEPT2", "DEC1", "JUN", "MAY", "INF2",
        "HLA-DRB1", "C1orf43", " ACTB ", "ENTREZ:4735",
    ],
)
def test_valid_symbols_accepted(value: object) -> None:
    assert is_valid_gene_symbol(value)


@pytest.mark.parametrize(
    "value",
    [
        datetime.datetime(2022, 9, 2),   # SEPT2 autocorrect artifact
        datetime.datetime(2022, 3, 1),   # MARCH1
        datetime.date(2022, 12, 1),      # DEC1 (date, not datetime)
        "2022-09-02", "2022-09-02 00:00:00", "2022-09-02T00:00:00",
        "2022-3-1 0:00", "3/1/2022", "01-03-22", "2022/03/01", "1.3.2022",
        "1-Mar", "02-Sep", "2-Sep-22", "1 March 2022",
        "Mar-01", "March 1", "Sep-22", "September 2",
    ],
)
def test_date_artifacts_rejected(value: object) -> None:
    assert is_date_artifact(value)
    assert not is_valid_gene_symbol(value)


@pytest.mark.parametrize(
    "value",
    [
        None, "", "   ", "-", "--", "n/a", "NA", "nan", "None", "null",
        "unknown", "?", "0",
        4735, 4735.0, True, "4735", "4735.0", "1e5", "inf",
        "X" * 65,
    ],
)
def test_non_symbols_rejected(value: object) -> None:
    assert not is_valid_gene_symbol(value)


@pytest.mark.parametrize("value", ["MARCH1", "SEPT2", "DEC1", 4735, 4735.0, None, ""])
def test_non_dates_are_not_date_artifacts(value: object) -> None:
    assert not is_date_artifact(value)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (4735, "4735"),
        (4735.0, "4735"),
        ("4735", "4735"),
        ("4735.0", "4735"),
        ("  4735  ", "4735"),
        ("0004735", "4735"),
        (0, None),
        (-1, None),
        (4.5, None),
        (True, None),
        (None, None),
        ("", None),
        ("-", None),
        ("abc", None),
    ],
)
def test_normalize_entrez_id(value: object, expected: str | None) -> None:
    assert normalize_entrez_id(value) == expected


def test_entrez_sentinel_is_valid_symbol() -> None:
    s = entrez_sentinel("4735")
    assert s == "ENTREZ:4735"
    assert s.startswith(ENTREZ_SENTINEL_PREFIX)
    assert is_valid_gene_symbol(s)
