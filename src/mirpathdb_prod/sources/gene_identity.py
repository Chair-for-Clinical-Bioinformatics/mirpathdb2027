"""Gene-symbol identity validation for Excel-derived target sources (W1).

Why this module exists
----------------------
The upstream miRTarBase MTI workbooks ship with ``Target Gene`` cells that Excel
has autocorrected into **date objects** -- the classic MARCH1/SEPT2/DEC1
gene-symbol corruption (``SEPT2`` -> ``datetime(2022, 9, 2)``, ``MARCH1`` ->
``datetime(2022, 3, 1)``). Verified in the cached production workbooks:

* ``mirtarbase_hsa_MTI.xlsx`` -- 1,863 date-typed Target Gene cells (24 distinct dates)
* ``mirtarbase_mmu_MTI.xlsx`` --   397 date-typed cells
* ``mirtarbase_rno_MTI.xlsx`` --    18 date-typed cells

Every one of those cells carries a valid ``Target Gene (Entrez ID)`` value, which
is the authoritative recovery path. A naive ``str(cell)`` turns the datetime into
``"2022-09-02 00:00:00"`` and propagates it as a gene "symbol" -- exactly what the
backend build did before this module existed (2026-08-24 audit, defect (b): 750
date-corrupted ``target_gene`` rows in the target-pairs-2026-08-16 deliverable).

This module is the single backend authority on whether a value may be used as a
gene symbol. :func:`mirpathdb_prod.sources.mirtarbase.load_mti` consults it for
every ``Target Gene`` cell and routes failures through Entrez-based recovery or
an explicit, counted skip -- never a silent drop and never a silently stored
timestamp. It deliberately mirrors the webserver's importer-side twin
(``webserver/mirpathdb/gene_identity.py``) so the generator and the importer
accept the SAME identity space; keep the two in sync.

Rules enforced by :func:`is_valid_gene_symbol`
----------------------------------------------
Rejected:

* ``datetime``/``date`` objects (the raw autocorrect artifact)
* stringified dates and timestamps (``"2022-03-01 00:00:00"``, ``"2022-09-02"``,
  ``"3/1/2022"``, ``"1-Mar"``, ``"Mar-01"``, ``"2-Sep-22"`` ...)
* empty / whitespace-only values and obvious metadata placeholders
  (``"nan"``, ``"N/A"``, ``"-"``, ``"none"``, ``"null"``, ``"unknown"``, ``"?"``)
* purely numeric values, including float renderings (``"4735"``, ``"4735.0"``, ``"1e5"``)
* values longer than 64 characters (the webserver ``Gene.symbol`` column width)

Accepted: everything else, notably real symbols that merely *look* month-like as
plain strings -- ``MARCH1``, ``SEPT2``, ``DEC1``, ``JUN``, ``MAY``, ``INF2`` --
because the date patterns all require a digit/separator structure a bare symbol
never has.

Known, accepted edge: a hypothetical symbol shaped exactly like an Excel date
rendering (e.g. ``"SEP-2"``) would be rejected; no current HGNC/MGI/RGD symbol
has that shape.
"""
from __future__ import annotations

import datetime
import re

#: Prefix used when an invalid Target Gene cell is recovered via its Entrez ID but the
#: workbook never uses a valid symbol for that ID (so no symbol can be resolved without
#: fabricating one). ``ENTREZ:4735`` is honest, unique per gene, and passes
#: :func:`is_valid_gene_symbol`.
ENTREZ_SENTINEL_PREFIX = "ENTREZ:"

_MONTHS = r"(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"

#: String shapes that are date artifacts, not gene symbols. Every pattern requires a
#: digit AND a separator, so bare symbols like ``MARCH1`` / ``SEPT2`` never match.
_DATE_STRING_RES = (
    # ISO date, optionally with a time part: 2022-03-01, 2022-03-01 00:00:00, ...T00:00:00.0
    re.compile(r"^\d{4}-\d{1,2}-\d{1,2}([ T]\d{1,2}:\d{2}(:\d{2})?(\.\d+)?)?$"),
    # numeric date with -, / or . separators: 3/1/2022, 01-03-22, 2022/03/01, 1.3.2022
    re.compile(r"^\d{1,4}[-/.]\d{1,2}[-/.]\d{1,4}$"),
    # day + month name, optional year: 1-Mar, 02-Sep, 2-Sep-22, 1 March 2022
    re.compile(rf"^\d{{1,2}}[-/ ]{_MONTHS}[a-z]*([-/ ]\d{{1,4}})?$", re.IGNORECASE),
    # month name + day/year: Mar-01, March 1, Sep-22, September 2
    re.compile(rf"^{_MONTHS}[a-z]*[-/ ]\d{{1,4}}$", re.IGNORECASE),
)

#: Placeholder strings that mean "no data", not a gene.
_METADATA_VALUES = frozenset({
    "", "-", "--", "n/a", "na", "nan", "none", "null", "unknown", "?", "0",
})

#: Mirrors the webserver's ``Gene.symbol`` = ``CharField(max_length=64)`` so both
#: ends accept the same identity space.
_MAX_SYMBOL_LEN = 64


def is_date_artifact(value: object) -> bool:
    """True iff ``value`` is an Excel date artifact (typed or stringified).

    Covers ``datetime``/``date`` objects and strings matching any of the
    date/timestamp renderings in :data:`_DATE_STRING_RES`. Non-string, non-date
    values (ints, floats, ...) are NOT date artifacts -- numeric rejection is a
    separate concern handled by :func:`is_valid_gene_symbol`.
    """
    if isinstance(value, (datetime.datetime, datetime.date)):
        return True
    if isinstance(value, str):
        s = value.strip()
        return bool(s) and any(rx.match(s) for rx in _DATE_STRING_RES)
    return False


def is_valid_gene_symbol(value: object) -> bool:
    """True iff ``value`` may be used as a gene symbol.

    Accepts the *raw* cell value (any type). Rejects date/datetime objects,
    date-shaped strings, empty/metadata placeholders, purely numeric values, and
    over-long strings. See the module docstring for the full rule list and
    rationale.
    """
    if value is None or isinstance(value, (datetime.datetime, datetime.date)):
        return False
    if isinstance(value, (int, float)):  # bool is an int subclass -> also rejected
        return False
    s = str(value).strip()
    if not s or s.lower() in _METADATA_VALUES:
        return False
    if len(s) > _MAX_SYMBOL_LEN:
        return False
    # Purely numeric in any rendering (int, float, scientific) is an identifier or
    # an artifact, never a symbol. float() also swallows "nan"/"inf" -- both fine
    # to reject.
    try:
        float(s)
        return False
    except ValueError:
        pass
    return not is_date_artifact(s)


def normalize_entrez_id(value: object) -> str | None:
    """Normalize a ``Target Gene (Entrez ID)`` cell to a canonical digit string.

    Handles the renderings openpyxl produces (``4735``, ``4735.0``, ``"4735"``,
    ``"4735.0"``, padded whitespace). Returns ``None`` for anything that is not a
    positive integer -- empty cells, ``"-"``, zero, negatives, non-numeric junk --
    so callers can treat ``None`` as "no recovery path".
    """
    if value is None:
        return None
    if isinstance(value, bool):  # bool is an int subclass; never an Entrez ID
        return None
    if isinstance(value, int):
        return str(value) if value > 0 else None
    if isinstance(value, float):
        return str(int(value)) if value > 0 and value.is_integer() else None
    s = str(value).strip()
    if s.endswith(".0"):
        s = s[:-2]
    if s.isdigit() and int(s) > 0:
        return str(int(s))  # int() round-trip strips leading zeros
    return None


def entrez_sentinel(entrez_id: str) -> str:
    """The symbol used when recovery resolves an Entrez ID but no usable symbol."""
    return f"{ENTREZ_SENTINEL_PREFIX}{entrez_id}"


#: Canonical column order of the release ``exclusion_ledger.tsv`` (pinned to the
#: mirpathdb-data-2026-08-25 release schema; the ledger ships header-only when a
#: build has zero exclusions). Producers: :func:`mirpathdb_prod.sources.mirtarbase.load_mti`
#: (loader lane, via its ``ledger`` argument) and ``scripts/export_target_pairs.py``
#: (export lane: sentinel-policy exclusions). Keep in sync with the webserver twin.
EXCLUSION_LEDGER_COLUMNS: tuple[str, ...] = (
    "record_id",
    "organism",
    "mirna",
    "raw_value",
    "entrez",
    "multiplicity",
    "n_source_rows_without_entrez",
    "annotation_symbols",
    "n_annotation_symbols",
    "recovered_symbol",
    "atomic_evidence",
    "date_form_check",
    "disposition",
    "reason",
    "support_types",
    "workbook_rows",
    "mirtarbase_ids",
)
