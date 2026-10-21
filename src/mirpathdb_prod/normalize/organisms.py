"""Organism normalization (pure).

Turns raw organism records (from :mod:`mirpathdb_prod.sources`) into canonical
:class:`mirpathdb_prod.models.Organism` instances. Deterministic and free of
I/O — identity resolution leans on :mod:`mirpathdb_prod.normalize.taxonomy`.

Identity is established by ``canonical_tax_id`` first. When a record carries a
taxid, it is trusted as authoritative and the canonical scientific name is
backfilled from the curated taxonomy (the original label is preserved in
``source_name`` and any divergence is recorded in ``normalization_notes``).
When a record lacks a taxid, the taxid is resolved from the source name/code via
the curated synonym map. Records that cannot be resolved to a taxid are
surfaced explicitly (a ``None`` taxid in the returned :class:`Organism`) rather
than silently dropped, so validation can flag them as ambiguous/unmapped.
"""

from __future__ import annotations

from typing import Any

from mirpathdb_prod.models import Organism
from mirpathdb_prod.normalize.taxonomy import (
    canonical_scientific_name,
    normalize_name,
    resolve_tax_id,
)


def _coerce_taxid(value: Any) -> int | None:
    """Best-effort coercion of a raw taxid field to int (None if blank/invalid)."""
    if value is None:
        return None
    if isinstance(value, bool):  # guard: bool is a subclass of int
        return None
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def normalize_organism(record: dict[str, Any], *, source: str, source_version: str) -> Organism:
    """Normalize a single raw organism record into an :class:`Organism`.

    The raw record may use any of the canonical-schema column names
    (``canonical_tax_id``, ``scientific_name``, ``source_name``, ``common_name``,
    ``source_record_id``, ``normalization_notes``) and/or the prototype-style
    keys (``code``/``display_name``/``common_name``). Identity is taxid-first:

    * If a taxid is present and resolvable, it is authoritative and the canonical
      scientific name is taken from the curated taxonomy.
    * Otherwise the taxid is resolved from the source name or code via the
      curated synonym map.

    Original source labels are always preserved in ``source_name`` (audit field)
    and never used as the basis of equality.

    Args:
        record: Raw record from a source adapter.
        source: Originating source identifier (e.g. "old_prototype", "mirpathdb2").
        source_version: Version/snapshot of that source.

    Returns:
        A canonical :class:`Organism`. ``canonical_tax_id`` is ``0`` (a sentinel
        meaning "unresolved") when no taxid could be established — callers/
        validation must treat a ``0`` taxid as an ambiguous/unmapped record.
    """
    # The label as it appeared in the source: prefer an explicit source_name,
    # then a prototype-style code, then the scientific/display name. This is an
    # audit field and is preserved verbatim.
    source_name = str(
        record.get("source_name")
        or record.get("code")
        or record.get("scientific_name")
        or record.get("display_name")
        or ""
    ).strip()

    raw_sci = str(record.get("scientific_name") or record.get("display_name") or "").strip()
    common_name = str(record.get("common_name") or "").strip()
    source_record_id = str(record.get("source_record_id") or record.get("code") or "").strip()
    incoming_notes = str(record.get("normalization_notes") or "").strip()

    notes: list[str] = []

    # --- taxid resolution (identity-first) ---
    tax_id = _coerce_taxid(record.get("canonical_tax_id") or record.get("tax_id"))
    if tax_id is None:
        # Resolve from the most specific available label: scientific name first,
        # then the source code/name.
        for candidate in (raw_sci, source_name, source_record_id):
            if candidate:
                tax_id = resolve_tax_id(candidate)
                if tax_id is not None:
                    notes.append(f"taxid backfilled from label {candidate!r}")
                    break

    # --- canonical scientific name resolution ---
    canonical = canonical_scientific_name(tax_id) if tax_id is not None else None
    if canonical is not None:
        scientific_name = canonical
        if raw_sci and normalize_name(raw_sci) != normalize_name(canonical):
            notes.append(
                f"scientific_name normalized {raw_sci!r} -> {canonical!r} (taxid {tax_id})"
            )
    else:
        # Unknown/unmapped taxid: fall back to whatever name the source gave,
        # and flag it for validation.
        scientific_name = raw_sci
        if tax_id is None:
            notes.append(f"UNRESOLVED: no canonical taxid for source_name {source_name!r}")
        else:
            notes.append(f"taxid {tax_id} not in curated taxonomy")

    if incoming_notes:
        notes.append(f"source_note: {incoming_notes}")

    return Organism(
        canonical_tax_id=tax_id if tax_id is not None else 0,
        scientific_name=scientific_name,
        source_name=source_name,
        common_name=common_name,
        source=source,
        source_version=source_version,
        source_record_id=source_record_id,
        normalization_notes="; ".join(notes),
    )


def normalize_organisms(
    records: list[dict[str, Any]], *, source: str, source_version: str
) -> list[Organism]:
    """Normalize a list of raw records, preserving input order."""
    return [
        normalize_organism(r, source=source, source_version=source_version) for r in records
    ]
