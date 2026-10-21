"""Core data models for the production build.

These are plain dataclasses (no ORM, no Django) so the database layer stays
pure and testable. The :class:`Organism` model carries the canonical fields the
briefing requires for organism concordance: identity is established by
``canonical_tax_id`` first, then ``scientific_name`` + curated synonyms — never
the display label alone.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Organism:
    """A single organism record from one source, in canonical form.

    Concordance between sources (prototype vs. new production vs. miRPathDB2) is
    decided primarily on ``canonical_tax_id``, falling back to
    ``scientific_name`` matched against curated synonyms. ``source_name`` /
    ``common_name`` / display labels are informational and must NOT be the basis
    of equality.

    Attributes:
        canonical_tax_id: NCBI Taxonomy ID — the primary identity key
            (e.g. 9606 for Homo sapiens). The single most authoritative field.
        scientific_name: Canonical binomial scientific name
            (e.g. "Homo sapiens").
        source_name: The organism name exactly as it appeared in the source
            (may differ from ``scientific_name`` before normalization).
        common_name: Common/vernacular name (e.g. "Human"); informational only.
        source: Identifier of the originating source
            (e.g. "prototype", "mirpathdb2", "legacy").
        source_version: Version/snapshot of that source (e.g. "v2", "v3-proto").
        source_record_id: The source's own record identifier for this organism
            (e.g. a three-letter code like "hsa", or a source row id).
        normalization_notes: Free-text notes on how this record was normalized
            (e.g. synonym applied, taxid backfilled). Aids auditability.
    """

    canonical_tax_id: int
    scientific_name: str
    source_name: str
    common_name: str = ""
    source: str = ""
    source_version: str = ""
    source_record_id: str = ""
    normalization_notes: str = ""

    def identity_key(self) -> int:
        """Return the primary identity key used for cross-source comparison."""
        return self.canonical_tax_id
