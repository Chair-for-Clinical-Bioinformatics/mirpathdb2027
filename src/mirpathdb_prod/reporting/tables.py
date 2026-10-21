"""Tabular output rendering (deterministic).

Renders canonical models / validation results as stable TSV or Markdown tables
suitable for committing and diffing (e.g. the organism tables).
"""

from __future__ import annotations

import io
from collections.abc import Iterable

from mirpathdb_prod.models import Organism
from mirpathdb_prod.sources import CANONICAL_ORGANISM_COLUMNS


def _organism_row(org: Organism) -> list[str]:
    """Map an :class:`Organism` to a row in canonical-column order."""
    return [
        str(org.canonical_tax_id),
        org.scientific_name,
        org.source_name,
        org.common_name,
        org.source,
        org.source_version,
        org.source_record_id,
        org.normalization_notes,
    ]


def render_organism_tsv(organisms: Iterable[Organism]) -> str:
    """Render organisms as a deterministic TSV string.

    Uses the canonical 8-column header shared by every organism table in the
    repo. Rows are emitted in the order given by the caller (callers that need a
    stable file should pre-sort, e.g. by ``canonical_tax_id``). A trailing
    newline is always written. Pure: no I/O.

    Args:
        organisms: Canonical organism records.

    Returns:
        TSV text with a stable header and the canonical column order.
    """
    buf = io.StringIO()
    buf.write("\t".join(CANONICAL_ORGANISM_COLUMNS))
    buf.write("\n")
    for org in organisms:
        buf.write("\t".join(_organism_row(org)))
        buf.write("\n")
    return buf.getvalue()
