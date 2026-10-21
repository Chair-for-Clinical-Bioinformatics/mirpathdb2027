"""Taxonomy resolution (pure).

Maps source organism names/codes to a canonical NCBI taxid + scientific name,
applying curated synonyms. Comparison and identity are taxid-first; this module
provides the curated synonym lookup that makes that resolution deterministic.

The curated maps below are intentionally small and explicit: they cover the
organisms this production build touches — the 4 production species (Homo
sapiens, Mus musculus, Rattus norvegicus, and Macaca mulatta — macaque is
restored as a PREDICTED-DOMINANT production organism: TargetScan feeds LOW (293
miRNAs), with a lone DIANA-TarBase v9 validated MTI marked into HIGH), plus the
evaluated-and-excluded dog (Canis lupus familiaris) and chicken (Gallus gallus)
kept resolvable for the availability-matrix reference rows — plus the obvious
aliases (three-letter codes, common names, and spelling/whitespace variants).
They are NOT a general NCBI Taxonomy mirror — a name that is not curated resolves
to ``None`` so callers can flag it as an ambiguous/unmapped record rather than
silently guessing.
"""

from __future__ import annotations

# Canonical scientific name per taxid. The single authoritative spelling used
# downstream regardless of how a source labelled the organism.
_CANONICAL_NAME: dict[int, str] = {
    9606: "Homo sapiens",
    10090: "Mus musculus",
    10116: "Rattus norvegicus",
    # Macaca mulatta is a PRODUCTION organism (predicted-dominant: TargetScan feeds
    # LOW, with a lone DIANA-TarBase v9 validated MTI marked into HIGH; miRTarBase +
    # miRDB do not cover it).
    9544: "Macaca mulatta",
    # Canis lupus familiaris (dog) + Gallus gallus (chicken) are NOT production
    # organisms (evaluated + EXCLUDED 2026-06-24 — an owner choice to omit, even
    # though TargetScan + miRDB both cover them). Kept resolvable only so the
    # availability matrix can name the excluded reference rows canonically.
    9615: "Canis lupus familiaris",
    9031: "Gallus gallus",
}

# Curated synonym -> canonical taxid. Keys are matched after normalization
# (see :func:`normalize_name`), so casing / surrounding whitespace / internal
# whitespace runs do not matter and underscores are treated as spaces.
_SYNONYM_TO_TAXID: dict[str, int] = {
    # Homo sapiens
    "homo sapiens": 9606,
    "hsa": 9606,
    "human": 9606,
    "h. sapiens": 9606,
    "h sapiens": 9606,
    # Mus musculus
    "mus musculus": 10090,
    "mmu": 10090,
    "mouse": 10090,
    "house mouse": 10090,
    "m. musculus": 10090,
    "m musculus": 10090,
    # Rattus norvegicus
    "rattus norvegicus": 10116,
    "rno": 10116,
    "rat": 10116,
    "norway rat": 10116,
    "brown rat": 10116,
    "r. norvegicus": 10116,
    "r norvegicus": 10116,
    # Canis lupus familiaris (dog) — evaluated + EXCLUDED; resolvable for the
    # availability-matrix reference row only.
    "canis lupus familiaris": 9615,
    "canis familiaris": 9615,
    "cfa": 9615,
    "dog": 9615,
    "domestic dog": 9615,
    "c. lupus familiaris": 9615,
    "c lupus familiaris": 9615,
    # Gallus gallus (chicken) — evaluated + EXCLUDED; resolvable for the
    # availability-matrix reference row only.
    "gallus gallus": 9031,
    "gga": 9031,
    "chicken": 9031,
    "red junglefowl": 9031,
    "g. gallus": 9031,
    "g gallus": 9031,
    # Macaca mulatta (production organism — predicted-dominant; lone DIANA validated MTI marked)
    "macaca mulatta": 9544,
    "mml": 9544,
    "rhesus macaque": 9544,
    "rhesus monkey": 9544,
    "rhesus": 9544,
    "m. mulatta": 9544,
    "m mulatta": 9544,
}


def normalize_name(name: str) -> str:
    """Normalize a free-text organism name/code for synonym lookup.

    Lower-cases, converts underscores to spaces, and collapses internal
    whitespace runs to a single space. Pure and deterministic.

    Args:
        name: Raw organism name, code, or synonym.

    Returns:
        The normalized lookup key (may be empty if ``name`` was blank).
    """
    return " ".join(name.replace("_", " ").split()).lower()


def resolve_tax_id(name: str) -> int | None:
    """Resolve an organism name/synonym to a canonical NCBI taxid.

    Args:
        name: A source organism name, code, or curated synonym.

    Returns:
        The canonical NCBI taxid, or ``None`` if unresolved.
    """
    if not name:
        return None
    return _SYNONYM_TO_TAXID.get(normalize_name(name))


def canonical_scientific_name(tax_id: int) -> str | None:
    """Return the canonical scientific name for a given NCBI taxid.

    Args:
        tax_id: NCBI Taxonomy ID.

    Returns:
        The canonical scientific name, or ``None`` if unknown.
    """
    return _CANONICAL_NAME.get(tax_id)


def is_known_taxid(tax_id: int) -> bool:
    """Return whether ``tax_id`` is in the curated taxonomy."""
    return tax_id in _CANONICAL_NAME
