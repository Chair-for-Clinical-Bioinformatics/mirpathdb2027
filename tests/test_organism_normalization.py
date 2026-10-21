"""Tests for organism normalization (taxid/synonym compare; audit fields).

These verify the pure normalization layer: identity is taxid-first, the
canonical scientific name is backfilled from the curated taxonomy, curated
synonyms resolve to the right taxid, and the original source labels are
preserved as audit fields (never used as the basis of equality).
"""

from __future__ import annotations

from mirpathdb_prod.normalize.organisms import normalize_organism, normalize_organisms
from mirpathdb_prod.normalize.taxonomy import (
    canonical_scientific_name,
    normalize_name,
    resolve_tax_id,
)


def test_resolve_tax_id_synonyms() -> None:
    # Codes, common names, and binomials all resolve to the same taxid.
    assert resolve_tax_id("hsa") == 9606
    assert resolve_tax_id("Homo sapiens") == 9606
    assert resolve_tax_id("human") == 9606
    assert resolve_tax_id("MMU") == 10090  # case-insensitive
    assert resolve_tax_id("Mouse") == 10090
    assert resolve_tax_id("mml") == 9544
    assert resolve_tax_id("Rhesus macaque") == 9544
    assert resolve_tax_id("Macaca  mulatta") == 9544  # collapsed whitespace


def test_resolve_tax_id_unknown_is_none() -> None:
    assert resolve_tax_id("") is None
    assert resolve_tax_id("Drosophila melanogaster") is None
    assert resolve_tax_id("dme") is None


def test_normalize_name_collapses_and_lowercases() -> None:
    assert normalize_name("  Homo   sapiens ") == "homo sapiens"
    assert normalize_name("Mus_musculus") == "mus musculus"


def test_canonical_scientific_name() -> None:
    assert canonical_scientific_name(9606) == "Homo sapiens"
    assert canonical_scientific_name(10090) == "Mus musculus"
    assert canonical_scientific_name(9544) == "Macaca mulatta"
    assert canonical_scientific_name(7227) is None


def test_normalize_preserves_audit_fields() -> None:
    # A record with an explicit taxid: identity is the taxid; the original code
    # ("hsa") and common name are preserved as audit fields, NOT used for identity.
    rec = {
        "canonical_tax_id": 9606,
        "scientific_name": "Homo sapiens",
        "source_name": "hsa",
        "common_name": "Human",
        "source_record_id": "queries.py::ORGANISMS[0]",
    }
    org = normalize_organism(rec, source="old_prototype", source_version="vX")
    assert org.canonical_tax_id == 9606
    assert org.identity_key() == 9606
    assert org.scientific_name == "Homo sapiens"
    assert org.source_name == "hsa"  # audit field preserved verbatim
    assert org.common_name == "Human"
    assert org.source == "old_prototype"
    assert org.source_version == "vX"
    assert org.source_record_id == "queries.py::ORGANISMS[0]"


def test_normalize_backfills_taxid_from_code() -> None:
    # No taxid given: resolved from the source code, canonical name backfilled.
    rec = {"code": "mml", "display_name": "Macaca mulatta", "common_name": "Rhesus macaque"}
    org = normalize_organism(rec, source="prototype", source_version="v3")
    assert org.canonical_tax_id == 9544
    assert org.scientific_name == "Macaca mulatta"
    assert org.source_name == "mml"  # original code preserved
    assert "backfilled" in org.normalization_notes


def test_normalize_canonicalizes_divergent_name_but_keeps_original() -> None:
    # A divergent display label is canonicalized for identity, but the original
    # label is preserved and the divergence is noted for audit.
    rec = {"canonical_tax_id": 10090, "scientific_name": "mus musculus", "source_name": "mmu"}
    org = normalize_organism(rec, source="s", source_version="v")
    assert org.canonical_tax_id == 10090
    assert org.scientific_name == "Mus musculus"  # canonical spelling
    assert org.source_name == "mmu"
    # name differed only in case -> normalized names match -> no divergence note
    assert "normalized" not in org.normalization_notes


def test_normalize_flags_unresolved_record() -> None:
    rec = {"source_name": "Drosophila melanogaster", "scientific_name": "Drosophila melanogaster"}
    org = normalize_organism(rec, source="s", source_version="v")
    assert org.canonical_tax_id == 0  # sentinel: unresolved
    assert "UNRESOLVED" in org.normalization_notes


def test_normalize_organisms_preserves_order() -> None:
    recs = [
        {"code": "mmu"},
        {"code": "hsa"},
        {"code": "mml"},
    ]
    orgs = normalize_organisms(recs, source="s", source_version="v")
    assert [o.canonical_tax_id for o in orgs] == [10090, 9606, 9544]
