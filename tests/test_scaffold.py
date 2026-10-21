"""Smoke tests for the production-build scaffold.

These check the foundation (package imports, CLI wiring, config loading, the
Organism model's canonical fields). They deliberately do NOT test parity logic
— that is owned by later phases.
"""

from __future__ import annotations

import dataclasses

from mirpathdb_prod.cli import build_parser
from mirpathdb_prod.config import load_config
from mirpathdb_prod.models import Organism

# Canonical Organism fields required by the briefing.
EXPECTED_ORGANISM_FIELDS = {
    "canonical_tax_id",
    "scientific_name",
    "source_name",
    "common_name",
    "source",
    "source_version",
    "source_record_id",
    "normalization_notes",
}

# Production organism set, keyed by canonical NCBI taxid (LOCKED DECISION
# 2026-06-24: EXACTLY 4 species — hsa/mmu floor + rno (rat, validated+predicted) +
# mml (macaque, PREDICTED-DOMINANT with a lone validated MTI marked) extensions.
# Dog (9615) and chicken (9031) were evaluated + EXCLUDED — owner-omitted, even
# though TargetScan + miRDB cover them).
EXPECTED_TAXIDS = {9606, 10090, 10116, 9544}


def test_organism_has_canonical_fields() -> None:
    fields = {f.name for f in dataclasses.fields(Organism)}
    assert fields == EXPECTED_ORGANISM_FIELDS


def test_organism_identity_is_taxid() -> None:
    org = Organism(canonical_tax_id=9606, scientific_name="Homo sapiens", source_name="hsa")
    assert org.identity_key() == 9606


def test_config_loads_production_organisms() -> None:
    cfg = load_config()
    taxids = {o["canonical_tax_id"] for o in cfg.organisms}
    assert taxids == EXPECTED_TAXIDS
    by_code = {o["code"]: o for o in cfg.organisms}
    # macaque (mml) is RESTORED as a production organism (predicted-dominant; the
    # lone DIANA-TarBase v9 validated MTI is marked, not buried).
    assert "mml" in by_code
    # dog (cfa) and chicken (gga) were evaluated and EXCLUDED — not production.
    assert "cfa" not in by_code
    assert "gga" not in by_code
    # in_mirpathdb2 is true only for the floor (hsa, mmu).
    assert by_code["hsa"]["in_mirpathdb2"] is True
    assert by_code["mmu"]["in_mirpathdb2"] is True
    assert by_code["rno"]["in_mirpathdb2"] is False
    assert by_code["mml"]["in_mirpathdb2"] is False
    # depth classes carry the target-evidence depth.
    assert by_code["hsa"]["depth_class"] == "full"
    assert by_code["mmu"]["depth_class"] == "full"
    assert by_code["rno"]["depth_class"] == "validated_predicted"
    assert by_code["mml"]["depth_class"] == "predicted_dominant"
    # native_sources is a non-empty list per organism.
    for o in cfg.organisms:
        assert isinstance(o["native_sources"], list) and o["native_sources"]


def test_config_macaque_is_predicted_with_lone_validated() -> None:
    # Macaque (mml, 9544) is a PREDICTED-DOMINANT production organism: its
    # substantive evidence is TargetScan-PREDICTED (LOW, 293 miRNAs). DIANA-TarBase
    # v9 contributes EXACTLY ONE validated MTI (mml-miR-125b-5p -> ANGPT2), which is
    # MARKED (surfaced, not buried) — HIGH = that single MTI, MEDIUM = none, LOW =
    # TargetScan. miRTarBase and miRDB do NOT cover macaque.
    cfg = load_config()
    by_code = {o["code"]: o for o in cfg.organisms}
    mml = by_code["mml"]
    assert mml["canonical_tax_id"] == 9544
    assert mml["scientific_name"] == "Macaca mulatta"
    assert mml["depth_class"] == "predicted_dominant"
    assert [s.lower() for s in mml["native_sources"]] == ["diana_tarbase", "targetscan"]
    # The lone validated MTI is MARKED, not hidden behind a "validated" label.
    assert mml["lone_validated_mti"] == "mml-miR-125b-5p -> ANGPT2"
    assert mml["validated_mti_count"] == 1
    # mml is NOT a documented removal (no dropped_organisms section anymore) and
    # NOT in evaluated_excluded.
    assert not cfg.raw.get("dropped_organisms")
    excluded = {e["canonical_tax_id"] for e in cfg.raw.get("evaluated_excluded", [])}
    assert 9544 not in excluded


def test_config_documents_cfa_gga_exclusion() -> None:
    # Dog + chicken are evaluated + EXCLUDED, recorded with a rationale (and in
    # NEITHER the approved_extensions nor approved_removals allowlist).
    cfg = load_config()
    excluded = cfg.raw.get("evaluated_excluded", [])
    by_tax = {e["canonical_tax_id"]: e for e in excluded}
    assert set(by_tax) == {9615, 9031}, "config must record dog + chicken as evaluated/excluded"
    for taxid, code in ((9615, "cfa"), (9031, "gga")):
        e = by_tax[taxid]
        assert e["code"] == code
        assert e.get("in_mirpathdb2") is False
        assert e.get("in_prototype") is False
        assert "exclud" in e["rationale"].lower()
    op = cfg.raw.get("validation", {}).get("organism_parity", {})
    ext = {int(x) for x in op.get("approved_extensions", [])}
    rem = {int(x) for x in op.get("approved_removals", [])}
    # cfa/gga belong in NEITHER allowlist.
    assert 9615 not in ext and 9615 not in rem
    assert 9031 not in ext and 9031 not in rem


def test_config_parity_policy_floor_expand() -> None:
    cfg = load_config()
    op = cfg.raw.get("validation", {}).get("organism_parity", {})
    assert sorted(int(x) for x in op["hard_floor"]) == [9606, 10090]
    # rat (10116) + macaque (9544) extensions now.
    assert sorted(int(x) for x in op["approved_extensions"]) == [9544, 10116]
    # no removals — macaque is restored as an approved extension.
    assert [int(x) for x in op["approved_removals"]] == []


def test_config_resolves_reference_paths() -> None:
    cfg = load_config()
    for name in ("prototype_dir", "dragon_backup_v2_dir", "paper_dir", "provenance_dir"):
        assert name in cfg.paths
        assert cfg.path(name).is_absolute()


def test_cli_exposes_all_subcommands() -> None:
    parser = build_parser()
    subactions = [
        a
        for a in parser._actions
        if getattr(a, "choices", None) and "inventory-prototype" in a.choices
    ]
    assert subactions, "no subparsers registered"
    commands = set(subactions[0].choices)
    expected = {
        "inventory-prototype",
        "extract-legacy-organisms",
        "build-mirpathdb2-reference",
        "build-new-organisms",
        "validate-organism-parity",
        "validate-census",
        "validate-diana-methods",
    }
    assert expected <= commands


def test_inventory_prototype_still_stub(capsys) -> None:
    # `inventory-prototype` is out of scope for the organism-parity phase and is
    # still deliberately unimplemented. It now exits cleanly (rc 2 + stderr message)
    # rather than raising an unguarded NotImplementedError. (The organism commands
    # ARE implemented and are exercised by the dedicated parity tests.)
    parser = build_parser()
    args = parser.parse_args(["inventory-prototype"])
    rc = args.func(args)
    assert rc == 2
    err = capsys.readouterr().err
    assert "not yet implemented" in err.lower()
    assert "Traceback" not in err
