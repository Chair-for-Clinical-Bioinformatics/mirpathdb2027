"""Tests for the GMT parser + family classifier (pathways/gmt.py).

Covers: family classification from v2-style filenames, standard GMT line
parsing (name / description / members), term-name de-dup/union, directory
splitting by family, and the TERM2GENE long-form flattener. Uses synthetic
fixtures only (no R, no network).
"""

from __future__ import annotations

import pytest

from mirpathdb_prod.pathways import gmt


def _write_gmt(path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_classify_family_v2_filenames() -> None:
    assert gmt.classify_family("GO_BP_validated_miRTarBase_strong.gmt") == gmt.FAMILY_GO_BP
    assert gmt.classify_family("GO_CC_validated_miRTarBase_all.gmt") == gmt.FAMILY_GO_CC
    assert gmt.classify_family("GO_MF_prediction_union.gmt") == gmt.FAMILY_GO_MF
    assert gmt.classify_family("REACTOME_validated_miRTarBase_strong.gmt") == gmt.FAMILY_REACTOME
    assert (
        gmt.classify_family("WIKIPATHWAYS_validated_miRTarBase_all.gmt")
        == gmt.FAMILY_WIKIPATHWAYS
    )
    assert gmt.classify_family("KEGG_validated_miRTarBase_strong.gmt") == gmt.FAMILY_KEGG


def test_classify_family_prefers_longest_token() -> None:
    # "GO_BP" must win over a bare "GO"; full path is fine (basename inspected).
    assert gmt.classify_family("/some/dir/GO_BP_x.gmt") == gmt.FAMILY_GO_BP


def test_classify_family_rejects_unknown() -> None:
    with pytest.raises(ValueError, match="cannot classify"):
        gmt.classify_family("PANTHER_validated.gmt")
    # A family must be followed by a separator: "GO_BPX" is not GO_BP.
    with pytest.raises(ValueError):
        gmt.classify_family("GO_BPX_strong.gmt")


def test_parse_gmt_members_and_description(tmp_path) -> None:
    f = tmp_path / "GO_MF_validated_miRTarBase_strong.gmt"
    _write_gmt(
        f,
        [
            "ATP binding\thttp://example/ATP\tGeneA\tGeneB\tGeneC",
            "DNA binding\thttp://example/DNA\tGeneB",
            "\t\t",  # blank-ish line: dropped (no name)
        ],
    )
    parsed = gmt.parse_gmt(f)
    assert parsed["ATP binding"] == {"GeneA", "GeneB", "GeneC"}
    assert parsed["DNA binding"] == {"GeneB"}
    assert "" not in parsed  # the empty line did not create a term


def test_parse_gmt_unions_repeated_terms(tmp_path) -> None:
    f = tmp_path / "GO_BP_validated_miRTarBase_all.gmt"
    _write_gmt(
        f,
        [
            "apoptosis\tdesc1\tG1\tG2",
            "apoptosis\tdesc2\tG2\tG3",  # same term repeated -> union
        ],
    )
    parsed = gmt.parse_gmt(f)
    assert parsed["apoptosis"] == {"G1", "G2", "G3"}


def test_parse_gmt_drops_blank_member_fields(tmp_path) -> None:
    f = tmp_path / "WIKIPATHWAYS_validated_miRTarBase_strong.gmt"
    # trailing tab => empty member field, must be dropped.
    _write_gmt(f, ["pathway X\tdesc\tGeneA\t\tGeneB\t"])
    parsed = gmt.parse_gmt(f)
    assert parsed["pathway X"] == {"GeneA", "GeneB"}


def test_parse_gmt_dir_splits_by_family(tmp_path) -> None:
    _write_gmt(tmp_path / "GO_BP_validated_miRTarBase_strong.gmt", ["t1\td\tA\tB"])
    _write_gmt(tmp_path / "GO_MF_validated_miRTarBase_strong.gmt", ["t2\td\tC"])
    _write_gmt(tmp_path / "WIKIPATHWAYS_validated_miRTarBase_all.gmt", ["t3\td\tD\tE"])
    # an unclassifiable file is skipped, not an error
    _write_gmt(tmp_path / "README_notes.gmt", ["x\td\tZ"])

    by_family = gmt.parse_gmt_dir(tmp_path)
    assert set(by_family) == {gmt.FAMILY_GO_BP, gmt.FAMILY_GO_MF, gmt.FAMILY_WIKIPATHWAYS}
    assert by_family[gmt.FAMILY_GO_BP]["t1"] == {"A", "B"}
    assert by_family[gmt.FAMILY_GO_MF]["t2"] == {"C"}
    assert by_family[gmt.FAMILY_WIKIPATHWAYS]["t3"] == {"D", "E"}


def test_parse_gmt_dir_missing_dir(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        gmt.parse_gmt_dir(tmp_path / "nope")


def test_iter_gmt_missing_file(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        list(gmt.iter_gmt(tmp_path / "nope.gmt"))


def test_term2gene_rows_deterministic() -> None:
    rows = gmt.term2gene_rows({"T2": {"b", "a"}, "T1": {"z"}})
    # terms sorted, genes sorted within a term
    assert rows == [("T1", "z"), ("T2", "a"), ("T2", "b")]


def test_publishable_families_exclude_kegg() -> None:
    assert gmt.FAMILY_KEGG in gmt.KNOWN_FAMILIES
    assert gmt.FAMILY_KEGG not in gmt.PUBLISHABLE_FAMILIES
