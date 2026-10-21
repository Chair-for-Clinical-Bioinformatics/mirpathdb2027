"""Tests for the shared label normalization, repair map, and residue gates (DATA-001)."""

from __future__ import annotations

import re

import pytest

from mirpathdb_prod.pathways.label_normalization import (
    DEFAULT_REPAIRS_PATH,
    LabelIntegrityError,
    PRESERVED_PUNCTUATION,
    apply_name_repairs,
    assert_labels_clean,
    label_violations,
    load_name_repairs,
    normalize_label,
    scan_labels,
)

# The seven shipped-corrupt labels of DATA-001 and their canonical GPML titles.
DATA_001_CASES = {
    "WP2059": ("Alzheimer 39 s disease and miRNA effects", "Alzheimer's disease and miRNA effects"),
    "WP2075": ("Alzheimer 39 s disease", "Alzheimer's disease"),
    "WP2371": ("Parkinson 39 s disease pathway", "Parkinson's disease pathway"),
    "WP3638": ("Parkinson 39 s disease", "Parkinson's disease"),
    "WP3853": ("ERK pathway in Huntington 39 s disease", "ERK pathway in Huntington's disease"),
    "WP5124": ("Alzheimer 39 s disease", "Alzheimer's disease"),
    "WP5470": (
        "Effect of omega 3 PUFA on Huntington 39 s disease pathways",
        "Effect of omega-3 PUFA on Huntington's disease pathways",
    ),
}


class TestNormalizeLabel:
    def test_decodes_numeric_entity(self):
        assert normalize_label("Parkinson&#39;s disease") == "Parkinson's disease"

    def test_decodes_double_encoded_entity(self):
        assert normalize_label("Parkinson&amp;#39;s disease") == "Parkinson's disease"

    def test_decodes_named_entities(self):
        assert normalize_label("A &amp; B &quot;x&quot;") == 'A & B "x"'

    def test_nfc_normalization(self):
        # e + COMBINING ACUTE -> precomposed
        assert normalize_label("café") == "café"

    def test_collapses_whitespace(self):
        assert normalize_label("  a  b\t c  ") == "a b c"

    def test_preserves_meaningful_punctuation(self):
        s = f"x{PRESERVED_PUNCTUATION}x"
        assert normalize_label(s) == s

    def test_idempotent(self):
        for raw, _canon in DATA_001_CASES.values():
            once = normalize_label(raw)
            assert normalize_label(once) == once

    def test_does_not_invent_repairs(self):
        # Normalization alone CANNOT fix the upstream mangling (that is the
        # repair map's job) — but it must not make it worse either.
        assert normalize_label("Parkinson 39 s disease") == "Parkinson 39 s disease"


class TestResidueGates:
    @pytest.mark.parametrize(
        ("label", "expected"),
        [
            ("Parkinson&#39;s disease", "numeric_entity"),
            ("x&#x27;y", "hex_entity"),
            ("A &amp; B", "named_entity"),
            ("Parkinson 39 s disease pathway", "mangled_apostrophe_39_s"),
            ("double  space", "double_space"),
            ("ctrl\x01char", "control_char"),
            ("bad � char", "replacement_char"),
        ],
    )
    def test_detects(self, label, expected):
        assert expected in label_violations(label)

    def test_clean_labels_pass(self):
        for _pid, (_raw, canon) in DATA_001_CASES.items():
            assert label_violations(canon) == []

    def test_assert_raises_with_ids(self):
        with pytest.raises(LabelIntegrityError) as exc:
            assert_labels_clean(
                [("WP2371", "Parkinson 39 s disease pathway"), ("GO:1", "fine")],
                context="unit test",
            )
        assert "WP2371" in str(exc.value)
        assert "unit test" in str(exc.value)

    def test_scan_returns_only_violations(self):
        out = scan_labels([("a", "clean name"), ("b", "bad  name")])
        assert set(out) == {"b"}

    def test_briefing_regex_equivalence(self):
        # The briefing's literal gate patterns must all be covered.
        briefing = re.compile(r"&#[0-9]+;|&#x[0-9A-Fa-f]+;|&amp;|&quot;|\b39\s+s\b")
        for probe in ("&#39;", "&#x27;", "&amp;", "&quot;", "x 39 s y"):
            assert briefing.search(probe)
            assert label_violations(probe)


class TestRepairMap:
    def test_checked_in_map_covers_all_seven(self):
        repairs = load_name_repairs(DEFAULT_REPAIRS_PATH, release="20260610")
        assert set(repairs) == set(DATA_001_CASES)
        for pid, (_raw, canon) in DATA_001_CASES.items():
            assert repairs[pid] == canon
            assert label_violations(canon) == []

    def test_apply_fixes_and_counts(self):
        term2name = {pid: raw for pid, (raw, _c) in DATA_001_CASES.items()}
        term2name["WP9999"] = "Untouched pathway"
        repairs = load_name_repairs(DEFAULT_REPAIRS_PATH, release="20260610")
        fixed, n = apply_name_repairs(term2name, repairs)
        assert n == len(DATA_001_CASES)
        assert fixed["WP9999"] == "Untouched pathway"
        assert_labels_clean(fixed.items(), context="repaired map")

    def test_apply_is_idempotent(self):
        repairs = load_name_repairs(DEFAULT_REPAIRS_PATH, release="20260610")
        base = {pid: canon for pid, (_r, canon) in DATA_001_CASES.items()}
        fixed, n = apply_name_repairs(base, repairs)
        assert n == 0
        assert fixed == base

    def test_release_filter(self):
        assert load_name_repairs(DEFAULT_REPAIRS_PATH, release="19990101") == {}

    def test_missing_file_is_empty(self):
        assert load_name_repairs("/nonexistent/nowhere.tsv") == {}
