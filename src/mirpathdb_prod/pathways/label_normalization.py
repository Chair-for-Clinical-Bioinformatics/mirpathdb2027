"""Shared pathway display-label normalization, repair, and integrity gates.

Every pathway display name the pipeline ships (GO, Reactome, inferred Reactome,
WikiPathways) passes through :func:`normalize_label`, and the final label export
is gated by :func:`assert_labels_clean` so an entity-mangled or whitespace-broken
name fails the build loudly instead of shipping (DATA-001).

Root cause of DATA-001 (byte-verified 2026-08-31): the upstream **WikiPathways
GMT export itself** mangles HTML entities in pathway titles — the pinned dated
archive ``data.wikipathways.org/20260610/gmt`` ships ``Parkinson&#39;s disease
pathway`` as ``Parkinson 39 s disease pathway`` (``&``/``#``/``;``/``'`` were
replaced by whitespace upstream, hyphens too: ``omega-3`` ships as ``omega 3``).
The pipeline's downloader and parser are byte-faithful, so the corruption cannot
be repaired by entity decoding alone: the affected ids are repaired from the SAME
release's GPML titles (XML ``Name`` attributes decode correctly) via the
checked-in map ``data/reference/wikipathways_name_repairs.tsv``, regenerable with
``scripts/build_wp_name_repairs.py``. See
``docs/DATA-001_wikipathways_label_root_cause.md`` for the full trace.

The gates run on FINAL display names. If a future upstream release introduces a
new mangling, the build stops here and the repair map must be extended — never
ship a residue, never patch the production database by hand.
"""

from __future__ import annotations

import csv
import html
import re
import unicodedata
from collections.abc import Iterable, Mapping
from pathlib import Path

#: Punctuation that carries meaning in pathway titles and MUST survive
#: normalization (straight/curly apostrophe, hyphen/en-dash, slash, colon,
#: parentheses).
PRESERVED_PUNCTUATION = "'’-–/:()"

#: Default location of the checked-in WikiPathways name-repair map.
DEFAULT_REPAIRS_PATH = (
    Path(__file__).resolve().parents[3] / "data" / "reference" / "wikipathways_name_repairs.tsv"
)

# Horizontal whitespace (incl. NBSP and friends) collapsed to a single space.
_HORIZONTAL_WS = re.compile(r"[ \t\u00a0\u2007\u2009\u202f]+")

#: Residue patterns that must never appear in a shipped display name. Keys are
#: stable identifiers used in violation reports and tests.
LABEL_RESIDUE_PATTERNS: dict[str, re.Pattern[str]] = {
    "numeric_entity": re.compile(r"&#[0-9]+;"),
    "hex_entity": re.compile(r"&#x[0-9A-Fa-f]+;"),
    "named_entity": re.compile(r"&(?:amp|quot|apos|lt|gt|nbsp);"),
    # The upstream GMT mangling signature: "&#39;s" flattened to "39 s".
    "mangled_apostrophe_39_s": re.compile(r"\b39\s+s\b"),
    "double_space": re.compile(r"  "),
    "control_char": re.compile(r"[\x00-\x1f\x7f]"),
    "replacement_char": re.compile("�"),
}


class LabelIntegrityError(ValueError):
    """A shipped-bound display name failed the residue/whitespace gates."""


def normalize_label(raw: str) -> str:
    """Normalize a display label without destroying meaningful punctuation.

    Order (DATA-001 remedy, fixed by contract):

    1. HTML entity decoding (two passes, so double-encoded ``&amp;#39;`` also
       resolves; a second pass is a no-op on clean text);
    2. Unicode NFC normalization;
    3. horizontal-whitespace collapse to single spaces + strip.

    The characters in :data:`PRESERVED_PUNCTUATION` pass through untouched.
    This function is idempotent.

    Args:
        raw: The label as read from a source or intermediate artifact.

    Returns:
        The normalized display label.
    """
    text = html.unescape(html.unescape(raw))
    text = unicodedata.normalize("NFC", text)
    text = _HORIZONTAL_WS.sub(" ", text).strip()
    return text


def label_violations(label: str) -> list[str]:
    """Return the residue-gate violations for one final display name."""
    hits = [key for key, pat in LABEL_RESIDUE_PATTERNS.items() if pat.search(label)]
    if label != label.strip():
        hits.append("leading_trailing_whitespace")
    if not label:
        hits.append("empty_label")
    return hits


def scan_labels(labels: Iterable[tuple[str, str]]) -> dict[str, list[str]]:
    """Scan ``(pathway_id, display_name)`` pairs; return ``{id: violations}``."""
    out: dict[str, list[str]] = {}
    for pid, name in labels:
        hits = label_violations(name)
        if hits:
            out[pid] = hits
    return out


def assert_labels_clean(labels: Iterable[tuple[str, str]], *, context: str) -> None:
    """Fail the build if any display name carries a residue.

    Args:
        labels: ``(pathway_id, display_name)`` pairs.
        context: Human-readable description of what is being gated (goes into
            the error message so the failing stage is unambiguous).

    Raises:
        LabelIntegrityError: naming every offending id, its violations, and the
            offending text.
    """
    bad = scan_labels(labels)
    if bad:
        lines = "\n".join(f"  {pid}: {', '.join(v)}" for pid, v in sorted(bad.items()))
        raise LabelIntegrityError(
            f"label integrity gate FAILED ({context}): {len(bad)} display name(s) "
            f"carry residue patterns — refusing to ship. Extend "
            f"data/reference/wikipathways_name_repairs.tsv (see "
            f"scripts/build_wp_name_repairs.py) or fix the source.\n{lines}"
        )


def load_name_repairs(
    path: str | Path = DEFAULT_REPAIRS_PATH, *, release: str | None = None
) -> dict[str, str]:
    """Load the checked-in ``pathway_id -> canonical display name`` repair map.

    Args:
        path: The repair TSV (columns: release, pathway_id, corrupt_name,
            display_name, source, source_sha256).
        release: If given, keep only rows pinned to this WikiPathways release.

    Returns:
        Mapping of pathway id to the canonical (GPML-derived) display name.
        Missing file returns an empty map — the gates still protect the build.
    """
    out: dict[str, str] = {}
    try:
        fh = Path(path).open(encoding="utf-8", newline="")
    except OSError:
        return out
    with fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if release is not None and row.get("release") != release:
                continue
            out[row["pathway_id"]] = row["display_name"]
    return out


def apply_name_repairs(
    term2name: Mapping[str, str], repairs: Mapping[str, str]
) -> tuple[dict[str, str], int]:
    """Apply the repair map to a ``id -> name`` mapping.

    A repaired id always takes the canonical name (idempotent: if upstream later
    fixes its exporter the repair becomes a no-op equality).

    Returns:
        ``(repaired mapping, number of names actually changed)``.
    """
    out = dict(term2name)
    changed = 0
    for pid, canonical in repairs.items():
        if pid in out and out[pid] != canonical:
            out[pid] = canonical
            changed += 1
    return out, changed
