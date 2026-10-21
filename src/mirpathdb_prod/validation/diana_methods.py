"""DIANA-TarBase experimental-method census + strict release gate.

DIANA-TarBase v9 rows carry an ``experimental_method`` and a ``confidence``. The
production HIGH/MEDIUM split (:func:`mirpathdb_prod.sources.tarbase.rollup_tier`)
is a *hybrid* rule: a row is HIGH iff ``confidence == "1"`` AND its method is
non-empty AND is in neither :data:`~mirpathdb_prod.sources.tarbase.CLIP_METHODS`
nor :data:`~mirpathdb_prod.sources.tarbase.KNOWN_OTHER_METHODS`. The rule is
deliberately **fail-open** for a method outside every vocabulary set: an unlisted
gene-specific assay is not silently demoted, but by the same token an unlisted
*non-pair-specific* method would silently inflate the strongly-validated tier
(exactly what ``pSILAC`` / ``SILAC`` / ``Biotin-Seq`` / ... did until the
2026-08-30 vocabulary revision classified them). ``validate_vocabulary`` only
*warns* on such a method, and warnings do not block a release.

This module turns that warning into an auditable artefact and a hard gate:

* **Census** -- every distinct *normalized* method observed in the cached v9 files
  becomes a :class:`MethodRecord`: its vocabulary class
  (:func:`~mirpathdb_prod.sources.tarbase.classify_method`), the raw spellings that
  collapsed onto it, the confidence values and effective tiers seen, and row / unique
  pair / miRNA / gene counts, overall and per organism.
* **Gate** -- ``fail`` iff at least one ``unknown`` method contributes HIGH rows;
  ``warn`` iff unknown methods exist but none is HIGH; else ``pass``. Nothing is
  reclassified automatically and the hybrid rule is NOT replaced by an allowlist:
  the gate makes the owner add each unknown method to the right vocabulary set
  deliberately before a release ships.

The census classifies rows with the SAME functions the band builders use
(``normalize_method`` / ``normalize_confidence`` / ``rollup_tier`` /
``classify_method`` from :mod:`mirpathdb_prod.sources.tarbase`, and the
``_clean`` identifier backstop from
:mod:`mirpathdb_prod.sources.production_target_sets` for pair counting), so it
cannot disagree with production about which rows are HIGH. Rows are loaded with
the production loaders (``load_mti_for_taxid``), so the Other_species species
slice and the vocabulary guards apply exactly as in a build.

Artefacts (all deterministic: sorted keys / rows, no timestamps) under ``out_dir``::

    diana_method_audit.tsv           one row per normalized method
    diana_method_audit.json          {"summary": ..., "methods": [...]}
    diana_method_audit_summary.json  the summary block alone
    diana_method_audit.md            human-readable report (unknown methods first)

plus, unless suppressed, ``<provenance_dir>/diana_method_audit.json`` -- the
provenance block that pins the source sha256s, the classifier rule, the vocabulary
sets and the report's sha256 next to the gate verdict.

The empty method string is recorded as its own census entry (class ``unknown``,
since it sits in no vocabulary set) so the inventory is complete, but it is NOT
counted as an unknown method by the gate, the ``unknown_methods`` lists or the
``n_unknown*`` counters: production's ``validate_vocabulary`` deliberately skips
empty methods (``if method:``), it cannot be curated into a vocabulary set, and
production bands it MEDIUM unconditionally, so it could never trip ``fail``
anyway. It is reported separately as ``n_empty_method`` (0 or 1) so the class
counts partition the census exactly::

    n_normalized_methods == n_known_high + n_known_clip + n_known_other
                            + n_unknown + n_empty_method

``n_pairs_high_from_unknown`` is the size of the UNION of ``(organism, mirna,
gene)`` pairs over every unknown HIGH row (a pair supported by two unknown
methods counts once). Only the accumulator can compute it -- the per-method
records cannot dedupe across methods -- so :func:`summarize` records ``None``
(JSON ``null``) when the accumulator's union is not passed in, never a
per-method sum that over-counts.

The TSV joins list-valued cells with a bare ``|`` (the spec'd encoding); a
method or spelling containing ``|`` therefore does not round-trip from the TSV
-- the JSON artefact is the lossless one. No DIANA v9 vocabulary entry contains
``|`` today.
"""

from __future__ import annotations

import csv
import json
import warnings
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from mirpathdb_prod.config import load_config
from mirpathdb_prod.sources import production_target_sets, tarbase
from mirpathdb_prod.validation.provenance import (
    file_sha256,
    pinned_sha256_index,
    verify_cached_sources,
)

# Artefact basenames (under ``out_dir``) and the provenance block basename.
AUDIT_TSV = "diana_method_audit.tsv"
AUDIT_JSON = "diana_method_audit.json"
SUMMARY_JSON = "diana_method_audit_summary.json"
AUDIT_MD = "diana_method_audit.md"
PROVENANCE_JSON = "diana_method_audit.json"

# Gate verdicts recorded in the summary / provenance block.
GATE_PASS = "pass"
GATE_WARN = "warn"
GATE_FAIL = "fail"

# ``source_release`` label: the resource name + the adapter's pinned version.
SOURCE_RELEASE_LABEL = "DIANA-TarBase"
SOURCE_RELEASE = f"{SOURCE_RELEASE_LABEL} {tarbase.SOURCE_VERSION}"

# Canonical processing order of the production organisms (config.yaml order).
ORGANISM_ORDER: tuple[str, ...] = ("hsa", "mmu", "rno", "mml")

# Column order of ``diana_method_audit.tsv``. List-valued cells are ``|``-joined
# with NO escaping of ``|`` inside an item (tab / quote / backslash are escaped by
# the csv writer); see the module docstring.
AUDIT_TSV_COLUMNS: tuple[str, ...] = (
    "method",
    "vocabulary_class",
    "contributes_to_high",
    "effective_tiers",
    "confidences",
    "n_rows",
    "n_rows_high",
    "n_rows_medium",
    "n_pairs",
    "n_pairs_high",
    "n_mirnas",
    "n_genes",
    "organisms",
    "raw_spellings",
)

# DIANA column accessors (same cells the band builders read).
_METHOD_COL = "experimental_method"
_CONFIDENCE_COL = "confidence"
_MIRNA_COL = "mirna_name"
_GENE_COL = "gene_name"

_HIGH = "HIGH"


@dataclass
class MethodRecord:
    """Census entry for ONE normalized ``experimental_method`` value.

    Pair / miRNA / gene counts apply the production ``_clean`` backstop to the
    ``mirna_name`` / ``gene_name`` cells: a row whose miRNA or gene cleans to
    ``None`` (empty, or an Excel-date artifact such as ``"2022-09-02"``) counts as a
    row but NOT as a pair -- exactly like ``diana_bands``. Identifiers are counted
    within an organism; the record-level ``n_pairs`` / ``n_mirnas`` / ``n_genes``
    are the sums over organisms, i.e. unique ``(organism, mirna, gene)`` /
    ``(organism, mirna)`` / ``(organism, gene)`` keys.

    Attributes:
        method: The normalized method (``tarbase.normalize_method``).
        raw_spellings: Sorted distinct raw cell values that normalize to ``method``.
        vocabulary_class: One of ``tarbase.VOCABULARY_CLASSES``.
        confidences: Sorted distinct normalized ``confidence`` values observed.
        effective_tiers: Sorted distinct ``rollup_tier`` outcomes (HIGH / MEDIUM).
        contributes_to_high: True iff at least one row rolled up HIGH.
        n_rows / n_rows_high / n_rows_medium: Row counts.
        n_pairs / n_pairs_high: Unique cleaned (mirna, gene) pairs, all rows / HIGH rows.
        n_mirnas / n_genes: Unique cleaned miRNAs / genes among the counted pairs.
        organisms: Sorted organism codes with >= 1 row.
        per_organism: ``org -> {n_rows, n_rows_high, n_pairs, n_pairs_high,
            n_mirnas, n_genes}``.
    """

    method: str
    raw_spellings: list[str]
    vocabulary_class: str
    confidences: list[str]
    effective_tiers: list[str]
    contributes_to_high: bool
    n_rows: int
    n_rows_high: int
    n_rows_medium: int
    n_pairs: int
    n_pairs_high: int
    n_mirnas: int
    n_genes: int
    organisms: list[str]
    per_organism: dict[str, dict[str, int]] = field(default_factory=dict)


@dataclass
class AuditResult:
    """Outcome of :func:`run_diana_method_audit`.

    Attributes:
        records: The finalized census, sorted by normalized method.
        summary: The summary block (also written to ``diana_method_audit_summary.json``).
        artefacts: ``{artefact name: written path}``.
        exit_code: 0 (pass / warn, or any report-only run) or 1 (strict fail).
    """

    records: list[MethodRecord]
    summary: dict[str, Any]
    artefacts: dict[str, Path]
    exit_code: int


class MethodAccumulator:
    """Streaming per-(method, organism) tallies behind :func:`audit_rows`.

    Keeps the unique pair sets in memory (hsa carries ~4.7M rows; sets of
    2-tuples per (method, organism) are well within reach) so the census can
    report unique pairs / miRNAs / genes rather than raw row counts alone.
    """

    def __init__(self) -> None:
        self._raw_spellings: dict[str, set[str]] = defaultdict(set)
        self._confidences: dict[str, set[str]] = defaultdict(set)
        self._tiers: dict[str, set[str]] = defaultdict(set)
        self._n_rows: Counter[tuple[str, str]] = Counter()
        self._n_rows_high: Counter[tuple[str, str]] = Counter()
        self._pairs: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
        self._pairs_high: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
        self._methods_by_organism: dict[str, set[str]] = defaultdict(set)
        # Union of (organism, mirna, gene) over every HIGH row whose method is
        # unknown -- the spec'd ``n_pairs_high_from_unknown`` (a pair backed by
        # two unknown methods counts once; the per-method records cannot dedupe).
        self._unknown_high_pairs: set[tuple[str, str, str]] = set()
        self.rows_by_organism: Counter[str] = Counter()

    def add(self, row: Mapping[str, Any], *, organism: str) -> None:
        """Tally one raw DIANA row (keyed by ``tarbase.TARBASE_COLUMNS``)."""
        raw = row.get(_METHOD_COL)
        method = tarbase.normalize_method(raw)
        confidence = tarbase.normalize_confidence(row.get(_CONFIDENCE_COL))
        # The production classifier, called on the raw row exactly as diana_bands
        # does (rollup_tier only reads .get) -- never re-implemented, never copied.
        tier = tarbase.rollup_tier(row)  # type: ignore[arg-type]
        key = (method, organism)

        self._raw_spellings[method].add("" if raw is None else str(raw))
        self._confidences[method].add(confidence)
        self._tiers[method].add(tier)
        self._methods_by_organism[organism].add(method)
        self._n_rows[key] += 1
        self.rows_by_organism[organism] += 1
        is_high = tier == _HIGH
        if is_high:
            self._n_rows_high[key] += 1

        # Pair counting mirrors ``diana_bands``: both identifiers must survive the
        # production ``_clean`` backstop, else the row is a row but not a pair.
        mirna = production_target_sets._clean(row.get(_MIRNA_COL))
        gene = production_target_sets._clean(row.get(_GENE_COL))
        if mirna and gene:
            pair = (mirna, gene)
            self._pairs[key].add(pair)
            if is_high:
                self._pairs_high[key].add(pair)
                # Same predicate as the gate / summary counters, so the union can
                # never disagree with n_unknown_contributing_high about what an
                # "unknown" is (today the empty method is never HIGH anyway).
                if _is_gate_unknown_method(method):
                    self._unknown_high_pairs.add((organism, mirna, gene))

    def methods(self) -> list[str]:
        """Sorted normalized methods seen so far (all organisms)."""
        return sorted(self._raw_spellings)

    def methods_for(self, organism: str) -> list[str]:
        """Sorted normalized methods seen in ONE organism."""
        return sorted(self._methods_by_organism.get(organism, set()))

    def pairs_for(self, organism: str, *, high_only: bool = False) -> set[tuple[str, str]]:
        """Union over methods of the cleaned (mirna, gene) pairs of one organism.

        With ``high_only`` this equals the pair set of ``diana_bands(rows)["HIGH"]``
        for the same rows; without it, HIGH | MEDIUM (used by the mirror tests).
        """
        source = self._pairs_high if high_only else self._pairs
        out: set[tuple[str, str]] = set()
        for (_method, org), pairs in source.items():
            if org == organism:
                out |= pairs
        return out

    def n_unknown_high_pairs(self) -> int:
        """Unique ``(organism, mirna, gene)`` pairs contributed HIGH by unknown methods."""
        return len(self._unknown_high_pairs)

    def record(self, method: str) -> MethodRecord:
        """Materialise the :class:`MethodRecord` for one normalized method."""
        organisms = sorted({org for (m, org) in self._n_rows if m == method})
        per_organism: dict[str, dict[str, int]] = {}
        for org in organisms:
            key = (method, org)
            pairs = self._pairs.get(key, set())
            per_organism[org] = {
                "n_rows": self._n_rows[key],
                "n_rows_high": self._n_rows_high[key],
                "n_pairs": len(pairs),
                "n_pairs_high": len(self._pairs_high.get(key, set())),
                "n_mirnas": len({m for m, _g in pairs}),
                "n_genes": len({g for _m, g in pairs}),
            }
        n_rows = sum(v["n_rows"] for v in per_organism.values())
        n_rows_high = sum(v["n_rows_high"] for v in per_organism.values())
        return MethodRecord(
            method=method,
            raw_spellings=sorted(self._raw_spellings[method]),
            vocabulary_class=tarbase.classify_method(method),
            confidences=sorted(self._confidences[method]),
            effective_tiers=sorted(self._tiers[method]),
            contributes_to_high=n_rows_high > 0,
            n_rows=n_rows,
            n_rows_high=n_rows_high,
            n_rows_medium=n_rows - n_rows_high,
            n_pairs=sum(v["n_pairs"] for v in per_organism.values()),
            n_pairs_high=sum(v["n_pairs_high"] for v in per_organism.values()),
            n_mirnas=sum(v["n_mirnas"] for v in per_organism.values()),
            n_genes=sum(v["n_genes"] for v in per_organism.values()),
            organisms=organisms,
            per_organism=per_organism,
        )


def audit_rows(
    rows: Iterable[Mapping[str, Any]], *, organism: str, acc: MethodAccumulator
) -> int:
    """Stream one organism's raw DIANA rows into ``acc``; return the row count."""
    n = 0
    for row in rows:
        acc.add(row, organism=organism)
        n += 1
    return n


def finalize(acc: MethodAccumulator) -> list[MethodRecord]:
    """Materialise the census, sorted by normalized method (plain ``sorted``)."""
    return [acc.record(m) for m in acc.methods()]


def _is_gate_unknown_method(method: str) -> bool:
    """The ONE definition of "unknown for gate purposes": unknown class AND non-empty.

    Mirrors production's ``validate_vocabulary`` (``if method:``): the empty
    method is inventoried but is neither curatable nor ever HIGH, so it is not a
    gate-relevant unknown (see the module docstring). Used by the accumulator's
    unknown-HIGH pair union and by :func:`is_gate_unknown` alike.
    """
    return bool(method) and tarbase.is_unknown_method(method)


def is_gate_unknown(rec: MethodRecord) -> bool:
    """True iff ``rec`` is an unknown method the gate must account for."""
    return rec.vocabulary_class == tarbase.VOCAB_UNKNOWN and _is_gate_unknown_method(rec.method)


def _fmt_count(value: int | None) -> str:
    """Render an optional count (``None`` -> ``n/a``) for the MD / verdict lines."""
    return "n/a" if value is None else str(value)


def gate_status(records: Sequence[MethodRecord]) -> str:
    """``fail`` iff an unknown method is HIGH; ``warn`` iff unknowns exist; else ``pass``."""
    unknown = [r for r in records if is_gate_unknown(r)]
    if any(r.contributes_to_high for r in unknown):
        return GATE_FAIL
    if unknown:
        return GATE_WARN
    return GATE_PASS


def summarize(
    records: Sequence[MethodRecord],
    *,
    source_files: Sequence[Mapping[str, Any]],
    n_pairs_high_from_unknown: int | None = None,
) -> dict[str, Any]:
    """Build the summary block over a finalized census.

    ``n_rows_high_from_unknown`` counts rows. ``n_pairs_high_from_unknown`` is the
    number of unique ``(organism, mirna, gene)`` pairs contributed HIGH by unknown
    methods; pass :meth:`MethodAccumulator.n_unknown_high_pairs` (the driver
    does). Records alone cannot dedupe a pair backed by two unknown methods, so
    when it is not supplied the key is ``None`` (JSON ``null``) rather than an
    over-counting per-method sum. The empty method is not an unknown here (see
    :func:`is_gate_unknown`) but is reported as ``n_empty_method`` so the four
    class counts plus it sum to ``n_normalized_methods``. ``organisms`` is the
    union of the organisms seen in the records and those named by
    ``source_files``, whose entries are ``{organisms, path, sha256, pinned_sha256,
    bytes}`` (see :func:`run_diana_method_audit`).
    """
    by_class: Counter[str] = Counter(r.vocabulary_class for r in records)
    unknown = [r for r in records if is_gate_unknown(r)]
    unknown_high = [r for r in unknown if r.contributes_to_high]
    raw_strings: set[str] = set()
    for r in records:
        raw_strings.update(r.raw_spellings)
    organisms_with_unknown_high = sorted(
        {
            org
            for r in unknown_high
            for org, counts in r.per_organism.items()
            if counts["n_rows_high"] > 0
        }
    )
    organisms = sorted(
        {org for r in records for org in r.organisms}
        | {org for sf in source_files for org in sf.get("organisms", [])}
    )
    per_organism: dict[str, dict[str, int]] = {}
    for org in organisms:
        recs = [r for r in records if org in r.per_organism]
        per_organism[org] = {
            "n_rows": sum(r.per_organism[org]["n_rows"] for r in recs),
            "n_rows_high": sum(r.per_organism[org]["n_rows_high"] for r in recs),
            "n_methods": len(recs),
            "n_unknown_methods": sum(1 for r in recs if is_gate_unknown(r)),
            "n_unknown_high_methods": sum(
                1 for r in recs if is_gate_unknown(r) and r.per_organism[org]["n_rows_high"] > 0
            ),
        }
    # The empty method (at most one record) is neither known nor a gate unknown.
    n_empty_method = sum(1 for r in records if not r.method)
    return {
        "n_raw_method_strings": len(raw_strings),
        "n_normalized_methods": len(records),
        "n_known_high": by_class.get(tarbase.VOCAB_KNOWN_HIGH, 0),
        "n_known_clip": by_class.get(tarbase.VOCAB_KNOWN_CLIP, 0),
        "n_known_other": by_class.get(tarbase.VOCAB_KNOWN_OTHER, 0),
        "n_unknown": len(unknown),
        "n_empty_method": n_empty_method,
        "n_unknown_contributing_high": len(unknown_high),
        "n_rows_high_from_unknown": sum(r.n_rows_high for r in unknown_high),
        "n_pairs_high_from_unknown": (
            None if n_pairs_high_from_unknown is None else int(n_pairs_high_from_unknown)
        ),
        "organisms_with_unknown_high": organisms_with_unknown_high,
        "unknown_methods": sorted(r.method for r in unknown),
        "unknown_high_methods": sorted(r.method for r in unknown_high),
        "classifier_rule": tarbase.DIANA_TIER_RULE,
        "vocabulary_sizes": {
            "high": len(tarbase.HIGH_METHODS),
            "clip": len(tarbase.CLIP_METHODS),
            "known_other": len(tarbase.KNOWN_OTHER_METHODS),
        },
        "gate_status": gate_status(records),
        "source_release": SOURCE_RELEASE,
        "source_files": [dict(sf) for sf in sorted(source_files, key=lambda s: s["path"])],
        "organisms": organisms,
        "per_organism": per_organism,
    }


# --------------------------------------------------------------------------- #
# Writers (deterministic: sorted keys / rows, trailing newline, no timestamps)
# --------------------------------------------------------------------------- #


def _dump_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )


def _tsv_row(rec: MethodRecord) -> list[str]:
    return [
        rec.method,
        rec.vocabulary_class,
        "true" if rec.contributes_to_high else "false",
        "|".join(rec.effective_tiers),
        "|".join(rec.confidences),
        str(rec.n_rows),
        str(rec.n_rows_high),
        str(rec.n_rows_medium),
        str(rec.n_pairs),
        str(rec.n_pairs_high),
        str(rec.n_mirnas),
        str(rec.n_genes),
        "|".join(rec.organisms),
        "|".join(rec.raw_spellings),
    ]


def write_audit_tsv(records: Sequence[MethodRecord], path: Path) -> None:
    """Write one row per method, columns :data:`AUDIT_TSV_COLUMNS` (``|``-joined lists)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(
            fh,
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            lineterminator="\n",
        )
        w.writerow(AUDIT_TSV_COLUMNS)
        for rec in records:
            w.writerow(_tsv_row(rec))


def write_audit_json(
    records: Sequence[MethodRecord], summary: Mapping[str, Any], path: Path
) -> None:
    """Write ``{"summary": ..., "methods": [asdict(record), ...]}``."""
    _dump_json({"summary": dict(summary), "methods": [asdict(r) for r in records]}, path)


def write_summary_json(summary: Mapping[str, Any], path: Path) -> None:
    """Write the summary block alone."""
    _dump_json(dict(summary), path)


def _md_cell(text: str) -> str:
    """Escape a value for a Markdown table cell."""
    return text.replace("|", "\\|").replace("\n", " ")


def _md_list(items: Sequence[str]) -> str:
    """Render a list-valued cell with a literal, escaped ``|`` separator."""
    return "\\|".join(_md_cell(x) for x in items) if items else "(empty)"


def _md_spellings(rec: MethodRecord) -> str:
    # JSON-quote raw spellings so stray whitespace is visible in the report.
    return ", ".join(_md_cell(json.dumps(s, ensure_ascii=False)) for s in rec.raw_spellings)


def write_audit_md(
    records: Sequence[MethodRecord], summary: Mapping[str, Any], path: Path
) -> None:
    """Write the human-readable report: summary, UNKNOWN methods first, then everything."""
    verdict = str(summary["gate_status"]).upper()
    lines: list[str] = [
        f"# {summary['source_release']} experimental-method audit",
        "",
        f"**Gate: {verdict}**"
        + (" (report-only; NOT release-gating)" if summary.get("report_only") else ""),
        "",
        f"- classifier rule: `{summary['classifier_rule']}`",
        f"- organisms: {', '.join(summary['organisms']) or '(none)'}",
        f"- raw method strings: {summary['n_raw_method_strings']}; normalized methods: "
        f"{summary['n_normalized_methods']} (known_high {summary['n_known_high']}, "
        f"known_clip {summary['n_known_clip']}, known_other {summary['n_known_other']}, "
        f"unknown {summary['n_unknown']}, empty-method {summary['n_empty_method']})",
        f"- unknown methods contributing HIGH: {summary['n_unknown_contributing_high']} "
        f"(rows {summary['n_rows_high_from_unknown']}, pairs "
        f"{_fmt_count(summary['n_pairs_high_from_unknown'])}; organisms: "
        f"{', '.join(summary['organisms_with_unknown_high']) or 'none'})",
        f"- vocabulary sizes: high {summary['vocabulary_sizes']['high']}, clip "
        f"{summary['vocabulary_sizes']['clip']}, known_other "
        f"{summary['vocabulary_sizes']['known_other']}",
        "",
        "## Source files",
        "",
        "| organisms | file | bytes | sha256 | pinned sha256 |",
        "|---|---|---:|---|---|",
    ]
    for sf in summary["source_files"]:
        lines.append(
            f"| {', '.join(sf['organisms'])} | `{sf['path']}` | {sf['bytes']} | "
            f"`{sf['sha256']}` | `{sf['pinned_sha256'] or '(unpinned)'}` |"
        )
    lines.append("")
    for org, counts in summary.get("per_organism", {}).items():
        lines.append(
            f"- {org}: rows {counts['n_rows']}, HIGH rows {counts['n_rows_high']}, "
            f"methods {counts['n_methods']}, unknown {counts['n_unknown_methods']} "
            f"(HIGH-contributing {counts['n_unknown_high_methods']})"
        )
    loader_warnings = summary.get("loader_warnings") or {}
    if loader_warnings:
        lines += ["", "## Loader warnings", ""]
        for org in sorted(loader_warnings):
            for msg in loader_warnings[org]:
                lines.append(f"- {org}: {msg}")

    unknown = [r for r in records if is_gate_unknown(r)]
    lines += [
        "",
        f"## Unknown methods ({len(unknown)})",
        "",
    ]
    if unknown:
        lines += [
            "| method | HIGH contributor | confidences | tiers | n_rows | n_rows_high | "
            "n_pairs_high | organisms | raw spellings |",
            "|---|---|---|---|---:|---:|---:|---|---|",
        ]
        for r in unknown:
            flag = "**YES**" if r.contributes_to_high else "no"
            lines.append(
                f"| `{_md_cell(r.method)}` | {flag} | {_md_list(r.confidences)} | "
                f"{_md_list(r.effective_tiers)} | {r.n_rows} | {r.n_rows_high} | "
                f"{r.n_pairs_high} | {', '.join(r.organisms)} | {_md_spellings(r)} |"
            )
    else:
        lines.append("(none)")

    lines += [
        "",
        f"## Full inventory ({len(records)} methods)",
        "",
        "| method | class | HIGH contributor | tiers | confidences | n_rows | n_rows_high | "
        "n_pairs | n_pairs_high | n_mirnas | n_genes | organisms | raw spellings |",
        "|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for r in records:
        lines.append(
            f"| `{_md_cell(r.method)}` | {r.vocabulary_class} | "
            f"{'yes' if r.contributes_to_high else 'no'} | {_md_list(r.effective_tiers)} | "
            f"{_md_list(r.confidences)} | {r.n_rows} | {r.n_rows_high} | {r.n_pairs} | "
            f"{r.n_pairs_high} | {r.n_mirnas} | {r.n_genes} | {', '.join(r.organisms)} | "
            f"{_md_spellings(r)} |"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def provenance_block(
    summary: Mapping[str, Any], *, report_sha256: str, report_path: str
) -> dict[str, Any]:
    """The provenance record written next to the run manifests (no timestamp)."""
    return {
        "classifier_rule": summary["classifier_rule"],
        "source_release": summary["source_release"],
        "source_sha256": {sf["path"]: sf["sha256"] for sf in summary["source_files"]},
        "source_files": [dict(sf) for sf in summary["source_files"]],
        "normalized_method_count": summary["n_normalized_methods"],
        "unknown_method_count": summary["n_unknown"],
        "empty_method_count": summary["n_empty_method"],
        "unknown_high_method_count": summary["n_unknown_contributing_high"],
        "unknown_high_row_count": summary["n_rows_high_from_unknown"],
        "unknown_high_pair_count": summary["n_pairs_high_from_unknown"],
        "report_sha256": report_sha256,
        "report_path": report_path,
        "gate_status": summary["gate_status"],
        "report_only": bool(summary.get("report_only", False)),
        "vocabulary": {
            "high_methods": sorted(tarbase.HIGH_METHODS),
            "clip_methods": sorted(tarbase.CLIP_METHODS),
            "known_other_methods": sorted(tarbase.KNOWN_OTHER_METHODS),
        },
        "organisms": list(summary["organisms"]),
    }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def _ordered(organisms: Sequence[str]) -> list[str]:
    """Canonical processing order (hsa, mmu, rno, mml; unknown codes last, sorted)."""
    rank = {code: i for i, code in enumerate(ORGANISM_ORDER)}
    return sorted(dict.fromkeys(organisms), key=lambda c: (rank.get(c, len(rank)), c))


def config_taxids() -> dict[str, int]:
    """``code -> canonical_tax_id`` from the repo ``config.yaml`` organisms."""
    return {str(o["code"]): int(o["canonical_tax_id"]) for o in load_config().organisms}


def run_diana_method_audit(
    *,
    cache_dir: Path,
    manifest_path: Path | None,
    organisms: Sequence[str],
    out_dir: Path,
    provenance_dir: Path | None,
    report_only: bool,
    taxids: Mapping[str, int] | None = None,
    report_path: str | None = None,
    log: Callable[[str], None] | None = None,
) -> AuditResult:
    """Load the cached DIANA files, run the census, write the artefacts, gate.

    Flow: (1) if ``manifest_path`` exists, sha256-verify the DIANA files present
    under ``cache_dir`` (:func:`verify_cached_sources`; a mismatch raises
    ``RuntimeError`` -- a gate failure the CLI turns into exit 1); a missing DIANA
    file raises ``FileNotFoundError``. (2) For each organism, in canonical order,
    load the rows with the production loader ``tarbase.load_mti_for_taxid`` (its
    vocabulary ``UserWarning``\\ s are captured into ``summary["loader_warnings"]``;
    warnings of any other category are re-emitted, not recorded, so an
    interpreter-specific ``DeprecationWarning`` cannot change the deterministic
    artefacts) and tally them. The loader's own fail-loud ``ValueError``\\ s
    (header drift, 0 strong rows, an empty species slice) propagate. (3) Finalize,
    summarize, write the four artefacts under ``out_dir``. (4) Write the
    provenance block if ``provenance_dir`` is given.

    Args:
        cache_dir: ``paths.source_cache_dir`` holding the cached v9 downloads.
        manifest_path: ``provenance/manifest.yaml`` (skipped when absent / None).
        organisms: Organism codes to audit (e.g. ``["hsa", "mmu", "rno", "mml"]``).
        out_dir: Where the four artefacts are written.
        provenance_dir: Where ``diana_method_audit.json`` (provenance) is written;
            ``None`` skips it. Must differ from ``out_dir`` (the audit JSON and
            the provenance block share a basename).
        report_only: Diagnostic mode -- exit code 0 regardless, with
            ``summary["report_only"] = True`` and ``strict_gate_would_fail`` set.
        taxids: ``code -> canonical NCBI taxid``; default: the repo ``config.yaml``
            organisms (the CLI passes the ``-c`` config's).
        report_path: The repo-relative label of the TSV recorded in the
            provenance block. Default: the TSV's absolute posix path -- stable
            per host but NOT byte-identical across hosts, the one provenance
            field that may differ for identical inputs when the report lives
            outside the repo (the CLI passes a repo-relative label whenever the
            output is under the repo).
        log: Progress sink (one line per organism); ``None`` for silence.

    Returns:
        :class:`AuditResult`.
    """
    cache_dir = Path(cache_dir)
    out_dir = Path(out_dir)
    if provenance_dir is not None and Path(provenance_dir).resolve() == out_dir.resolve():
        raise ValueError(
            f"diana method audit: provenance_dir must differ from out_dir ({out_dir}); "
            f"both write a {PROVENANCE_JSON} and the provenance block would overwrite "
            f"the full audit JSON"
        )
    if taxids is None:
        taxids = config_taxids()
    codes = _ordered(organisms)
    if not codes:
        raise ValueError("diana method audit: no organisms requested")
    missing_taxid = [c for c in codes if c not in taxids]
    if missing_taxid:
        raise ValueError(
            f"diana method audit: no canonical_tax_id for organism(s) {missing_taxid}"
        )
    emit = log or (lambda _msg: None)

    # (1) Provenance: sha256 gate on the DIANA files actually used, then identity.
    files: dict[str, list[str]] = defaultdict(list)  # basename -> organisms
    for code in codes:
        files[tarbase.filename_for_taxid(int(taxids[code]))].append(code)
    pinned: dict[str, str] = {}
    if manifest_path is not None and Path(manifest_path).exists():
        verify_cached_sources(cache_dir, Path(manifest_path), names=sorted(files))
        pinned = pinned_sha256_index(Path(manifest_path))
    source_files: list[dict[str, Any]] = []
    for basename in sorted(files):
        path = cache_dir / basename
        if not path.exists():
            raise FileNotFoundError(
                f"DIANA-TarBase v9 file not cached: {path} (needed for "
                f"{', '.join(sorted(files[basename]))})"
            )
        source_files.append(
            {
                "organisms": sorted(files[basename]),
                "path": basename,
                "sha256": file_sha256(path),
                "pinned_sha256": pinned.get(basename),
                "bytes": path.stat().st_size,
            }
        )

    # (2) Load with the production loaders; capture (do not silence) their warnings.
    acc = MethodAccumulator()
    loader_warnings: dict[str, list[str]] = {}
    for code in codes:
        taxid = int(taxids[code])
        path = cache_dir / tarbase.filename_for_taxid(taxid)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            rows = tarbase.load_mti_for_taxid(path, taxid=taxid)
        # Only the loader's vocabulary UserWarnings belong in the deterministic
        # artefacts; anything else (Deprecation/Resource warnings that vary by
        # interpreter or library version) is re-emitted at its origin instead.
        loader_warnings[code] = [
            str(w.message) for w in caught if issubclass(w.category, UserWarning)
        ]
        for w in caught:
            if not issubclass(w.category, UserWarning):
                warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
        n_before = len(acc.methods())
        n_rows = audit_rows(rows, organism=code, acc=acc)
        del rows
        emit(
            f"[diana-audit] {code}: rows={n_rows} methods={len(acc.methods_for(code))} "
            f"(total={len(acc.methods())}, new={len(acc.methods()) - n_before}) "
            f"file={path.name} warnings={len(loader_warnings[code])}"
        )

    # (3) Finalize + summarize + write.
    records = finalize(acc)
    summary = summarize(
        records,
        source_files=source_files,
        n_pairs_high_from_unknown=acc.n_unknown_high_pairs(),
    )
    summary["loader_warnings"] = loader_warnings
    summary["report_only"] = bool(report_only)
    summary["strict_gate_would_fail"] = summary["gate_status"] == GATE_FAIL

    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = out_dir / AUDIT_TSV
    write_audit_tsv(records, tsv_path)
    write_audit_json(records, summary, out_dir / AUDIT_JSON)
    write_summary_json(summary, out_dir / SUMMARY_JSON)
    write_audit_md(records, summary, out_dir / AUDIT_MD)
    artefacts: dict[str, Path] = {
        "tsv": tsv_path,
        "json": out_dir / AUDIT_JSON,
        "summary": out_dir / SUMMARY_JSON,
        "md": out_dir / AUDIT_MD,
    }

    # (4) Provenance block (pins the report + sources next to the run manifests).
    if provenance_dir is not None:
        block = provenance_block(
            summary,
            report_sha256=file_sha256(tsv_path),
            report_path=report_path or tsv_path.resolve().as_posix(),
        )
        prov_path = Path(provenance_dir) / PROVENANCE_JSON
        _dump_json(block, prov_path)
        artefacts["provenance"] = prov_path

    # (5) Exit code: strict -> fail is 1; report-only never gates.
    exit_code = 0 if report_only or summary["gate_status"] != GATE_FAIL else 1
    return AuditResult(records=records, summary=summary, artefacts=artefacts, exit_code=exit_code)


def verdict_lines(result: AuditResult) -> list[str]:
    """Human-readable summary + verdict lines the CLI prints after a run."""
    s = result.summary
    lines = [
        f"[diana-audit] source: {s['source_release']}; organisms: {', '.join(s['organisms'])}",
        f"[diana-audit] raw method strings={s['n_raw_method_strings']} "
        f"normalized={s['n_normalized_methods']} known_high={s['n_known_high']} "
        f"known_clip={s['n_known_clip']} known_other={s['n_known_other']} "
        f"unknown={s['n_unknown']} empty_method={s['n_empty_method']}",
        f"[diana-audit] unknown contributing HIGH: {s['n_unknown_contributing_high']} "
        f"method(s), rows={s['n_rows_high_from_unknown']} "
        f"pairs={_fmt_count(s['n_pairs_high_from_unknown'])} "
        f"organisms={', '.join(s['organisms_with_unknown_high']) or 'none'}",
    ]
    if s["unknown_methods"]:
        lines.append(f"[diana-audit] unknown methods: {', '.join(s['unknown_methods'])}")
    for name, path in sorted(result.artefacts.items()):
        lines.append(f"[diana-audit] wrote {name}: {path}")
    status = s["gate_status"]
    if status == GATE_FAIL and s.get("report_only"):
        verdict = "WOULD FAIL (report-only, NOT release-gating)"
    else:
        verdict = status.upper()
    lines.append(f"diana method gate: {verdict}")
    by_method = {r.method: r for r in result.records}
    for method in s["unknown_high_methods"]:
        r = by_method[method]
        lines.append(
            f"  unknown method contributes HIGH: {method!r} "
            f"n_rows_high={r.n_rows_high} n_pairs_high={r.n_pairs_high} "
            f"organisms={','.join(r.organisms)}"
        )
    return lines
