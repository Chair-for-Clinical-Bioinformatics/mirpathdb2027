"""Markdown report rendering (deterministic)."""

from __future__ import annotations

import io
from typing import Any


def _taxid_list(values: list[int]) -> str:
    return ", ".join(str(v) for v in values) if values else "(none)"


def _detail_table(rows: list[dict[str, Any]]) -> list[str]:
    """Render an additions/removals detail block as a Markdown table."""
    out = [
        "| taxid | scientific_name | source_name | in_floor | in_mirpathdb2 "
        "| in_old_prototype | in_production |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in sorted(rows, key=lambda d: d.get("canonical_tax_id", 0)):
        out.append(
            f"| {r.get('canonical_tax_id', '')} | {r.get('scientific_name', '')} "
            f"| {r.get('source_name', '')} | {r.get('in_floor', '')} "
            f"| {r.get('in_mirpathdb2', '')} | {r.get('in_old_prototype', '')} "
            f"| {r.get('in_production', '')} |"
        )
    if not rows:
        out.append("| (none) | — | — | — | — | — | — |")
    return out


def render_parity_report(result: dict) -> str:
    """Render an organism-parity validation result as a Markdown report.

    Accepts the structured result from
    :func:`mirpathdb_prod.validation.organism_parity.evaluate_organism_parity`
    (the "no-loss FLOOR + expand" gate).

    Args:
        result: The structured validation result.

    Returns:
        A deterministic Markdown document.
    """
    buf = io.StringIO()
    w = buf.write

    verdict = result.get("verdict", "PASS" if result.get("passed") else "FAIL")
    w("# Organism Parity / Concordance Report\n\n")
    w(f"**Overall verdict: {verdict}**\n\n")
    w(
        "Gate model: **no-loss FLOOR + expand**. Comparison is by canonical NCBI "
        "taxid first, then normalized scientific name + curated synonyms — never "
        "the display label alone. The miRPathDB2 set is the primary external "
        "anchor for the floor.\n\n"
    )

    for line in result.get("summary", []):
        w(f"- {line}\n")
    w("\n")

    sizes = result.get("sizes", {})
    taxids = result.get("taxids", {})
    w("## Set sizes & taxids\n\n")
    w(
        f"- old (prototype baseline): {sizes.get('old', '?')} organism(s) "
        f"— {_taxid_list(taxids.get('old', []))}\n"
    )
    w(
        f"- new (production): {sizes.get('new', '?')} organism(s) "
        f"— {_taxid_list(taxids.get('new', []))}\n"
    )
    w(
        f"- miRPathDB2 (external anchor): {sizes.get('mirpathdb2', '?')} organism(s) "
        f"— {_taxid_list(taxids.get('mirpathdb2', []))}\n\n"
    )

    w("## Gate policy (config-driven)\n\n")
    w(f"- hard_floor (must be present): {_taxid_list(result.get('hard_floor', []))}\n")
    w(
        "- approved_extensions (allowed additions): "
        f"{_taxid_list(result.get('approved_extensions', []))}\n"
    )
    w(
        "- approved_removals (documented removals): "
        f"{_taxid_list(result.get('approved_removals', []))}\n\n"
    )

    # --- 1) floor preserved ---
    w("## 1. Floor preserved (no organism loss vs. miRPathDB2 floor)\n\n")
    floor_missing = result.get("floor_missing", [])
    w(f"- result: **{'PASS' if not floor_missing else 'FAIL'}**\n")
    w(f"- floor present in production: {_taxid_list(result.get('floor_present', []))}\n")
    w(f"- floor MISSING from production (FAIL if any): {_taxid_list(floor_missing)}\n")
    floor_not_v2 = result.get("floor_not_in_mirpathdb2", [])
    if floor_not_v2:
        w(
            "- note: floor taxid(s) not found in the miRPathDB2 anchor "
            f"(allowlist drift): {_taxid_list(floor_not_v2)}\n"
        )
    w("\n")

    # --- 2) approved additions ---
    w("## 2. Approved additions (expansion beyond the floor)\n\n")
    unapproved_add = result.get("unapproved_additions", [])
    w(f"- result: **{'PASS' if not unapproved_add else 'FAIL'}**\n")
    additions = _taxid_list(result.get("production_additions", []))
    w(f"- production additions beyond floor: {additions}\n")
    w(f"- approved additions: {_taxid_list(result.get('approved_additions', []))}\n")
    w(f"- UNAPPROVED additions (FAIL if any): {_taxid_list(unapproved_add)}\n\n")
    for line in _detail_table(result.get("approved_additions_detail", [])):
        w(line + "\n")
    w("\n")

    # --- 3) approved removals ---
    w("## 3. Approved removals (documented, not accidental loss)\n\n")
    unapproved_rm = result.get("unapproved_removals", [])
    w(f"- result: **{'PASS' if not unapproved_rm else 'FAIL'}**\n")
    w(f"- prototype organisms dropped: {_taxid_list(result.get('prototype_dropped', []))}\n")
    w(f"- approved removals applied: {_taxid_list(result.get('approved_removals_applied', []))}\n")
    w(f"- UNAPPROVED removals (FAIL if any): {_taxid_list(unapproved_rm)}\n\n")
    for line in _detail_table(result.get("approved_removals_detail", [])):
        w(line + "\n")
    w("\n")

    # --- production set classification ---
    w("## Production organism set (classification)\n\n")
    w("| taxid | scientific_name | source_name | common_name | classification | ok |\n")
    w("|---|---|---|---|---|---|\n")
    for c in sorted(
        result.get("production_classifications", []),
        key=lambda d: d.get("canonical_tax_id", 0),
    ):
        w(
            f"| {c.get('canonical_tax_id', '')} | {c.get('scientific_name', '')} "
            f"| {c.get('source_name', '')} | {c.get('common_name', '')} "
            f"| {c.get('classification', '')} | {c.get('ok', '')} |\n"
        )
    w("\n")

    # --- unresolved mappings ---
    unresolved = result.get("unresolved", [])
    w("### Ambiguous / unresolved mappings\n\n")
    if unresolved:
        for a in unresolved:
            w(f"- {a.get('source_name')} (no canonical taxid)\n")
    else:
        w("(none — every production record resolved to a canonical taxid)\n")
    w("\n")

    if result.get("failures"):
        w("## Failures\n\n")
        for f in result["failures"]:
            w(f"- {f}\n")
        w("\n")

    return buf.getvalue()
