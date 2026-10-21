"""Organism parity / concordance validation — "no-loss FLOOR + expand" model.

Compares organism sets by ``canonical_tax_id`` first, then scientific name +
curated synonyms — never the display label alone. Emits an explicit pass/fail
result.

The production organism set is 4 species: the miRPathDB2 floor (human 9606,
mouse 10090) plus two approved native-coverage additions — rat (10116; native
validated + predicted) and macaque (Macaca mulatta, 9544; predicted-dominant —
TargetScan LOW, with a lone DIANA-TarBase v9 validated MTI marked). Dog (cfa,
9615) and chicken (gga, 9031) were evaluated and
EXCLUDED (owner 2026-06-24 — an owner choice to omit, even though TargetScan +
miRDB cover them); they appear in neither allowlist. There are no approved
removals: the prototype carried mml, and production re-adds it.

Acceptance policy ("no-loss FLOOR + expand"):
  * ``hard_floor`` taxids MUST be present in production. Missing one => FAIL.
  * A production organism that is neither in ``hard_floor`` nor in
    ``approved_extensions`` => FAIL (unapproved addition).
  * An old-prototype organism missing from production AND not in
    ``approved_removals`` => FAIL (organism loss).
  * The miRPathDB2 reference is the primary external anchor for the floor; any
    floor organism that miRPathDB2 carries but production drops is also caught
    by the floor check.

``compare_organism_sets`` remains the generic taxid-first two-set comparison
used for audit cross-checks (and the floor-vs-production legs below).
"""

from __future__ import annotations

from collections.abc import Iterable

from mirpathdb_prod.models import Organism
from mirpathdb_prod.normalize.taxonomy import normalize_name


def _index_by_taxid(organisms: Iterable[Organism]) -> dict[int, Organism]:
    """Index organisms by canonical taxid (last record wins on collision)."""
    return {o.canonical_tax_id: o for o in organisms}


def _org_brief(org: Organism) -> dict:
    """A compact, JSON-serializable view of an organism (audit fields preserved)."""
    return {
        "canonical_tax_id": org.canonical_tax_id,
        "scientific_name": org.scientific_name,
        "source_name": org.source_name,
        "common_name": org.common_name,
        "source": org.source,
    }


def compare_organism_sets(
    expected: Iterable[Organism],
    actual: Iterable[Organism],
    *,
    approved_additions: Iterable[int] = (),
) -> dict:
    """Compare two organism sets by canonical taxid and report parity.

    This is the generic two-set comparison. It is reused by
    :func:`evaluate_organism_parity` for audit cross-checks (e.g. confirming the
    production set contains the miRPathDB2 floor).

    Args:
        expected: The reference organism set.
        actual: The organism set under test.
        approved_additions: Canonical taxids allowed to be present in ``actual``
            but absent from ``expected`` without failing.

    Returns:
        A structured, JSON-serializable result mapping with an explicit
        ``passed`` flag.
    """
    approved = set(approved_additions)

    expected_by_tax = _index_by_taxid(expected)
    actual_by_tax = _index_by_taxid(actual)
    expected_taxids = set(expected_by_tax)
    actual_taxids = set(actual_by_tax)

    # --- comparison by taxid (the authoritative key) ---
    missing_taxids = sorted(expected_taxids - actual_taxids)  # in expected, not in actual
    extra_taxids = sorted(actual_taxids - expected_taxids)  # in actual, not in expected
    common_taxids = sorted(expected_taxids & actual_taxids)

    # Approved vs unapproved additions (extras).
    approved_extras = sorted(t for t in extra_taxids if t in approved)
    unapproved_extras = sorted(t for t in extra_taxids if t not in approved)

    # --- comparison by normalized scientific name (audit / synonym cross-check) ---
    # Detect "synonym-only" differences: same taxid, but the scientific names
    # disagree even after normalization (would indicate a labelling drift).
    synonym_only_diffs = []
    for t in common_taxids:
        e_name = expected_by_tax[t].scientific_name
        a_name = actual_by_tax[t].scientific_name
        if normalize_name(e_name) != normalize_name(a_name):
            synonym_only_diffs.append(
                {
                    "canonical_tax_id": t,
                    "expected_scientific_name": e_name,
                    "actual_scientific_name": a_name,
                }
            )

    # --- ambiguous / unresolved mappings (sentinel taxid 0) ---
    ambiguous = [
        _org_brief(o)
        for o in list(expected) + list(actual)
        if o.canonical_tax_id == 0
    ]
    # de-dup while preserving determinism
    seen = set()
    ambiguous_unique = []
    for a in ambiguous:
        key = (a["source"], a["source_name"])
        if key not in seen:
            seen.add(key)
            ambiguous_unique.append(a)

    # --- pass/fail ---
    failures: list[str] = []
    if missing_taxids:
        failures.append(
            "organism(s) present in expected but missing from actual: "
            + ", ".join(str(t) for t in missing_taxids)
        )
    if unapproved_extras:
        failures.append(
            "organism(s) present in actual but absent from expected and not approved: "
            + ", ".join(str(t) for t in unapproved_extras)
        )
    if ambiguous_unique:
        failures.append(
            f"{len(ambiguous_unique)} organism record(s) could not be resolved to a canonical taxid"
        )

    passed = not failures

    return {
        "passed": passed,
        "expected_taxids": sorted(expected_taxids),
        "actual_taxids": sorted(actual_taxids),
        "common_taxids": common_taxids,
        "missing_taxids": missing_taxids,
        "extra_taxids": extra_taxids,
        "approved_additions": sorted(approved),
        "approved_extras": approved_extras,
        "unapproved_extras": unapproved_extras,
        "synonym_only_diffs": synonym_only_diffs,
        "ambiguous_mappings": ambiguous_unique,
        "failures": failures,
        "expected": [_org_brief(o) for o in expected_by_tax.values()],
        "actual": [_org_brief(o) for o in actual_by_tax.values()],
    }


def _classify_production(
    org: Organism,
    *,
    hard_floor: set[int],
    approved_extensions: set[int],
) -> dict:
    """Classify one production organism against the floor + expand policy."""
    t = org.canonical_tax_id
    if t == 0:
        klass = "UNRESOLVED"
        ok = False
    elif t in hard_floor:
        klass = "floor"
        ok = True
    elif t in approved_extensions:
        klass = "approved_extension"
        ok = True
    else:
        klass = "unapproved_addition"
        ok = False
    return {
        "canonical_tax_id": t,
        "scientific_name": org.scientific_name,
        "source_name": org.source_name,
        "common_name": org.common_name,
        "classification": klass,
        "ok": ok,
    }


def evaluate_organism_parity(
    old: Iterable[Organism],
    new: Iterable[Organism],
    mirpathdb2: Iterable[Organism],
    *,
    hard_floor: Iterable[int] = (),
    approved_extensions: Iterable[int] = (),
    approved_removals: Iterable[int] = (),
) -> dict:
    """Run the organism-concordance gate under the "no-loss FLOOR + expand" model.

    Three policy checks, all config-driven:

    1. **Floor preserved** — every ``hard_floor`` taxid (the miRPathDB2 floor:
       human + mouse) MUST be present in ``new``. A missing floor organism is a
       FAIL.
    2. **Additions approved** — every ``new`` organism that is not in the floor
       must be in ``approved_extensions`` (each justified by the organism
       availability matrix). An unapproved addition is a FAIL.
    3. **Removals approved** — every ``old`` prototype organism that is absent
       from ``new`` must be in ``approved_removals`` (a documented removal). An
       undocumented loss is a FAIL.

    miRPathDB2 stays the primary external anchor: the floor is exactly the
    miRPathDB2 set, and the report cross-checks that the floor is a subset of
    both miRPathDB2 and ``new``.

    Args:
        old: Old prototype organism set (historical baseline for the removal
            audit; hsa/mmu/mml).
        new: New production organism set (under test).
        mirpathdb2: miRPathDB2 reference set (primary external anchor; hsa/mmu).
        hard_floor: Canonical taxids that MUST be present in production.
        approved_extensions: Canonical taxids allowed beyond the floor.
        approved_removals: Old-prototype taxids deliberately not carried forward.

    Returns:
        A structured, JSON-serializable result with an overall ``passed`` flag,
        per-organism classifications, the documented additions/removals, and a
        human-readable ``summary``.
    """
    old = list(old)
    new = list(new)
    mirpathdb2 = list(mirpathdb2)
    floor = set(hard_floor)
    extensions = set(approved_extensions)
    removals = set(approved_removals)

    old_taxids = {o.canonical_tax_id for o in old}
    new_taxids = {o.canonical_tax_id for o in new}
    v2_taxids = {o.canonical_tax_id for o in mirpathdb2}

    new_by_tax = _index_by_taxid(new)

    # --- per-organism classification of the production set ---
    classifications = [
        _classify_production(
            o, hard_floor=floor, approved_extensions=extensions
        )
        for o in sorted(new, key=lambda x: x.canonical_tax_id)
    ]

    # --- 1) floor preserved ---
    floor_missing = sorted(floor - new_taxids)
    floor_present = sorted(floor & new_taxids)

    # --- 2) additions approved ---
    additions = sorted(new_taxids - floor)  # production organisms beyond the floor
    approved_addition_taxids = sorted(t for t in additions if t in extensions)
    unapproved_additions = sorted(t for t in additions if t not in extensions)

    # --- 3) removals approved ---
    dropped = sorted(old_taxids - new_taxids)  # old-prototype organisms not in production
    approved_removal_taxids = sorted(t for t in dropped if t in removals)
    unapproved_removals = sorted(t for t in dropped if t not in removals)

    # --- unresolved records in production (sentinel taxid 0) ---
    unresolved = [c for c in classifications if c["classification"] == "UNRESOLVED"]

    # --- detail blocks (audit-friendly, name-resolved) ---
    def _detail(taxids: list[int]) -> list[dict]:
        out = []
        for t in taxids:
            org = new_by_tax.get(t)
            old_org = next((o for o in old if o.canonical_tax_id == t), None)
            out.append(
                {
                    "canonical_tax_id": t,
                    "scientific_name": (org or old_org).scientific_name
                    if (org or old_org)
                    else "",
                    "source_name": (org or old_org).source_name
                    if (org or old_org)
                    else "",
                    "in_floor": t in floor,
                    "in_mirpathdb2": t in v2_taxids,
                    "in_old_prototype": t in old_taxids,
                    "in_production": t in new_taxids,
                }
            )
        return out

    approved_additions_detail = _detail(approved_addition_taxids)
    approved_removals_detail = _detail(approved_removal_taxids)

    # --- floor cross-check against the miRPathDB2 anchor ---
    floor_not_in_v2 = sorted(floor - v2_taxids)  # floor should be a subset of miRPathDB2

    # --- pass/fail ---
    failures: list[str] = []
    if floor_missing:
        failures.append(
            "hard_floor organism(s) missing from production: "
            + ", ".join(str(t) for t in floor_missing)
        )
    if unapproved_additions:
        failures.append(
            "production organism(s) neither floor nor approved_extension: "
            + ", ".join(str(t) for t in unapproved_additions)
        )
    if unapproved_removals:
        failures.append(
            "old-prototype organism(s) dropped from production without an approved removal: "
            + ", ".join(str(t) for t in unapproved_removals)
        )
    if unresolved:
        failures.append(
            f"{len(unresolved)} production organism record(s) could not be resolved "
            "to a canonical taxid"
        )
    if floor_not_in_v2:
        # The hard_floor must be a subset of the miRPathDB2 anchor (the floor is
        # DEFINED as the miRPathDB2-covered species). A floor taxid absent from v2
        # means the floor/anchor are misconfigured -- previously computed and
        # reported but never gated; now a hard failure.
        failures.append(
            "hard_floor organism(s) not present in the miRPathDB2 anchor "
            "(floor must be a subset of miRPathDB2): "
            + ", ".join(str(t) for t in floor_not_in_v2)
        )

    overall_passed = not failures

    summary_lines: list[str] = []
    if overall_passed:
        summary_lines.append(
            "Organism concordance PASS (no-loss FLOOR + expand): "
            f"floor preserved ({', '.join(str(t) for t in floor_present)}), "
            f"{len(approved_addition_taxids)} approved addition(s), "
            f"{len(approved_removal_taxids)} approved removal(s)."
        )
        summary_lines.append(
            "Floor preserved: " + ", ".join(str(t) for t in floor_present) + "."
        )
        for d in approved_additions_detail:
            summary_lines.append(
                f"Approved addition: taxid {d['canonical_tax_id']} "
                f"({d['scientific_name']}) — native target coverage beyond the "
                "miRPathDB2 floor."
            )
        for d in approved_removals_detail:
            summary_lines.append(
                f"Approved removal: taxid {d['canonical_tax_id']} "
                f"({d['scientific_name']}) — documented removal, not an accidental loss."
            )
    else:
        summary_lines.append("Organism concordance FAIL (no-loss FLOOR + expand).")
        summary_lines.extend(failures)

    return {
        "passed": overall_passed,
        "verdict": "PASS" if overall_passed else "FAIL",
        "model": "floor+expand",
        "hard_floor": sorted(floor),
        "approved_extensions": sorted(extensions),
        "approved_removals": sorted(removals),
        # floor check
        "floor_present": floor_present,
        "floor_missing": floor_missing,
        "floor_not_in_mirpathdb2": floor_not_in_v2,
        # additions
        "production_additions": additions,
        "approved_additions": approved_addition_taxids,
        "approved_additions_detail": approved_additions_detail,
        "unapproved_additions": unapproved_additions,
        # removals
        "prototype_dropped": dropped,
        "approved_removals_applied": approved_removal_taxids,
        "approved_removals_detail": approved_removals_detail,
        "unapproved_removals": unapproved_removals,
        # production set
        "production_classifications": classifications,
        "unresolved": unresolved,
        "failures": failures,
        "summary": summary_lines,
        "sizes": {
            "old": len(old),
            "new": len(new),
            "mirpathdb2": len(mirpathdb2),
        },
        "taxids": {
            "old": sorted(old_taxids),
            "new": sorted(new_taxids),
            "mirpathdb2": sorted(v2_taxids),
        },
    }
