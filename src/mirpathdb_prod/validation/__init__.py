"""Validation layer.

Produces EXPLICIT pass/fail results. Validators consume canonical models and
emit structured outcomes (no side effects beyond returning results / writing to
the configured ``validation/`` dir via reporting). The headline validator is
organism concordance, under the "no-loss FLOOR + expand" model:

  * Floor preserved: the miRPathDB2 floor (human 9606 + mouse 10090) MUST be
    present in production. A missing floor organism is a FAIL.
  * Approved additions: production may expand beyond the floor only with
    ``approved_extensions`` (rat 10116 — native validated + predicted; macaque
    9544 — predicted-dominant: TargetScan LOW, with a lone DIANA-TarBase v9
    validated MTI marked). An unapproved addition is
    a FAIL. Dog (9615) + chicken (9031) were evaluated and EXCLUDED (owner-omitted,
    even though TargetScan + miRDB cover them), so they are in neither allowlist.
  * Approved removals: an old-prototype organism dropped from production must be
    an ``approved_removals`` entry (a DOCUMENTED removal, not an accidental loss).
    The allowlist is EMPTY now: the prototype carried mml, and production re-adds
    it as an approved extension, so nothing is removed.
"""
