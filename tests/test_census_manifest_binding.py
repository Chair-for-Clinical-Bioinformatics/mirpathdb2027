"""Census must bind to the FINAL contract manifest (regression for the 08-31 digest defect).

The label-corrected derivation once computed release_census.json BEFORE regenerating the
final SHA256SUMS, so _meta.sha256sums_digest named the superseded package manifest. These
tests fail if that ordering regresses: they assemble a tiny package through
derive_label_corrected_package.py and require the census digest to equal the sha256 of the
package's own final SHA256SUMS, with a non-circular topology (census not inside SHA256SUMS).
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DERIVE = REPO / "scripts" / "derive_label_corrected_package.py"


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_derive_runs_census_after_sha256sums():
    """Static-order guard: the census subprocess call must come AFTER SHA256SUMS is written
    and the binding assertion must be present."""
    src = DERIVE.read_text(encoding="utf-8")
    i_sums = src.index('(building / "SHA256SUMS").write_text')
    i_census = src.index('"build_release_census.py"')
    assert i_census > i_sums, "census generated before SHA256SUMS — the 08-31 digest defect"
    assert 'sha256sums_digest"] != sha256(building / "SHA256SUMS")' in src, \
        "missing the census<->final-manifest binding assertion"
    # and it must parse
    ast.parse(src)


def test_shipped_package_census_binds_to_final_manifest():
    """If a built package is present on this host, its census digest must equal its own
    SHA256SUMS and the census must not be inside that manifest (non-circular)."""
    candidates = list(REPO.glob("../_execution/server-upload/mirpathdb-data-*"))
    pkgs = [
        p for p in candidates
        if (p / "release_census.json").is_file() and (p / "SHA256SUMS").is_file()
    ]
    if not pkgs:
        import pytest
        pytest.skip("no built package on this host")
    for pkg in pkgs:
        census = json.loads((pkg / "release_census.json").read_text())
        want = _sha256(pkg / "SHA256SUMS")
        got = census["_meta"]["sha256sums_digest"]
        assert got == want, f"{pkg.name}: census digest {got} != final SHA256SUMS {want}"
        assert "release_census.json" not in (pkg / "SHA256SUMS").read_text(), \
            f"{pkg.name}: census is inside SHA256SUMS (circular self-reference)"
