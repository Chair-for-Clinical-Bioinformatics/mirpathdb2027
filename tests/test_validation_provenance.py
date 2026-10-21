"""Tests for the manifest-checksum validator (the reconciliation gate's verify half)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

from mirpathdb_prod.validation.provenance import file_sha256, validate_manifest


def _write(path: Path, data: bytes) -> str:
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def _manifest(tmp: Path, raw_files: list[dict]) -> Path:
    doc = {"sources": {"demo": {"role": "production_data_source", "raw_files": raw_files}}}
    mpath = tmp / "manifest.yaml"
    mpath.write_text(yaml.safe_dump(doc), encoding="utf-8")
    return mpath


def test_file_sha256_matches_hashlib(tmp_path: Path) -> None:
    f = tmp_path / "x.bin"
    expected = _write(f, b"hello world\n" * 1000)
    assert file_sha256(f) == expected


def test_pinned_file_verifies(tmp_path: Path) -> None:
    f = tmp_path / "src.gz"
    digest = _write(f, b"real source bytes")
    mpath = _manifest(tmp_path, [{"path": str(f), "bytes": f.stat().st_size, "sha256": digest}])
    res = validate_manifest(mpath, tmp_path)
    assert res["ok"] is True
    assert len(res["checked"]) == 1 and not res["failures"]


def test_sha256_mismatch_is_failure(tmp_path: Path) -> None:
    f = tmp_path / "src.gz"
    _write(f, b"real source bytes")
    mpath = _manifest(tmp_path, [{"path": str(f), "sha256": "deadbeef"}])
    res = validate_manifest(mpath, tmp_path)
    assert res["ok"] is False
    assert res["failures"][0]["error"] == "sha256_mismatch"


def test_missing_pinned_file_is_failure(tmp_path: Path) -> None:
    mpath = _manifest(tmp_path, [{"path": str(tmp_path / "absent.gz"), "sha256": "abc"}])
    res = validate_manifest(mpath, tmp_path)
    assert res["ok"] is False
    assert res["failures"][0]["error"] == "missing"


def test_on_demand_file_recorded_not_failed(tmp_path: Path) -> None:
    f = tmp_path / "big.gz"
    digest = _write(f, b"unpinned but present")
    mpath = _manifest(tmp_path, [{"path": str(f), "bytes": None, "checksum": "on-demand"}])
    res = validate_manifest(mpath, tmp_path)
    assert res["ok"] is True  # unpinned-present is not a failure
    rec = next(s for s in res["skipped"] if s["path"] == str(f))
    assert rec["sha256"] == digest and rec["bytes"] == f.stat().st_size
