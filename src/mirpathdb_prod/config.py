"""Configuration loading for the production build.

Loads ``config.yaml`` (with optional ``config.local.yaml`` override) and exposes
configurable paths and the production organism list. No filesystem paths are
hard-coded in source — they all come from the config file so the build is
portable across hosts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# Repository root = two levels up from this file
# (src/mirpathdb_prod/config.py -> src -> repo root).
REPO_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_CONFIG_NAME = "config.yaml"
LOCAL_CONFIG_NAME = "config.local.yaml"


@dataclass(frozen=True)
class Config:
    """Parsed build configuration.

    Attributes:
        paths: Mapping of logical path names to resolved absolute ``Path``s.
        organisms: The production organism list (raw dicts from ``config.yaml``).
        raw: The full parsed config document, for keys not promoted above.
        source_path: The config file that was loaded.
    """

    paths: dict[str, Path]
    organisms: list[dict[str, Any]]
    raw: dict[str, Any] = field(default_factory=dict)
    source_path: Path | None = None

    def path(self, name: str) -> Path:
        """Return a configured path by logical name (raises KeyError if absent)."""
        return self.paths[name]


def _resolve(value: str, base: Path) -> Path:
    """Resolve a configured path against the repo root if it is relative."""
    p = Path(value).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` into ``base`` (override wins)."""
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(config_path: str | Path | None = None) -> Config:
    """Load ``config.yaml`` (plus optional local override) into a :class:`Config`.

    Args:
        config_path: Explicit path to a config file. If ``None``, uses
            ``<repo_root>/config.yaml`` and, when present,
            ``<repo_root>/config.local.yaml`` as an override.

    Returns:
        A frozen :class:`Config` with paths resolved relative to the repo root.
    """
    if config_path is not None:
        primary = Path(config_path).expanduser().resolve()
        base = primary.parent
    else:
        primary = REPO_ROOT / DEFAULT_CONFIG_NAME
        base = REPO_ROOT

    with primary.open("r", encoding="utf-8") as fh:
        doc: dict[str, Any] = yaml.safe_load(fh) or {}

    local = base / LOCAL_CONFIG_NAME
    if local.exists():
        with local.open("r", encoding="utf-8") as fh:
            doc = _deep_merge(doc, yaml.safe_load(fh) or {})

    raw_paths: dict[str, Any] = doc.get("paths", {}) or {}
    resolved = {name: _resolve(str(value), base) for name, value in raw_paths.items()}

    organisms: list[dict[str, Any]] = list(doc.get("organisms", []) or [])

    return Config(
        paths=resolved,
        organisms=organisms,
        raw=doc,
        source_path=primary,
    )
