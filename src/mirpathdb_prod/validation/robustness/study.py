"""Study configuration, scope (strata + query sizes), manifest and stage orchestration."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from mirpathdb_prod.validation.robustness import STUDY_MODULE_VERSION
from mirpathdb_prod.validation.robustness.inputs import (
    ReleaseInputs,
    Stratum,
    load_manifest_cell_counts,
    percentile_sizes,
)

REACTOME_PREFIX = {"hsa": "R-HSA-", "mmu": "R-MMU-", "rno": "R-RNO-", "mml": "R-MML-"}


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
@dataclass
class StudyConfig:
    study_id: str
    seed: int
    replicates: int
    alphas: list[float]
    retention_alpha: float
    organisms: list[str]
    families: list[str]
    tiers: list[str]
    query_size_percentiles: list[float]
    cohort_min_size: int
    figure_examples: dict[str, dict[str, object]]
    figure_terms: dict[str, list[str]]
    degree_null: dict[str, object]
    r_kernel: dict[str, object]
    raw: dict[str, object] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: str | Path) -> StudyConfig:
        raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
        return cls.from_dict(raw)

    @classmethod
    def from_dict(cls, raw: dict[str, object]) -> StudyConfig:
        def req(key: str) -> object:
            if key not in raw:
                raise ValueError(f"study config: missing required key {key!r}")
            return raw[key]

        return cls(
            study_id=str(req("study_id")),
            seed=int(req("seed")),
            replicates=int(req("replicates")),
            alphas=[float(a) for a in req("alphas")],
            retention_alpha=float(raw.get("retention_alpha", 0.25)),
            organisms=[str(o) for o in req("organisms")],
            families=[str(f) for f in req("families")],
            tiers=[str(t) for t in req("tiers")],
            query_size_percentiles=[
                float(p) for p in raw.get("query_size_percentiles", [10, 50, 90])
            ],
            cohort_min_size=int(raw.get("cohort_min_size", 30)),
            figure_examples=dict(raw.get("figure_examples", {}) or {}),
            figure_terms={
                str(k): [str(t) for t in v] for k, v in (raw.get("figure_terms", {}) or {}).items()
            },
            degree_null=dict(raw.get("degree_null", {}) or {}),
            r_kernel=dict(raw.get("r_kernel", {}) or {}),
            raw=dict(raw),
        )

    # ---- figure helpers ----------------------------------------------------
    def figure_mirnas(self, organism: str) -> dict[str, str]:
        """``{mirna: example_name}`` of the prespecified figure miRNAs for one organism."""
        out: dict[str, str] = {}
        for name, spec in self.figure_examples.items():
            mirnas = spec.get("mirnas", {}) if isinstance(spec, dict) else {}
            if isinstance(mirnas, dict) and organism in mirnas and mirnas[organism]:
                out[str(mirnas[organism])] = name
        return out

    def figure_terms_for(self, organism: str, family: str) -> list[str]:
        """Prespecified figure terms in the organism's own identifier space."""
        ids = self.figure_terms.get(family, [])
        if family == "REACTOME":
            prefix = REACTOME_PREFIX.get(organism)
            if prefix is None:
                return []
            return [t if t.startswith("R-") else f"{prefix}{t}" for t in ids]
        return list(ids)


# --------------------------------------------------------------------------- #
# Scope
# --------------------------------------------------------------------------- #
@dataclass
class SizeConfig:
    organism: str
    family: str
    tier: str
    label: str  # p10 / p50 / p90 (deduplicated labels joined by '+')
    n: int

    @property
    def stratum_key(self) -> str:
        return f"{self.organism}_{self.family}_{self.tier}"

    @property
    def key(self) -> str:
        return f"{self.stratum_key}_{self.label}_n{self.n}"


@dataclass
class StratumScope:
    organism: str
    family: str
    tier: str
    available: bool
    N: int
    m_full: int
    n_eligible: int
    n_manifest: int | None
    n_empty: int
    n_raw_mirnas: int
    cohort_size: int
    size_source: str
    sizes: dict[str, int]
    normalization: dict[str, int]
    note: str = ""
    n_universe_file: int | None = None
    universe_equals_annotated_union: bool | None = None

    @property
    def key(self) -> str:
        return f"{self.organism}_{self.family}_{self.tier}"


def _sizes_for(
    stratum: Stratum, cohort: set[str], cfg: StudyConfig
) -> tuple[dict[str, int], str, int]:
    sizes = stratum.query_sizes()
    if len(cohort) >= cfg.cohort_min_size:
        base = [sizes[m] for m in sorted(cohort) if m in sizes]
        src = "shared_cohort"
    else:
        base = list(sizes.values())
        src = "all_eligible"
    return percentile_sizes(base, cfg.query_size_percentiles), src, len(base)


def build_scope(
    inputs: ReleaseInputs,
    cfg: StudyConfig,
    *,
    provenance_dir: Path | None = None,
    strata_filter: Iterable[str] | None = None,
    include_grid: bool = True,
    log=print,
) -> tuple[list[StratumScope], list[SizeConfig], dict[str, set[str]]]:
    """Enumerate available strata, the shared cohorts, and the pilot query sizes.

    Returns ``(strata, size_configs, cohorts)`` where ``cohorts`` maps
    ``"<org>_<family>"`` to the miRNAs eligible in EVERY requested tier of that
    organism x family (the shared cohort used for tier comparisons).
    """
    manifest_counts = load_manifest_cell_counts(provenance_dir) if provenance_dir else {}
    wanted = set(strata_filter) if strata_filter else None
    strata: list[StratumScope] = []
    configs: list[SizeConfig] = []
    cohorts: dict[str, set[str]] = {}
    for org in cfg.organisms:
        if org not in inputs.organisms():
            log(f"[scope] {org}: no enrichment_inputs directory -> skipped")
            continue
        for fam in cfg.families:
            if fam not in inputs.available_families(org):
                log(f"[scope] {org}/{fam}: not shipped for this organism -> skipped")
                strata.append(
                    StratumScope(
                        org,
                        fam,
                        "*",
                        False,
                        0,
                        0,
                        0,
                        None,
                        0,
                        0,
                        0,
                        "n/a",
                        {},
                        {},
                        "family not shipped",
                    )
                )
                continue
            # The shared cohort and the percentile sizes are STUDY DEFINITIONS: they are
            # always computed over every requested tier, whatever subset of strata a run
            # is restricted to, so a filtered run can never change a configuration key.
            if wanted is not None and not any(f"{org}_{fam}_{t}" in wanted for t in cfg.tiers):
                continue
            per_tier: dict[str, Stratum] = {
                tier: inputs.stratum(org, fam, tier) for tier in cfg.tiers
            }
            eligible_sets = [set(s.queries) for s in per_tier.values()]
            cohort = set.intersection(*eligible_sets) if eligible_sets else set()
            cohorts[f"{org}_{fam}"] = cohort
            for tier, st in per_tier.items():
                if wanted is not None and f"{org}_{fam}_{tier}" not in wanted:
                    inputs.drop_cached_stratum(org, fam, tier)
                    continue
                sizes, src, n_base = _sizes_for(st, cohort, cfg)
                ss = StratumScope(
                    organism=org,
                    family=fam,
                    tier=tier,
                    available=bool(st.queries),
                    N=st.space.N,
                    m_full=st.space.n_terms,
                    n_eligible=len(st.queries),
                    n_manifest=manifest_counts.get((org, fam, tier)),
                    n_empty=len(st.empty_mirnas),
                    n_raw_mirnas=len(st.n_raw),
                    cohort_size=len(cohort),
                    size_source=src,
                    sizes=sizes,
                    normalization=st.normalization,
                    note="" if st.queries else "no eligible miRNA (empty stratum)",
                    n_universe_file=st.space.provenance.get("n_universe_file"),
                    universe_equals_annotated_union=st.space.provenance.get(
                        "universe_equals_annotated_union"
                    ),
                )
                strata.append(ss)
                log(
                    f"[scope] {ss.key}: N={ss.N} m_full={ss.m_full} eligible={ss.n_eligible} "
                    f"(manifest {ss.n_manifest}) empty={ss.n_empty} cohort={ss.cohort_size} "
                    f"sizes={sizes} ({src})"
                )
                # Deduplicate equal sizes (tiny strata collapse p10/p50/p90 onto one n).
                by_n: dict[int, list[str]] = {}
                for label, n in sizes.items():
                    by_n.setdefault(n, []).append(label)
                for n, labels in sorted(by_n.items()):
                    if n >= 1:
                        configs.append(SizeConfig(org, fam, tier, "+".join(labels), n))
                # Optional expansion beyond the percentile pilot: an explicit size grid
                # per stratum (config key extra_query_sizes: {<org>_<FAMILY>_<TIER>: [n, ...]}),
                # labelled "grid" so it is never confused with the prespecified percentiles.
                if not include_grid:
                    continue
                extra = (cfg.raw.get("extra_query_sizes") or {}).get(ss.key, [])
                for n in sorted({int(x) for x in extra}):
                    if 1 <= n <= st.space.N and n not in by_n:
                        configs.append(SizeConfig(org, fam, tier, "grid", n))
    return strata, configs, cohorts


def write_scope(
    out_dir: Path,
    strata: Sequence[StratumScope],
    configs: Sequence[SizeConfig],
    cohorts: dict[str, set[str]],
    inputs: ReleaseInputs,
    cfg: StudyConfig,
) -> None:
    scope_dir = out_dir / "scope"
    scope_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for s in strata:
        row = asdict(s)
        sizes = row.pop("sizes")
        norm = row.pop("normalization")
        for lab in ("p10", "p50", "p90"):
            row[f"size_{lab}"] = sizes.get(lab)
        row["n_kept_pairs"] = norm.get("n_kept_pairs")
        row["n_dropped_pairs"] = norm.get("n_dropped_pairs")
        row["key"] = s.key
        rows.append(row)
    pd.DataFrame(rows).to_csv(scope_dir / "strata.tsv", sep="\t", index=False)
    pd.DataFrame([asdict(c) | {"key": c.key} for c in configs]).to_csv(
        scope_dir / "query_size_configs.tsv", sep="\t", index=False
    )
    with (scope_dir / "cohorts.json").open("w", encoding="utf-8") as fh:
        json.dump({k: sorted(v) for k, v in cohorts.items()}, fh, indent=1)
    # Per-stratum actual query sizes (after mapping, family restriction, cumulative union).
    q_dir = scope_dir / "query_sizes"
    q_dir.mkdir(exist_ok=True)
    for s in strata:
        if not s.available or s.tier == "*":
            continue
        st = inputs.stratum(s.organism, s.family, s.tier)
        sizes = st.query_sizes()
        cohort = cohorts.get(f"{s.organism}_{s.family}", set())
        recs = [
            {
                "mirna": m,
                "n_raw_cumulative": st.n_raw.get(m, 0),
                "n_in_universe": sizes.get(m, 0),
                "eligible": m in sizes,
                "in_shared_cohort": m in cohort,
            }
            for m in sorted(st.n_raw)
        ]
        pd.DataFrame(recs).to_csv(q_dir / f"{s.key}.tsv", sep="\t", index=False)
    (scope_dir / "config_used.yaml").write_text(
        yaml.safe_dump(cfg.raw, sort_keys=False), encoding="utf-8"
    )


# --------------------------------------------------------------------------- #
# Manifest
# --------------------------------------------------------------------------- #
def _git(cmd: list[str], cwd: Path) -> str:
    try:
        return subprocess.run(
            ["git", *cmd], cwd=cwd, capture_output=True, text=True, check=False, timeout=30
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def module_fingerprint(package_dir: Path) -> dict[str, str]:
    """SHA-256 of every source file of the validation module (uncommitted code is still pinned)."""
    out: dict[str, str] = {}
    for p in sorted(package_dir.glob("*.py")):
        out[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def software_versions() -> dict[str, str]:
    import matplotlib
    import scipy

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "pyyaml": getattr(yaml, "__version__", "?"),
        "robustness_module": STUDY_MODULE_VERSION,
    }


def write_manifest(
    out_dir: Path,
    *,
    cfg: StudyConfig,
    inputs: ReleaseInputs,
    strata: Sequence[StratumScope],
    configs: Sequence[SizeConfig],
    repo_root: Path,
    limitations: Sequence[str],
    extra: dict[str, object] | None = None,
) -> Path:
    """Write ``manifest.json``: inputs + hashes, versions, seeds, actual scope, limitations."""
    used_files: list[dict[str, object]] = []
    seen: set[str] = set()
    for s in strata:
        if not s.available or s.tier == "*":
            continue
        space = inputs.family_space(s.organism, s.family)
        for kind in ("universe", "term2gene"):
            path = str(space.provenance.get(f"{kind}_path"))
            if path in seen:
                continue
            seen.add(path)
            used_files.append(
                {
                    "role": f"enrichment_inputs/{kind}",
                    "organism": s.organism,
                    "family": s.family,
                    "path": path,
                    "sha256": space.provenance.get(f"{kind}_sha256"),
                    "bytes": Path(path).stat().st_size if Path(path).exists() else None,
                    "n_universe_file": space.provenance.get("n_universe_file"),
                    "n_universe_effective": space.N,
                    "universe_equals_annotated_union": space.provenance.get(
                        "universe_equals_annotated_union"
                    ),
                    "n_terms_in_file": space.provenance.get("n_terms_in_file"),
                    "n_terms_eligible": space.n_terms,
                }
            )
        tp = inputs.target_pairs_provenance(s.organism)
        if tp["path"] not in seen:
            seen.add(str(tp["path"]))
            used_files.append({"role": "target_pairs", "organism": s.organism, **tp})
    manifest: dict[str, object] = {
        "study_id": cfg.study_id,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "release_info": inputs.release_info(),
        "inputs_dir": str(inputs.root),
        "input_files": used_files,
        "software": software_versions(),
        "code": {
            "repo_root": str(repo_root),
            "git_head": _git(["rev-parse", "HEAD"], repo_root),
            "git_describe": _git(["describe", "--tags", "--always", "--dirty"], repo_root),
            "git_status_short": _git(["status", "--short"], repo_root),
            "validation_module_sha256": module_fingerprint(Path(__file__).parent),
        },
        "seed": cfg.seed,
        "replicates": cfg.replicates,
        "alphas": cfg.alphas,
        "retention_alpha": cfg.retention_alpha,
        "config": cfg.raw,
        "scope": {
            "strata": [asdict(s) | {"key": s.key} for s in strata],
            "size_configs": [asdict(c) | {"key": c.key} for c in configs],
        },
        "limitations": list(limitations),
    }
    if extra:
        manifest.update(extra)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "manifest.json"
    # Keep what earlier stages recorded on top of the base manifest (production-kernel
    # provenance, stage wall times, report index) when a later run rewrites the base.
    if path.exists():
        try:
            previous = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            previous = {}
        for key, value in previous.items():
            if key not in manifest:
                manifest[key] = value
    path.write_text(json.dumps(manifest, indent=1, default=str), encoding="utf-8")
    return path


def update_manifest(out_dir: Path, **fields: object) -> None:
    """Merge fields into an existing manifest (e.g. R provenance after a kernel run)."""
    path = out_dir / "manifest.json"
    data = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    for k, v in fields.items():
        if isinstance(v, list) and isinstance(data.get(k), list):
            data[k] = list(data[k]) + [x for x in v if x not in data[k]]
        elif isinstance(v, dict) and isinstance(data.get(k), dict):
            data[k] = {**data[k], **v}
        else:
            data[k] = v
    path.write_text(json.dumps(data, indent=1, default=str), encoding="utf-8")
