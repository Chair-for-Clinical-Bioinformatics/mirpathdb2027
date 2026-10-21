"""The paired old-versus-candidate summary is built from the replicate tables of both kernels."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mirpathdb_prod.validation.robustness.report import paired_table

REP_HEADER = (
    "rep\tn\tm_overlap\tmin_p\tmin_padj_overlap\tmin_padj_fixed\t"
    "rej_overlap_0.005\trej_fixed_0.005\trej_overlap_0.05\trej_fixed_0.05\t"
    "ret_overlap_0.25\tret_fixed_0.25\n"
)


def _write(out: Path, key: str) -> None:
    cal = out / "calibration"
    (cal / "r_kernel").mkdir(parents=True)
    (cal / f"{key}.json").write_text(
        json.dumps(
            {
                "key": key,
                "status": "complete",
                "organism": "hsa",
                "family": "GO_BP",
                "tier": "HIGH",
                "size_label": "p10",
                "n": 3,
                "m_full": 100,
            }
        ),
        encoding="utf-8",
    )
    # rep: (rej_o_005, rej_f_005, rej_o_05, rej_f_05, ret_o, ret_f)
    reps = [
        (0, 0, 4, 0, 4, 0),  # production rejects at 0.05, candidate does not
        (1, 1, 2, 1, 6, 2),  # both reject
        (0, 0, 0, 0, 1, 0),  # neither rejects; retained only under production
        (0, 0, 3, 0, 3, 0),  # production only
    ]
    body = "".join(
        f"{i}\t3\t10\t0.001\t0.01\t0.5\t{a}\t{b}\t{c}\t{d}\t{e}\t{f}\n"
        for i, (a, b, c, d, e, f) in enumerate(reps)
    )
    (cal / f"{key}_replicates.tsv").write_text(REP_HEADER + body, encoding="utf-8")
    (cal / "r_kernel" / f"{key}_replicates.tsv").write_text(REP_HEADER + body, encoding="utf-8")


def test_paired_table_reports_excess_with_interval_and_zero_reverse(tmp_path: Path) -> None:
    _write(tmp_path, "hsa_GO_BP_HIGH_p10_n3")
    df = paired_table(tmp_path)
    assert set(df.kernel) == {"independent_kernel", "production_kernel"}
    r = df[(df.kernel == "production_kernel") & (df.measure == "rejection") & (df.alpha == 0.05)]
    assert len(r) == 1
    r = r.iloc[0]
    assert r.production_any == 0.75 and r.candidate_any == 0.25
    assert r.excess_any == 0.5 and r.reverse_any == 0.0
    assert r.excess_ci95_low < 0.5 < r.excess_ci95_high
    assert r.candidate_never_exceeds_production
    assert r.mean_paired_difference == (4 + 1 + 0 + 3) / 4
    ret = df[(df.kernel == "independent_kernel") & (df.measure == "retention")].iloc[0]
    assert ret.alpha == 0.25 and ret.production_any == 1.0 and ret.candidate_any == 0.25
    assert ret.excess_any == 0.75


def test_paired_table_ignores_incomplete_or_unknown_configurations(tmp_path: Path) -> None:
    _write(tmp_path, "hsa_GO_BP_HIGH_p10_n3")
    (tmp_path / "calibration" / "orphan_replicates.tsv").write_text(REP_HEADER, encoding="utf-8")
    df = paired_table(tmp_path)
    assert set(df.key) == {"hsa_GO_BP_HIGH_p10_n3"}
    assert isinstance(df, pd.DataFrame)
