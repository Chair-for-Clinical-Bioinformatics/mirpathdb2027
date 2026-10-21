"""Statistical robustness checks for the production pathway ORA (validation only).

This sub-package implements a *separate* validation study over the frozen release
inputs. It never touches production outputs, the deployed services, or the
release statistics. Two questions are addressed:

1. **Calibration** of the production procedure (hypergeometric ORA with
   Benjamini-Hochberg correction over the pathways that share at least one gene
   with the query) under its stated uniform-gene null, compared with a
   fixed-family comparator (BH over every eligible term, P = 1 for zero-overlap
   terms).
2. **Sensitivity** of the biological reading to (a) the background choice
   (annotation universe vs. the tier-wide union of miRNA targets) and (b) unequal
   gene target frequency (a size-matched, degree-binned resampling null).

Production normalisation (:mod:`mirpathdb_prod.sources.production_target_sets`)
and the unchanged production kernel
(:func:`mirpathdb_prod.enrichment.clusterprofiler.run_ora_batch_with_provenance`)
are reused; the fast numpy/scipy kernel in :mod:`.kernel` is an *independent
implementation* whose agreement with the pinned R kernel is measured (not
assumed) by :mod:`.r_kernel`.

Entry point: ``scripts/validation/run_robustness_study.py`` with the checked-in
configuration ``config/robustness_study.yaml``.
"""

from __future__ import annotations

STUDY_MODULE_VERSION = "0.1.0"

__all__ = ["STUDY_MODULE_VERSION"]
