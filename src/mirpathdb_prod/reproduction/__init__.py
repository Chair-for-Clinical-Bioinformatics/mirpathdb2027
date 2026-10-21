"""Controlled v2 reproduction (engine-isolation comparison).

This subpackage glues the target-source builders
(:mod:`mirpathdb_prod.sources.target_sets`), the v2 gene-set reconstruction
(:mod:`mirpathdb_prod.pathways.v2_genesets`), and the clusterProfiler batch
kernel (:mod:`mirpathdb_prod.enrichment.clusterprofiler`) into the end-to-end
"build our (miRNA x pathway) p-value matrix the way v2 did, with the new engine"
pipeline used to validate engine parity against the frozen v2 matrices.
"""

from __future__ import annotations
