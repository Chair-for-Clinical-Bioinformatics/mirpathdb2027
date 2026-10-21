"""Normalization layer.

Pure, deterministic, testable transformations from raw source records into
canonical models (e.g. :class:`mirpathdb_prod.models.Organism`). No I/O, no
network, no pass/fail judgments — given the same input these functions always
return the same output.
"""
