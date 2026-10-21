"""Reporting layer.

Renders deterministic, human-readable outputs (Markdown reports, TSV/Markdown
tables) from validation results and canonical models. Outputs are stable given
the same input so they can be committed and diffed. Files land in the
configured ``validation/`` and ``docs/`` directories.
"""
