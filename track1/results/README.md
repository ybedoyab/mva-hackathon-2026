# Track 1 results

This directory is for **derived** Track 1 artifacts that might eventually be committed (tables, plots, shortlists in organizer format, and similar).

## Rules

- Raw patient-level genomic files never belong here (no FASTQ, BAM/CRAM, VCF/BCF, or equivalent).
- Files that still contain identifiable genotypes, sequences, or clinical detail should be written **outside** the repository instead.
- Every artifact must be reviewed against [../../docs/data_policy.md](../../docs/data_policy.md) before it is committed.
- Run `python scripts/check_repo_safety.py` from the repository root before committing.

The pipeline should support an output directory outside this repository (see `outputs.directory` in the example config). Using `track1/results/` is optional and only for artifacts that have passed review.
