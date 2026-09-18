# Reproducibility

Track 1 work should be recoverable by a teammate from this repository plus **local** controlled data. This document describes what we will record. The pipeline itself is not implemented yet.

## Software versions

Record:

- Python version (`python --version`);
- the project version from `pyproject.toml`;
- installed package versions (`python -m pip freeze` or an equivalent lock once we introduce one);
- versions of any external CLI tools (for example aligners or annotators) when they are added.

Prefer declaring dependencies in `pyproject.toml` and installing with `pip install -e ".[dev]"`.

## Reference genome and build

Once the local data inventory identifies the assembly, record the build name (for example GRCh38) and, when available, the exact FASTA filename, source, and checksum. Do not copy the reference genome into Git if it is large or redistributable only under restricted terms; store the identifier and checksum instead.

## Annotation database versions

For each annotator and frequency catalog we eventually use, record product name, version or date stamp, and genome build compatibility. Keep this next to the config that selected those sources.

## Command-line parameters

Every analysis command should be either:

- encoded in a script in this repository, or
- written into a log next to the outputs, including the working directory and the config path.

Do not rely on unrecorded interactive notebook state for a submission.

## Random seeds

The example config sets `project.random_seed: 2026`. Any stochastic step must read the seed from configuration and write the seed used into the run record. If a step is fully deterministic, say so.

## Configuration

Committed files under `configs/` are **examples** only. Local configs that contain real paths stay untracked.

A reproducible run record should store a copy of the resolved config (with secrets redacted) alongside outputs.

## Pipeline steps

Follow the Track 1 skeleton in [methodology.md](methodology.md). When a step is implemented, document:

- inputs;
- software;
- outputs;
- whether outputs may be committed (see [data_policy.md](data_policy.md)).

## Checksums for local inputs

For private FASTQ/BAM/VCF (and similar) files, record cryptographic hashes (for example SHA-256) in a local run log or an untracked manifest so we can detect silent file substitution.

Do **not** commit the input files themselves. Hashes do not replace access control; they only identify which local bytes were used.

## Outputs

Prefer writing pipeline outputs **outside** the repository. If a derived file is later committed under `track1/results/`, it must pass the data policy review and `python scripts/check_repo_safety.py`.
