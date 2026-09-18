# Rare Disease, Real Kid: MVA Hackathon 2026

Team workspace for the [2026 Rare Disease Real Kid MVA Hackathon](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026), an open research challenge organized around real genomic and clinical data from a child living with mosaic variegated aneuploidy (MVA).

This repository is for **methods, configuration templates, documentation, and non-identifying derived artifacts**. It is not a data dump.

## Hackathon tracks

The challenge has two tracks:

- **Track 1 — variant prioritization.** Rank candidate variants from the provided genomic (and related) data. Submissions are scored against a held-out clinically confirmed answer.
- **Track 2 — drug repurposing.** Propose evidence-based hypotheses about already-approved drugs that may merit further investigation, given the disrupted biology.

We intend to participate in both tracks. **All current development is Track 1 only.** Track 2 is a placeholder and must not receive implementation effort yet.

## Scientific goal (Track 1)

Build a reproducible pipeline that inventories the local data, quality-controls it, annotates and filters variants, models inheritance, and ranks a short list of candidates for submission.

This repository does **not** contain a diagnosis, a causal gene, or a causal variant. Nothing here should be read as a clinical claim.

## Reproducibility-first

Every analysis step should be recoverable from:

- pinned software and package versions;
- recorded reference/annotation versions;
- committed example configuration (with local paths kept private);
- scripts and command lines;
- documented random seeds.

See [docs/reproducibility.md](docs/reproducibility.md) and the methodology skeleton in [docs/methodology.md](docs/methodology.md).

## Controlled data are not in this repository

Hackathon genomic and patient-level files stay on local, private storage. They must never be committed, pushed, or stored with Git LFS.

- Policy: [docs/data_policy.md](docs/data_policy.md)
- Ignore rules: [.gitignore](.gitignore)
- Pre-commit scan: `python scripts/check_repo_safety.py`

`.gitignore` is necessary but not sufficient. Run the safety checker before every commit.

## Development setup

Python 3.11 or newer is required.

```bash
python -m venv .venv
```

Activate the environment:

- Windows (PowerShell): `.venv\Scripts\Activate.ps1`
- macOS / Linux: `source .venv/bin/activate`

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Useful commands:

```bash
mva-track1 --help
mva-track1 info
mva-track1 preflight
pytest
ruff check .
python scripts/check_repo_safety.py
python scripts/download_track1_data.py --output-dir <dir-outside-repo> --dry-run
```

Copy [configs/track1.example.yaml](configs/track1.example.yaml) to an untracked local config before pointing at private files. Copy [.env.example](.env.example) to `.env` if you need local secrets; never commit `.env`.

## Repository structure

```text
mva-hackathon-2026/
├── README.md
├── .gitignore
├── .env.example
├── pyproject.toml
├── scripts/
│   └── check_repo_safety.py
├── configs/
│   └── track1.example.yaml
├── docs/
│   ├── data_policy.md
│   ├── methodology.md
│   ├── reproducibility.md
│   └── track1_challenge_spec.md
├── shared/
│   └── utils/
├── track1/                  # active work area
│   ├── src/mva_track1/
│   ├── scripts/
│   ├── tests/
│   ├── results/
│   └── report/
└── track2/                  # placeholder only
```

## License and data access

Code in this repository is team work product for the hackathon. Access to the controlled dataset is governed by the hackathon organizers, not by this repo. Do not redistribute that data.
