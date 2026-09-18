# Track 1

This is the **active work area** for the 2026 Rare Disease Real Kid MVA Hackathon.

Track 1 is computational variant prioritization: rank candidate variants from the local genomic (and related) data for organizer scoring. No ranking method is implemented in this repository yet, and no causal variant is asserted here.

Track 2 lives in `../track2/` as a placeholder only.

## Layout

- `src/mva_track1/` — installable Python package and CLI (`mva-track1`)
- `scripts/` — Track 1 helper scripts (none yet)
- `tests/` — unit tests
- `results/` — optional publishable derived artifacts after review
- `report/` — written report materials

Controlled inputs do not belong in any of these directories. See [../docs/data_policy.md](../docs/data_policy.md).

## Provisional pipeline

The following sequence is a **provisional** planning outline. It is not a claim that this will be the final architecture, and every stage is unimplemented:

```text
data inventory
→ QC
→ variant normalization
→ annotation
→ rarity filtering
→ inheritance analysis
→ phenotype matching
→ pathogenicity evidence
→ read-level validation
→ candidate ranking
→ submission
```

Details: [../docs/methodology.md](../docs/methodology.md).

## Commands

```bash
mva-track1 --help
mva-track1 info
mva-track1 preflight
mva-track1 preflight --data-dir PATH
```

Example configuration (placeholders only): [../configs/track1.example.yaml](../configs/track1.example.yaml).

Challenge specification: [../docs/track1_challenge_spec.md](../docs/track1_challenge_spec.md).

Safe local inventory (does not read variant records or clinical text):

```bash
python scripts/inspect_track1_dataset.py --data-dir PATH --json-output PATH
```

Dataset download is **not** run from tests. Dry-run only unless `--confirm-download` is passed, and the destination must be outside this git tree:

```bash
python scripts/download_track1_data.py --output-dir PATH --dry-run
```
