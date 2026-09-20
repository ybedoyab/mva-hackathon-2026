"""Private local output paths. Do not write candidate tables unless explicitly asked."""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_WORK_ROOT = Path("D:/mva-hackathon-2026-work")
DEFAULT_DATA_ROOT = Path("D:/mva-hackathon-2026-data")
DEFAULT_REFERENCE_ROOT = Path("D:/mva-reference")
DEFAULT_VEP_DATA_ROOT = Path("D:/mva-vep-data")
DEFAULT_HPO_ROOT = Path("D:/mva-hpo-data")


def work_root() -> Path:
    raw = os.environ.get("MVA_OUTPUT_DIR", "").strip()
    return Path(raw) if raw else DEFAULT_WORK_ROOT


def data_root() -> Path:
    raw = os.environ.get("MVA_DATA_DIR", "").strip()
    return Path(raw) if raw else DEFAULT_DATA_ROOT


def private_root() -> Path:
    """Directory for future patient-level tables. Must remain outside Git."""
    return work_root() / "private"


def vcf_profile_root() -> Path:
    return work_root() / "vcf_profile"


def annotation_root() -> Path:
    return work_root() / "annotation"


def normalized_root() -> Path:
    return work_root() / "normalized"


def logs_root() -> Path:
    return work_root() / "logs"


def reference_root() -> Path:
    raw = os.environ.get("MVA_REFERENCE_DIR", "").strip()
    return Path(raw) if raw else DEFAULT_REFERENCE_ROOT


def vep_data_root() -> Path:
    raw = os.environ.get("MVA_VEP_DATA_DIR", "").strip()
    return Path(raw) if raw else DEFAULT_VEP_DATA_ROOT


def hpo_data_root() -> Path:
    raw = os.environ.get("MVA_HPO_DIR", "").strip()
    return Path(raw) if raw else DEFAULT_HPO_ROOT


def hpo_release_dir(release: str = "v2026-06-23") -> Path:
    return hpo_data_root() / release


def phenotype_root() -> Path:
    return work_root() / "phenotype"


def private_phenotype_root() -> Path:
    return private_root() / "phenotype"


def ranking_root() -> Path:
    return work_root() / "ranking"


def future_private_table_path(stem: str) -> Path:
    """Return a path for a future private table. Does not create or write the file."""
    name = Path(stem).name
    if not name.endswith(".csv"):
        name = f"{name}.csv"
    return private_root() / name


def ensure_private_layout() -> Path:
    """Create the private output tree without writing candidate tables."""
    root = private_root()
    root.mkdir(parents=True, exist_ok=True)
    (root / "README.txt").write_text(
        "Private Track 1 analytical outputs. Do not copy into Git.\n",
        encoding="utf-8",
    )
    private_phenotype_root().mkdir(parents=True, exist_ok=True)
    phenotype_root().mkdir(parents=True, exist_ok=True)
    ranking_root().mkdir(parents=True, exist_ok=True)
    return root
