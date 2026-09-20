"""Private local output paths. Do not write candidate tables unless explicitly asked."""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_WORK_ROOT = Path("D:/mva-hackathon-2026-work")
DEFAULT_DATA_ROOT = Path("D:/mva-hackathon-2026-data")
DEFAULT_REFERENCE_ROOT = Path("D:/mva-reference")
DEFAULT_VEP_DATA_ROOT = Path("D:/mva-vep-data")


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
        "Private Track 1 analytical outputs. Do not copy into Git.\n"
        "Candidate variant tables are not generated until a later iteration.\n",
        encoding="utf-8",
    )
    return root
