"""Repository path helpers."""

from __future__ import annotations

from pathlib import Path


def repo_root(start: Path | None = None) -> Path:
    """Return the repository root containing pyproject.toml and scripts/."""
    start = (start or Path.cwd()).resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "scripts").is_dir():
            return candidate
    return start


def default_config_path(start: Path | None = None) -> Path:
    return repo_root(start) / "configs" / "track1.example.yaml"
