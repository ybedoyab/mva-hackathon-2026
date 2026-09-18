"""Preflight checks that do not require scientific analysis."""

from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml
from rich.console import Console

from mva_track1 import __version__
from mva_track1.dataset_inspect import format_summary, inspect_dataset
from mva_track1.download import RECOMMENDED_FREE_BYTES, free_disk_bytes, huggingface_token_present
from mva_track1.paths import default_config_path, repo_root

OPTIONAL_EXECUTABLES = ("git",)


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    required: bool = True


@dataclass
class PreflightReport:
    checks: list[CheckResult] = field(default_factory=list)
    dataset_summary: str | None = None

    @property
    def failed(self) -> bool:
        return any(not item.ok and item.required for item in self.checks)


def _python_check() -> CheckResult:
    version = sys.version.split()[0]
    ok = sys.version_info >= (3, 11)
    return CheckResult("Python version", ok, f"{version} (require >=3.11)")


def _safety_check(root: Path) -> CheckResult:
    script = root / "scripts" / "check_repo_safety.py"
    if not script.is_file():
        return CheckResult("Repository safety checker", False, f"missing {script}")
    completed = subprocess.run(
        [sys.executable, str(script), "--root", str(root)],
        check=False,
        capture_output=True,
        text=True,
    )
    detail = "exit 0" if completed.returncode == 0 else f"exit {completed.returncode}"
    return CheckResult("Repository safety checker", completed.returncode == 0, detail)


def _genome_build_check(config_path: Path) -> CheckResult:
    if not config_path.is_file():
        return CheckResult("Configured genome build", False, f"missing config {config_path}")
    with config_path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    build = (loaded.get("data") or {}).get("genome_build")
    ok = str(build) == "GRCh38"
    return CheckResult("Configured genome build", ok, str(build) if build else "(unset)")


def _token_check() -> CheckResult:
    present = huggingface_token_present()
    detail = "present (value not printed)" if present else "not set"
    return CheckResult("HF_TOKEN", True, detail, required=False)


def _disk_check(path: Path) -> CheckResult:
    free = free_disk_bytes(path)
    free_gb = free / (1024**3)
    warn = free < RECOMMENDED_FREE_BYTES
    detail = f"{free_gb:.1f} GB free"
    if warn:
        detail += " (warning: under 150 GB recommended for a full download)"
    return CheckResult("Free disk space", True, detail, required=False)


def _optional_executables() -> list[CheckResult]:
    results: list[CheckResult] = []
    for name in OPTIONAL_EXECUTABLES:
        found = shutil.which(name)
        detail = found or "not found (optional)"
        results.append(CheckResult(f"Optional executable: {name}", True, detail, required=False))
    results.append(
        CheckResult(
            "Bioinformatics CLIs",
            True,
            "not required at this stage (no aligner/annotator mandated)",
            required=False,
        )
    )
    return results


def run_preflight(
    *,
    data_dir: Path | None = None,
    config_path: Path | None = None,
    root: Path | None = None,
) -> PreflightReport:
    base = root or repo_root()
    config = config_path or default_config_path(base)
    report = PreflightReport(
        checks=[
            CheckResult("mva-track1 version", True, __version__, required=False),
            _python_check(),
            _safety_check(base),
            _genome_build_check(config),
            _token_check(),
            _disk_check(base),
            *_optional_executables(),
        ]
    )
    if data_dir is not None:
        try:
            manifest = inspect_dataset(data_dir, config_path=config)
            report.dataset_summary = format_summary(data_dir, manifest)
            report.checks.append(
                CheckResult("Dataset directory", True, "inspected (safe metadata only)")
            )
        except FileNotFoundError as exc:
            report.checks.append(CheckResult("Dataset directory", False, str(exc)))
    return report


def render_preflight(report: PreflightReport, console: Console) -> None:
    console.print("mva-track1 preflight")
    for item in report.checks:
        mark = "OK" if item.ok else "FAIL"
        suffix = "" if item.required else " [info]"
        console.print(f"  [{mark}] {item.name}: {item.detail}{suffix}")
    if report.dataset_summary:
        console.print()
        console.print("Safe dataset inspection:")
        console.print(report.dataset_summary.rstrip())
    if report.failed:
        console.print("Preflight failed.")
    else:
        console.print("Preflight passed (scientific pipeline not run).")
