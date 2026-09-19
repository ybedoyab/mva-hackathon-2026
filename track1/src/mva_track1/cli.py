"""Command-line interface for Track 1.

Scientific pipeline commands are not implemented yet.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from mva_track1 import __version__
from mva_track1.paths import default_config_path
from mva_track1.preflight import render_preflight, run_preflight

console = Console()
app = typer.Typer(
    name="mva-track1",
    help="Track 1 CLI for the 2026 Rare Disease Real Kid MVA Hackathon.",
    no_args_is_help=True,
)


@app.callback()
def main() -> None:
    """Track 1 CLI for the 2026 Rare Disease Real Kid MVA Hackathon."""


@app.command()
def info() -> None:
    """Print package information. No analysis is performed."""
    console.print("mva-track1")
    console.print(f"version: {__version__}")
    console.print("active track: Track 1 (variant prioritization)")
    console.print("status: foundation plus challenge-spec tooling")
    console.print("scientific pipeline: not implemented")


@app.command()
def preflight(
    data_dir: Path | None = typer.Option(
        None,
        "--data-dir",
        exists=False,
        file_okay=False,
        dir_okay=True,
        help="Optional local dataset directory for safe metadata inspection only.",
    ),
    config: Path | None = typer.Option(
        None,
        "--config",
        help="Optional YAML config (defaults to configs/track1.example.yaml).",
    ),
) -> None:
    """Run non-scientific environment checks. Does not prioritize variants."""
    config_path = config or default_config_path()
    report = run_preflight(data_dir=data_dir, config_path=config_path)
    render_preflight(report, console)
    if report.failed:
        raise typer.Exit(code=1)


@app.command("profile-vcf")
def profile_vcf_cmd(
    data_dir: Path | None = typer.Option(
        None,
        "--data-dir",
        exists=False,
        file_okay=False,
        dir_okay=True,
        help="Local dataset directory containing the official VCF (outside git).",
    ),
    vcf: Path | None = typer.Option(
        None,
        "--vcf",
        exists=False,
        file_okay=True,
        dir_okay=False,
        help="Optional explicit local VCF path. Filename is never printed.",
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        help="Directory for aggregate JSON/Markdown (default: local work/vcf_profile).",
    ),
) -> None:
    """Stream a local VCF and print aggregate QC only. Never prints variant records."""
    from mva_track1.vcf_profile import main as profile_main

    argv: list[str] = []
    if data_dir is not None:
        argv.extend(["--data-dir", str(data_dir)])
    if vcf is not None:
        argv.extend(["--vcf", str(vcf)])
    if output_dir is not None:
        argv.extend(["--output-dir", str(output_dir)])
    raise typer.Exit(profile_main(argv))
