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
