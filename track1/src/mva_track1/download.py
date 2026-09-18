"""Safe download workflow for SageBio/mva-hackathon-2026-data.

Never hard-codes tokens. Never downloads into the git repository tree.
Does not run during tests unless explicitly invoked by a test with mocks.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

from mva_track1.paths import repo_root

DATASET_REPO_ID = "SageBio/mva-hackathon-2026-data"
DATASET_REPO_TYPE = "dataset"
RECOMMENDED_FREE_BYTES = 150 * 1024**3
EXPECTED_SIZE_NOTE = "~85 GB compressed; organizers recommend 100-150 GB free"


class DownloadError(ValueError):
    """Unsafe or incomplete download request."""


@dataclass(frozen=True)
class DownloadPlan:
    output_dir: Path
    repo_id: str = DATASET_REPO_ID
    token_present: bool = False
    free_bytes: int = 0
    inside_repo: bool = False


def huggingface_token_present() -> bool:
    return bool(os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN"))


def _token_value() -> str | None:
    token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    return token or None


def is_inside_repo(path: Path, root: Path | None = None) -> bool:
    resolved = path.resolve()
    base = (root or repo_root()).resolve()
    return resolved == base or resolved.is_relative_to(base)


def free_disk_bytes(path: Path) -> int:
    probe = path if path.exists() else path.parent
    if not probe.exists():
        probe = Path.cwd()
    return shutil.disk_usage(probe).free


def plan_download(output_dir: Path, root: Path | None = None) -> DownloadPlan:
    destination = output_dir.expanduser()
    base = root or repo_root()
    inside = is_inside_repo(destination, base)
    return DownloadPlan(
        output_dir=destination,
        token_present=huggingface_token_present(),
        free_bytes=free_disk_bytes(destination),
        inside_repo=inside,
    )


def format_plan(plan: DownloadPlan, *, dry_run: bool) -> str:
    free_gb = plan.free_bytes / (1024**3)
    mode = "DRY RUN (no download)" if dry_run else "LIVE DOWNLOAD"
    token_state = "present (value not printed)" if plan.token_present else "not set in environment"
    lines = [
        f"Track 1 dataset download: {mode}",
        f"Repository: {plan.repo_id} ({DATASET_REPO_TYPE})",
        f"Destination: {plan.output_dir.resolve()}",
        f"Expected size: {EXPECTED_SIZE_NOTE}",
        f"HF token in environment: {token_state}",
        f"Free disk space: {free_gb:.1f} GB",
        "Workflow: huggingface_hub.snapshot_download with official filenames preserved.",
        "Controlled challenge data must not be committed to git or uploaded to public remotes.",
    ]
    if plan.free_bytes < RECOMMENDED_FREE_BYTES:
        lines.append(
            "WARNING: fewer than 150 GB free; organizers recommend 100-150 GB."
        )
    if plan.inside_repo:
        lines.append("ERROR: destination is inside the git repository tree.")
    return "\n".join(lines) + "\n"


def run_download(plan: DownloadPlan, *, confirm: bool, dry_run: bool) -> None:
    if plan.inside_repo:
        raise DownloadError(
            "Refusing to download controlled data into the git repository tree. "
            "Choose --output-dir outside the repo."
        )
    print(format_plan(plan, dry_run=dry_run), end="")
    if dry_run:
        print("Dry run complete. No files were downloaded.")
        return
    if not confirm:
        raise DownloadError(
            "Refusing to download without --confirm-download (or use --dry-run)."
        )

    from huggingface_hub import snapshot_download

    snapshot_download(
        repo_id=plan.repo_id,
        repo_type=DATASET_REPO_TYPE,
        local_dir=str(plan.output_dir.resolve()),
        token=_token_value(),
    )
    print("Download finished. Keep this directory private and outside git.")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download SageBio/mva-hackathon-2026-data outside the git repository."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Destination directory outside the repo.",
    )
    parser.add_argument(
        "--confirm-download",
        action="store_true",
        help="Required to actually download. Omitted for dry-run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the destination and workflow without downloading.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        plan = plan_download(args.output_dir)
        run_download(plan, confirm=args.confirm_download, dry_run=args.dry_run)
    except DownloadError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    return 0
