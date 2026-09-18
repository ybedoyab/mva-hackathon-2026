#!/usr/bin/env python3
"""Fail if files that look unsafe would be included in the repository.

This script inspects Git-candidate files (tracked or untracked and not ignored).
It never prints the contents of detected secrets.

Exit codes:
    0 — no issues found
    1 — at least one issue found
    2 — the scan could not be completed
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

MAX_FILE_BYTES = 25 * 1024 * 1024
TEXT_SCAN_MAX_BYTES = 1 * 1024 * 1024

# Relative POSIX paths allowed to exceed MAX_FILE_BYTES.
LARGE_FILE_ALLOWLIST: frozenset[str] = frozenset()

GENOMIC_SUFFIXES: tuple[str, ...] = (
    ".g.vcf.gz",
    ".fastq.gz",
    ".fq.gz",
    ".vcf.gz",
    ".vcf.bgz",
    ".bed.gz",
    ".fastq",
    ".fq",
    ".bam",
    ".bai",
    ".sam",
    ".cram",
    ".crai",
    ".ubam",
    ".vcf",
    ".bcf",
    ".tbi",
    ".csi",
    ".gvcf",
    ".g.vcf",
    ".ped",
    ".fam",
    ".fast5",
    ".pod5",
    ".bgen",
    ".pgen",
    ".pvar",
    ".psam",
    ".sra",
)

CREDENTIAL_EXACT_NAMES: frozenset[str] = frozenset(
    {
        ".env",
        "credentials.json",
        "secrets.json",
        "token.json",
        "id_rsa",
        "id_rsa.pub",
        "id_dsa",
        "id_dsa.pub",
        "id_ecdsa",
        "id_ecdsa.pub",
        "id_ed25519",
        "id_ed25519.pub",
        ".netrc",
    }
)
CREDENTIAL_SUFFIXES: tuple[str, ...] = (".pem", ".key", ".p12", ".pfx")
ALLOWED_ENV_NAMES: frozenset[str] = frozenset({".env.example"})

HF_TOKEN_MARKER = "hf_"
GITHUB_TOKEN_PREFIXES: tuple[str, ...] = (
    "ghp_",
    "gho_",
    "ghu_",
    "ghs_",
    "ghr_",
    "github_pat_",
)
TOKEN_BODY_CHARS = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_")
MIN_TOKEN_BODY_LEN = 20

SKIP_DIRS: frozenset[str] = frozenset(
    {
        ".git",
        ".venv",
        "venv",
        "env",
        "ENV",
        "__pycache__",
        ".pytest_cache",
        ".mypy_cache",
        ".ruff_cache",
        ".idea",
        ".vscode",
        "data",
        "raw",
        "private",
        "downloads",
        "controlled_data",
        "tmp",
        "temp",
        "logs",
        "dist",
        "build",
        "node_modules",
        ".tox",
        ".nox",
        "htmlcov",
        ".ipynb_checkpoints",
        ".snakemake",
        ".nextflow",
        "cromwell-executions",
        "work",
        ".eggs",
        ".hypothesis",
    }
)


@dataclass(frozen=True)
class Issue:
    path: str
    reason: str


def default_root() -> Path:
    return Path(__file__).resolve().parent.parent


def relative_posix(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def _nul_split(payload: bytes) -> list[str]:
    if not payload:
        return []
    parts = payload.split(b"\0")
    return [item.decode("utf-8", errors="replace") for item in parts if item]


def git_candidate_paths(root: Path) -> list[Path] | None:
    try:
        tracked = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=root,
            capture_output=True,
            check=True,
        )
        untracked = subprocess.run(
            ["git", "ls-files", "-z", "--others", "--exclude-standard"],
            cwd=root,
            capture_output=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    names = set(_nul_split(tracked.stdout)) | set(_nul_split(untracked.stdout))
    paths: list[Path] = []
    for name in sorted(names):
        candidate = root / name
        if candidate.is_file():
            paths.append(candidate)
    return paths


def fallback_candidate_paths(root: Path) -> list[Path]:
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts[:-1]):
            continue
        files.append(path)
    return files


def iter_candidate_files(root: Path) -> list[Path]:
    git_paths = git_candidate_paths(root)
    if git_paths is not None:
        return git_paths
    return fallback_candidate_paths(root)


def is_credential_filename(name: str) -> bool:
    lower = name.lower()
    if lower in ALLOWED_ENV_NAMES:
        return False
    if lower in CREDENTIAL_EXACT_NAMES:
        return True
    if lower.startswith(".env."):
        return True
    if lower.endswith(CREDENTIAL_SUFFIXES):
        return True
    if "token" in lower or "secret" in lower:
        return True
    return False


def is_genomic_filename(name: str) -> bool:
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in GENOMIC_SUFFIXES)


def looks_like_text(path: Path) -> bool:
    try:
        chunk = path.read_bytes()[:8192]
    except OSError:
        return False
    if b"\x00" in chunk:
        return False
    try:
        chunk.decode("utf-8")
    except UnicodeDecodeError:
        return False
    return True


def _has_prefixed_token(text: str, prefix: str) -> bool:
    start = 0
    while True:
        index = text.find(prefix, start)
        if index == -1:
            return False
        body = text[index + len(prefix) :]
        taken = 0
        for char in body:
            if char not in TOKEN_BODY_CHARS:
                break
            taken += 1
        if taken >= MIN_TOKEN_BODY_LEN:
            return True
        start = index + 1


def secret_kind_labels(text: str) -> list[str]:
    labels: list[str] = []
    if _has_prefixed_token(text, HF_TOKEN_MARKER):
        labels.append("possible Hugging Face token in text file")
    if any(_has_prefixed_token(text, prefix) for prefix in GITHUB_TOKEN_PREFIXES):
        labels.append("possible GitHub personal access token in text file")
    return labels


def inspect_file(
    path: Path,
    root: Path,
    *,
    large_file_allowlist: frozenset[str] = LARGE_FILE_ALLOWLIST,
) -> list[Issue]:
    rel = relative_posix(path, root)
    issues: list[Issue] = []
    name = path.name

    if is_genomic_filename(name):
        issues.append(Issue(rel, "genomic data file extension"))
    if is_credential_filename(name):
        issues.append(Issue(rel, "credential-looking filename"))

    try:
        size = path.stat().st_size
    except OSError:
        issues.append(Issue(rel, "could not read file metadata"))
        return issues

    if size > MAX_FILE_BYTES and rel not in large_file_allowlist:
        issues.append(Issue(rel, f"file exceeds 25 MB ({size} bytes)"))

    if looks_like_text(path) and size <= TEXT_SCAN_MAX_BYTES:
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            issues.append(Issue(rel, "could not read text file"))
            return issues
        for label in secret_kind_labels(text):
            issues.append(Issue(rel, label))

    return issues


def scan_files(
    files: Iterable[Path],
    root: Path,
    *,
    large_file_allowlist: frozenset[str] = LARGE_FILE_ALLOWLIST,
) -> list[Issue]:
    issues: list[Issue] = []
    for path in files:
        issues.extend(inspect_file(path, root, large_file_allowlist=large_file_allowlist))
    return issues


def format_report(issues: list[Issue], scanned: int) -> str:
    if not issues:
        return f"Repository safety check: OK\nScanned {scanned} candidate file(s).\n"
    lines = [
        "Repository safety check: UNSAFE",
        f"Found {len(issues)} issue(s) in {scanned} candidate file(s):",
    ]
    for issue in issues:
        lines.append(f"  - {issue.reason}: {issue.path}")
    lines.append("Refusing to treat this tree as safe to commit.")
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check that Git-candidate files do not look like controlled data or secrets."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=default_root(),
        help="Repository root to scan (default: parent of scripts/).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    if not root.is_dir():
        print(f"Root does not exist or is not a directory: {root}", file=sys.stderr)
        return 2

    try:
        files = iter_candidate_files(root)
        issues = scan_files(files, root)
    except OSError as exc:
        print(f"Safety check failed: {exc}", file=sys.stderr)
        return 2

    report = format_report(issues, scanned=len(files))
    if issues:
        print(report, end="", file=sys.stderr)
        return 1
    print(report, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
