"""Run bcftools normalization via Docker. Logs stay private; stdout is aggregates only."""

from __future__ import annotations

import gzip
import os
import re
import subprocess
from pathlib import Path
from typing import Any

BCFTOOLS_IMAGE = os.environ.get("MVA_BCFTOOLS_IMAGE", "staphb/bcftools:1.21")
SAMTOOLS_IMAGE = os.environ.get("MVA_SAMTOOLS_IMAGE", "staphb/samtools:1.21")
MISMATCH_RE = re.compile(r"REF.*(mismatch|does not match)|does not match the reference", re.I)


def docker_available() -> bool:
    try:
        proc = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
            check=False,
        )
        return proc.returncode == 0
    except OSError:
        return False


def _win_to_docker(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":")
    rest = str(resolved)[len(resolved.drive) :].replace("\\", "/")
    return f"/{drive.lower()}{rest}"


def _host_mount(path: Path) -> str:
    """Docker Desktop on Windows accepts drive-colon paths with forward slashes."""
    resolved = path.resolve()
    return str(resolved).replace("\\", "/")


def run_docker_tool(
    image: str,
    args: list[str],
    *,
    mounts: list[Path],
    log_path: Path,
    check: bool = True,
    extra_docker_args: list[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["docker", "run", "--rm"]
    if extra_docker_args:
        cmd.extend(extra_docker_args)
    mapping: dict[Path, str] = {}
    for mount in mounts:
        root = mount.resolve()
        if not root.is_dir():
            root = root.parent
        if root not in mapping:
            dest = f"/mnt/d{len(mapping)}"
            mapping[root] = dest
            cmd.extend(["-v", f"{_host_mount(root)}:{dest}"])

    docker_args: list[str] = []
    for arg in args:
        replaced = arg
        candidate = Path(arg)
        if candidate.drive or (len(arg) > 2 and arg[1] == ":"):
            resolved = candidate.resolve()
            for root, dest in mapping.items():
                try:
                    rel = resolved.relative_to(root)
                    replaced = f"{dest}/{rel.as_posix()}"
                    break
                except ValueError:
                    continue
        docker_args.append(replaced)
    cmd.append(image)
    cmd.extend(docker_args)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("command_argc=" + str(len(cmd)) + "\n")
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=log,
            text=True,
            check=False,
        )
        if check and proc.returncode != 0:
            raise RuntimeError(f"docker tool failed with exit {proc.returncode}; see private log")
    return proc


def count_log_mismatches(log_path: Path) -> int:
    if not log_path.is_file():
        return 0
    text = log_path.read_text(encoding="utf-8", errors="replace")
    return sum(1 for line in text.splitlines() if MISMATCH_RE.search(line))


def count_vcf_records(path: Path) -> int:
    n = 0
    opener = gzip.open if path.name.lower().endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            n += 1
    return n


def run_bcftools(
    args: list[str],
    *,
    mounts: list[Path],
    log_path: Path,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    return run_docker_tool(BCFTOOLS_IMAGE, args, mounts=mounts, log_path=log_path, check=check)


def index_vcf(vcf_path: Path, log_path: Path) -> None:
    run_bcftools(
        ["bcftools", "index", "-t", "-f", str(vcf_path.resolve())],
        mounts=[vcf_path.parent],
        log_path=log_path,
    )


def faidx(fasta: Path, log_path: Path) -> None:
    run_docker_tool(
        SAMTOOLS_IMAGE,
        ["samtools", "faidx", str(fasta.resolve())],
        mounts=[fasta.parent],
        log_path=log_path,
    )


def write_sites_only_vcf(source: Path, dest: Path) -> None:
    """Drop sample/FORMAT columns; keep INFO including ORIG."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    src_open = gzip.open if source.name.lower().endswith((".gz", ".bgz")) else open
    dst_open = gzip.open if dest.name.lower().endswith((".gz", ".bgz")) else open
    with src_open(source, "rt", encoding="utf-8", errors="replace") as inf, dst_open(
        dest, "wt", encoding="utf-8"
    ) as out:
        for line in inf:
            if line.startswith("##"):
                if line.startswith("##FORMAT"):
                    continue
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write("\t".join(line.rstrip("\n").split("\t")[:8]) + "\n")
                continue
            out.write("\t".join(line.rstrip("\n").split("\t")[:8]) + "\n")


def summarize_norm(
    *,
    original_records: int,
    normalized_path: Path,
    log_path: Path,
) -> dict[str, Any]:
    normalized_records = count_vcf_records(normalized_path)
    return {
        "original_records": original_records,
        "normalized_records": normalized_records,
        "records_created_by_split": max(0, normalized_records - original_records),
        "ref_mismatch_log_lines": count_log_mismatches(log_path),
        "normalized_size_bytes": normalized_path.stat().st_size,
        "bcftools_image": BCFTOOLS_IMAGE,
    }
