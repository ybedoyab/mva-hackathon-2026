"""Run official Ensembl VEP in Docker, fully offline, with private logs.

Never prints variant records, coordinates, alleles, genes, or sample IDs.
"""

from __future__ import annotations

import os
import re
from collections import Counter
from pathlib import Path
from typing import Any

from mva_track1.vcf_norm import run_docker_tool

VEP_IMAGE = os.environ.get("MVA_VEP_IMAGE", "ensemblorg/ensembl-vep:release_116.2")
CHROM_NOT_FOUND_RE = re.compile(
    r"chromosome\s+'?(?P<chrom>[^\s']+)'?\s+not found|not found in cache",
    re.I,
)
WARNING_RE = re.compile(r"^\s*(WARNING|WARN|ERROR|MSG)\b[:\s]", re.I)


def vep_args(
    *,
    input_vcf: Path,
    output_vcf: Path,
    cache_dir: Path,
    warning_file: Path,
    stats_file: Path,
    forks: int,
    synonyms: Path | None,
) -> list[str]:
    args = [
        "./vep",
        "--offline",
        "--cache",
        "--cache_version",
        "116",
        "--assembly",
        "GRCh38",
        "--species",
        "homo_sapiens",
        "--dir_cache",
        str(cache_dir),
        "--input_file",
        str(input_vcf),
        "--output_file",
        str(output_vcf),
        "--vcf",
        "--compress_output",
        "bgzip",
        "--force_overwrite",
        "--symbol",
        "--biotype",
        "--canonical",
        "--mane",
        "--protein",
        "--ccds",
        "--numbers",
        "--domains",
        "--variant_class",
        "--sift",
        "b",
        "--polyphen",
        "b",
        "--check_existing",
        "--clin_sig_allele",
        "1",
        "--af",
        "--af_1kg",
        "--af_gnomade",
        "--af_gnomadg",
        "--max_af",
        "--gene_phenotype",
        "--allele_number",
        "--flag_pick_allele_gene",
        "--fork",
        str(forks),
        "--warning_file",
        str(warning_file),
        "--stats_file",
        str(stats_file),
    ]
    if synonyms is not None:
        args.extend(["--synonyms", str(synonyms)])
    return args


def run_offline_vep(
    *,
    input_vcf: Path,
    output_vcf: Path,
    cache_dir: Path,
    log_path: Path,
    warning_file: Path,
    stats_file: Path,
    forks: int = 4,
    synonyms: Path | None = None,
) -> dict[str, Any]:
    mounts = [
        input_vcf.parent,
        output_vcf.parent,
        cache_dir,
        log_path.parent,
        warning_file.parent,
        stats_file.parent,
    ]
    args = vep_args(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        cache_dir=cache_dir,
        warning_file=warning_file,
        stats_file=stats_file,
        forks=forks,
        synonyms=synonyms,
    )
    proc = run_docker_tool(
        VEP_IMAGE,
        args,
        mounts=mounts,
        log_path=log_path,
        check=False,
        extra_docker_args=[
            "--network",
            "none",
            "-u",
            "0",
            "-w",
            "/opt/vep/src/ensembl-vep",
        ],
    )
    return {
        "returncode": proc.returncode,
        "image": VEP_IMAGE,
        "forks": forks,
        "network": "none",
        "output_exists": output_vcf.is_file(),
        "output_size_bytes": output_vcf.stat().st_size if output_vcf.is_file() else 0,
    }


def summarize_vep_warnings(warning_file: Path, stderr_file: Path | None = None) -> dict[str, Any]:
    """Aggregate warning categories. Does not emit contig names or coordinates."""
    counts: Counter[str] = Counter()
    unrecognized_contigs = 0
    seen_contigs: set[str] = set()
    files = [warning_file]
    if stderr_file is not None:
        files.append(stderr_file)
    for path in files:
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            match = CHROM_NOT_FOUND_RE.search(line)
            if match:
                counts["chrom_not_found"] += 1
                token = match.group("chrom") or ""
                if token and token not in seen_contigs:
                    seen_contigs.add(token)
                    unrecognized_contigs += 1
                continue
            if WARNING_RE.search(line):
                low = line.lower()
                if "skip" in low:
                    counts["skipped"] += 1
                elif "deprecated" in low:
                    counts["deprecated"] += 1
                else:
                    counts["other_warning"] += 1
    return {
        "warning_categories": dict(counts),
        "unrecognized_contig_count": unrecognized_contigs,
        "does_not_list_contig_names": True,
    }
