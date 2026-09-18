"""Safe local Track 1 dataset inspection.

Inspects filenames, sizes, and VCF *headers only*. Never reads variant records,
FASTQ sequences, or clinical document text.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from mva_track1.paths import default_config_path

INFO_ID_RE = re.compile(r"##INFO=<ID=([^,>]+)")
FORMAT_ID_RE = re.compile(r"##FORMAT=<ID=([^,>]+)")
CONTIG_ID_RE = re.compile(r"##contig=<ID=([^,>]+)")
FASTQ_LANE_RE = re.compile(r"_L(\d+)_R([12])(?:_|\.)", re.IGNORECASE)
VCF_SUFFIXES = (".vcf.gz", ".vcf.bgz", ".vcf")
INDEX_SUFFIXES = (".tbi", ".csi")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
PHENOTYPE_SUFFIXES = (".docx",)

PHASING_UNKNOWN = "unknown until controlled local analysis"


@dataclass
class VcfHeaderSummary:
    readable: bool = False
    sample_count: int = 0
    contig_style: str | None = None
    info_ids: list[str] = field(default_factory=list)
    format_ids: list[str] = field(default_factory=list)
    gzip_ok: bool | None = None
    chrom_header_present: bool = False


@dataclass
class DatasetManifest:
    dataset_size_bytes: int = 0
    file_type_counts: dict[str, int] = field(default_factory=dict)
    fastq_pair_count: int = 0
    fastq_lane_count: int = 0
    vcf_present: bool = False
    vcf_index_present: bool = False
    phenotype_document_present: bool = False
    genome_build_expected: str = "GRCh38"
    vcf_header: VcfHeaderSummary = field(default_factory=VcfHeaderSummary)
    phasing: str = PHASING_UNKNOWN
    warnings: list[str] = field(default_factory=list)

    def to_json_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        # Keep the documented JSON contract first-class; extra safe fields remain.
        return {
            "dataset_size_bytes": payload["dataset_size_bytes"],
            "file_type_counts": payload["file_type_counts"],
            "fastq_pair_count": payload["fastq_pair_count"],
            "vcf_present": payload["vcf_present"],
            "vcf_index_present": payload["vcf_index_present"],
            "phenotype_document_present": payload["phenotype_document_present"],
            "genome_build_expected": payload["genome_build_expected"],
            "vcf_header": {
                "readable": payload["vcf_header"]["readable"],
                "sample_count": payload["vcf_header"]["sample_count"],
                "contig_style": payload["vcf_header"]["contig_style"],
                "info_ids": payload["vcf_header"]["info_ids"],
                "format_ids": payload["vcf_header"]["format_ids"],
            },
            "fastq_lane_count": payload["fastq_lane_count"],
            "phasing": payload["phasing"],
            "warnings": payload["warnings"],
        }


def _iter_files(data_dir: Path) -> list[Path]:
    return sorted(path for path in data_dir.rglob("*") if path.is_file())


def _file_type_key(path: Path) -> str:
    name = path.name.lower()
    for suffix in (
        ".fastq.gz",
        ".fq.gz",
        ".vcf.gz",
        ".vcf.bgz",
        ".g.vcf.gz",
        ".docx",
        ".tbi",
        ".csi",
    ):
        if name.endswith(suffix):
            return suffix
    return path.suffix.lower() or "no_suffix"


def _endswith_any(name: str, suffixes: tuple[str, ...]) -> bool:
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in suffixes)


def expected_genome_build(config_path: Path | None = None) -> str:
    path = config_path or default_config_path()
    if not path.is_file():
        return "GRCh38"
    with path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    build = (loaded.get("data") or {}).get("genome_build")
    return str(build) if build else "GRCh38"


def _open_vcf_text(path: Path):
    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def read_vcf_header_only(path: Path) -> VcfHeaderSummary:
    """Read ## meta-lines and the #CHROM header. Stop before any variant record."""
    summary = VcfHeaderSummary()
    name = path.name.lower()
    gzipped = name.endswith(".gz") or name.endswith(".bgz")
    try:
        with _open_vcf_text(path) as handle:
            summary.gzip_ok = True if gzipped else None
            for line in handle:
                if line.startswith("##"):
                    info_match = INFO_ID_RE.match(line)
                    if info_match:
                        summary.info_ids.append(info_match.group(1))
                    format_match = FORMAT_ID_RE.match(line)
                    if format_match:
                        summary.format_ids.append(format_match.group(1))
                    contig_match = CONTIG_ID_RE.match(line)
                    if contig_match and summary.contig_style is None:
                        contig_id = contig_match.group(1)
                        summary.contig_style = "chr" if contig_id.startswith("chr") else "no_chr"
                    continue
                if line.startswith("#CHROM"):
                    summary.chrom_header_present = True
                    columns = line.rstrip("\n").split("\t")
                    try:
                        format_index = columns.index("FORMAT")
                        # Count samples without retaining identifiers.
                        summary.sample_count = max(0, len(columns) - format_index - 1)
                    except ValueError:
                        summary.sample_count = 0
                    summary.readable = True
                    break
                # Any other line is a record (or garbage). Do not consume further.
                break
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        summary.readable = False
        if gzipped:
            summary.gzip_ok = False
        raise VcfHeaderError(f"could not read VCF header: {exc}") from exc
    return summary


class VcfHeaderError(RuntimeError):
    pass


def _fastq_pairs(files: list[Path]) -> tuple[int, int]:
    lanes: dict[str, set[str]] = {}
    for path in files:
        if not _endswith_any(path.name, FASTQ_SUFFIXES):
            continue
        match = FASTQ_LANE_RE.search(path.name)
        if not match:
            continue
        lanes.setdefault(match.group(1), set()).add(match.group(2))
    pair_count = sum(1 for reads in lanes.values() if "1" in reads and "2" in reads)
    return pair_count, len(lanes)


def inspect_dataset(data_dir: Path, config_path: Path | None = None) -> DatasetManifest:
    if not data_dir.is_dir():
        raise FileNotFoundError(f"dataset directory not found: {data_dir}")

    files = _iter_files(data_dir)
    manifest = DatasetManifest(genome_build_expected=expected_genome_build(config_path))
    manifest.dataset_size_bytes = sum(path.stat().st_size for path in files)
    manifest.file_type_counts = dict(Counter(_file_type_key(path) for path in files))

    vcf_files = [path for path in files if _endswith_any(path.name, VCF_SUFFIXES)]
    index_files = [path for path in files if _endswith_any(path.name, INDEX_SUFFIXES)]
    phenotype_files = [path for path in files if _endswith_any(path.name, PHENOTYPE_SUFFIXES)]
    fastq_files = [path for path in files if _endswith_any(path.name, FASTQ_SUFFIXES)]

    manifest.vcf_present = bool(vcf_files)
    manifest.vcf_index_present = bool(index_files)
    manifest.phenotype_document_present = bool(phenotype_files)
    manifest.fastq_pair_count, manifest.fastq_lane_count = _fastq_pairs(files)
    if len(fastq_files) and manifest.fastq_pair_count == 0:
        manifest.warnings.append(
            "FASTQ files present but Illumina-style lane/pair names were not detected"
        )

    manifest.phasing = PHASING_UNKNOWN

    if vcf_files:
        try:
            manifest.vcf_header = read_vcf_header_only(vcf_files[0])
        except VcfHeaderError as exc:
            manifest.warnings.append(str(exc))
            manifest.vcf_header.readable = False

    return manifest


def format_summary(data_dir: Path, manifest: DatasetManifest) -> str:
    size_gb = manifest.dataset_size_bytes / (1024**3)
    header = manifest.vcf_header
    fastq_count = sum(
        manifest.file_type_counts.get(ext, 0) for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    )
    info_ids = ", ".join(header.info_ids) if header.info_ids else "(none in header)"
    format_ids = ", ".join(header.format_ids) if header.format_ids else "(none in header)"
    lines = [
        f"Dataset directory: detected ({data_dir.resolve()})",
        f"VCF: {'present' if manifest.vcf_present else 'absent'}",
        f"VCF index: {'present' if manifest.vcf_index_present else 'absent'}",
        f"FASTQ files: {fastq_count}",
        f"Paired sequencing lanes: {manifest.fastq_pair_count}",
        (
            "Clinical phenotype document: "
            f"{'present' if manifest.phenotype_document_present else 'absent'}"
        ),
        f"Reference build declared in project config: {manifest.genome_build_expected}",
        f"Total dataset size: {size_gb:.3f} GB",
        f"VCF header readable: {'yes' if header.readable else 'no'}",
        f"VCF sample count: {header.sample_count}",
        f"VCF contig style: {header.contig_style or 'unknown'}",
        f"VCF INFO field IDs: {info_ids}",
        f"VCF FORMAT field IDs: {format_ids}",
        f"Phasing: {manifest.phasing}",
    ]
    if header.gzip_ok is False:
        lines.append("VCF gzip integrity: failed")
    elif header.gzip_ok is True:
        lines.append("VCF gzip integrity: ok")
    if manifest.warnings:
        lines.append("Warnings:")
        lines.extend(f"  - {item}" for item in manifest.warnings)
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inspect a local Track 1 dataset without reading variant or clinical contents."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Local dataset directory (outside git).",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="Optional JSON path (no sequences, genotypes, or clinical text).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional Track 1 YAML config used only for the expected genome build.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        manifest = inspect_dataset(args.data_dir, config_path=args.config)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(format_summary(args.data_dir, manifest), end="")
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(manifest.to_json_dict(), indent=2) + "\n"
        args.json_output.write_text(payload, encoding="utf-8")
        print(f"Wrote safe JSON manifest: {args.json_output}")
    return 0
