#!/usr/bin/env python3
"""Run offline Ensembl VEP 116.2 on the private annotation-only VCF."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


def _ensure_src_on_path() -> None:
    root = Path(__file__).resolve().parent.parent
    src = root / "track1" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Offline VEP; logs stay private.")
    parser.add_argument("--input-vcf", type=Path, required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--forks", type=int, default=4)
    parser.add_argument("--synonyms", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    _ensure_src_on_path()
    from mva_track1.vep_run import run_offline_vep, summarize_vep_warnings

    args = parse_args(argv)
    log_dir = args.log_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    summary = run_offline_vep(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        cache_dir=args.cache_dir,
        log_path=log_dir / "vep116.stderr.txt",
        warning_file=log_dir / "vep116.warnings.txt",
        stats_file=log_dir / "vep116.stats.html",
        forks=args.forks,
        synonyms=args.synonyms,
    )
    summary["elapsed_seconds"] = round(time.perf_counter() - started, 1)
    summary.update(
        summarize_vep_warnings(log_dir / "vep116.warnings.txt", log_dir / "vep116.stderr.txt")
    )
    text = (
        "VEP run summary (aggregates only)\n"
        f"returncode={summary['returncode']} forks={summary['forks']} "
        f"elapsed_s={summary['elapsed_seconds']}\n"
        f"output_exists={summary['output_exists']} "
        f"output_bytes={summary['output_size_bytes']}\n"
        f"warnings={summary.get('warning_categories')}\n"
        f"unrecognized_contig_count={summary.get('unrecognized_contig_count')}\n"
    )
    print(text, end="")
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0 if summary["returncode"] == 0 and summary["output_exists"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
