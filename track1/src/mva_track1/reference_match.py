"""Compare a VCF contig dictionary to public GRCh38 reference metadata.

Reads VCF ##contig header lines only. Never reads variant records.
Never dumps the complete contig list to stdout.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from mva_track1.vcf_profile import find_proband_vcf

CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)(?:,length=(\d+))?", re.IGNORECASE)
PRIMARY = tuple(str(i) for i in range(1, 23)) + ("X", "Y", "MT")
MITO_SYNONYMS = {"M": "MT", "MT": "MT", "CHRM": "MT", "CHRM T": "MT"}

CLASS_EXACT = "exact"
CLASS_NEAR_EXACT = "near_exact"
CLASS_INCOMPATIBLE = "incompatible"
CLASS_UNKNOWN = "unknown"

CLASS_DEFINITIONS = {
    CLASS_EXACT: (
        "100% of VCF contig names and lengths match the candidate with no name "
        "transformation, and the candidate contig count equals the VCF contig count."
    ),
    CLASS_NEAR_EXACT: (
        "100% of VCF contigs match name+length after a documented synonym "
        "(strip/add 'chr' and/or mitochondrial M/MT), and either the contig counts "
        "are equal or the VCF set is a complete subset of the candidate. Primary "
        "chromosome 1-22/X/Y/MT lengths must also match."
    ),
    CLASS_INCOMPATIBLE: (
        "Primary chromosome lengths disagree, or fewer than 90% of VCF contigs "
        "match name+length even after documented synonyms."
    ),
    CLASS_UNKNOWN: "Candidate metadata could not be loaded.",
}


@dataclass(frozen=True)
class ContigDict:
    name: str
    lengths: dict[str, int]
    source: str

    @property
    def n(self) -> int:
        return len(self.lengths)


def strip_chr(name: str) -> str:
    value = name.strip()
    if value.lower().startswith("chr"):
        return value[3:]
    return value


def canonical_mito(name: str) -> str:
    core = strip_chr(name).upper()
    if core in {"M", "MT", "MITO"}:
        return "MT"
    return strip_chr(name)


def identity(name: str) -> str:
    return name


def read_vcf_contig_header(vcf_path: Path) -> ContigDict:
    """Parse ##contig header entries. Stop before variant records."""
    lengths: dict[str, int] = {}
    opener: Callable[..., Any]
    name = vcf_path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        opener = gzip.open
    else:
        opener = open
    with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                match = CONTIG_RE.match(line.rstrip("\n"))
                if match and match.group(2):
                    lengths[match.group(1)] = int(match.group(2))
                continue
            break
    return ContigDict(name="vcf_header", lengths=lengths, source=str(vcf_path.parent))


def load_fai(path: Path, name: str) -> ContigDict:
    lengths: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        if raw.startswith("@SQ"):
            fields = dict(part.split(":", 1) for part in raw.split("\t")[1:] if ":" in part)
            lengths[fields["SN"]] = int(fields["LN"])
            continue
        if raw.startswith("@"):
            continue
        parts = raw.split("\t")
        lengths[parts[0]] = int(parts[1])
    return ContigDict(name=name, lengths=lengths, source=str(path))


def load_chrom_sizes(path: Path, name: str) -> ContigDict:
    lengths: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.startswith("#"):
            continue
        chrom, size = raw.split()[:2]
        lengths[chrom] = int(size)
    return ContigDict(name=name, lengths=lengths, source=str(path))


def load_assembly_report(path: Path, name: str, *, roles: set[str] | None = None) -> ContigDict:
    """Load NCBI assembly-report Sequence-Name / Sequence-Length pairs."""
    lengths: dict[str, int] = {}
    header: list[str] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith("# Sequence-Name"):
            header = [part.strip() for part in raw[2:].split("\t")]
            continue
        if raw.startswith("#") or not raw.strip() or header is None:
            continue
        parts = raw.split("\t")
        row = dict(zip(header, parts, strict=False))
        role = row.get("Sequence-Role", "")
        if roles is not None and role not in roles:
            continue
        seq_name = row.get("Sequence-Name")
        seq_len = row.get("Sequence-Length")
        if seq_name and seq_len and seq_len.isdigit():
            lengths[seq_name] = int(seq_len)
    return ContigDict(name=name, lengths=lengths, source=str(path))


def _primary_length_matches(vcf: ContigDict, ref: dict[str, int]) -> int:
    matched = 0
    for chrom in PRIMARY:
        v_len = None
        aliases = (chrom, f"chr{chrom}", "M" if chrom == "MT" else chrom)
        mito = "chrM" if chrom == "MT" else None
        for alias in (*aliases, mito):
            if alias and alias in vcf.lengths:
                v_len = vcf.lengths[alias]
                break
        r_len = None
        for alias in (*aliases, mito):
            if alias and alias in ref:
                r_len = ref[alias]
                break
            stripped = strip_chr(alias) if alias else None
            if stripped and stripped in ref:
                r_len = ref[stripped]
                break
        if v_len is not None and r_len is not None and v_len == r_len:
            matched += 1
    return matched


def compare_contig_dicts(
    vcf: ContigDict,
    ref: ContigDict,
    *,
    transform: Callable[[str], str] = identity,
    transform_name: str = "identity",
) -> dict[str, Any]:
    ref_t = {transform(name): length for name, length in ref.lengths.items()}
    name_matches = 0
    name_len_matches = 0
    missing = 0
    for name, length in vcf.lengths.items():
        key = transform(name)
        if key in ref_t:
            name_matches += 1
            if length == ref_t[key]:
                name_len_matches += 1
        else:
            missing += 1
    extra = max(0, len(ref_t) - name_matches)
    n_vcf = max(len(vcf.lengths), 1)
    frac = name_len_matches / n_vcf
    primary_ok = _primary_length_matches(vcf, {**ref.lengths, **ref_t})
    equal_count = len(vcf.lengths) == len(ref.lengths)
    classification = classify_match(
        frac=frac,
        equal_count=equal_count,
        missing=missing,
        primary_ok=primary_ok,
        transform_name=transform_name,
    )
    return {
        "candidate": ref.name,
        "source": ref.source,
        "transform": transform_name,
        "vcf_contig_count": len(vcf.lengths),
        "ref_contig_count": len(ref.lengths),
        "exact_name_matches": name_matches,
        "exact_name_and_length_matches": name_len_matches,
        "vcf_contigs_missing_from_ref": missing,
        "ref_contigs_not_in_vcf": extra,
        "primary_chrom_length_matches": primary_ok,
        "primary_chrom_denominator": 25,
        "percent_vcf_name_length_match": round(100 * frac, 4),
        "equal_contig_count": equal_count,
        "vcf_is_complete_subset": missing == 0,
        "classification": classification,
    }


def classify_match(
    *,
    frac: float,
    equal_count: bool,
    missing: int,
    primary_ok: int,
    transform_name: str,
) -> str:
    if frac >= 1.0 and equal_count and missing == 0 and transform_name == "identity":
        return CLASS_EXACT
    if frac >= 1.0 and missing == 0 and primary_ok == 25:
        return CLASS_NEAR_EXACT
    if primary_ok < 25 or frac < 0.90:
        return CLASS_INCOMPATIBLE
    return CLASS_INCOMPATIBLE


def default_candidates(metadata_dir: Path) -> list[tuple[ContigDict, str]]:
    """Return (dict, human_label) for each available public metadata file."""
    mapping = [
        (metadata_dir / "Homo_sapiens_assembly38.fasta.fai", "broad_gatk_assembly38"),
        (metadata_dir / "Homo_sapiens_assembly38.dict", "broad_gatk_assembly38_dict"),
        (
            metadata_dir / "no_alt_plus_hs38d1.fai",
            "ncbi_grch38_no_alt_plus_hs38d1",
        ),
        (metadata_dir / "no_alt_analysis_set.fai", "ncbi_grch38_no_alt_analysis_set"),
        (metadata_dir / "full_plus_hs38d1.fai", "ncbi_grch38_full_plus_hs38d1"),
        (metadata_dir / "hg38.p14.chrom.sizes", "ucsc_hg38_p14_chrom_sizes"),
        (
            metadata_dir / "GRCh38.p14_assembly_report.txt",
            "ensembl_gencode_primary_from_ncbi_report",
        ),
        (
            metadata_dir / "GRCh38.p14_assembly_report.txt",
            "ensembl_gencode_toplevel_from_ncbi_report",
        ),
    ]
    out: list[tuple[ContigDict, str]] = []
    seen: set[str] = set()
    for path, label in mapping:
        if not path.is_file() or label in seen:
            continue
        seen.add(label)
        if path.suffix == ".txt" and "assembly_report" in path.name:
            roles = (
                {"assembled-molecule", "unlocalized-scaffold", "unplaced-scaffold"}
                if "primary" in label
                else None
            )
            out.append((load_assembly_report(path, label, roles=roles), label))
        elif path.suffix == ".sizes" or path.name.endswith("chrom.sizes"):
            out.append((load_chrom_sizes(path, label), label))
        else:
            out.append((load_fai(path, label), label))
    return out


def compare_all(vcf: ContigDict, candidates: list[tuple[ContigDict, str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for contig_dict, _label in candidates:
        for transform, tname in (
            (identity, "identity"),
            (strip_chr, "strip_chr"),
            (canonical_mito, "strip_chr_and_mito"),
        ):
            rows.append(
                compare_contig_dicts(
                    vcf, contig_dict, transform=transform, transform_name=tname
                )
            )
    return rows


def select_reference(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Pick a dominant compatible candidate. Prefer equal-count near_exact/exact."""
    exact = [row for row in rows if row["classification"] == CLASS_EXACT]
    if exact:
        exact.sort(key=lambda row: (-row["percent_vcf_name_length_match"], row["ref_contig_count"]))
        return exact[0]
    near = [row for row in rows if row["classification"] == CLASS_NEAR_EXACT]
    if not near:
        return None
    near.sort(
        key=lambda row: (
            -int(row["equal_contig_count"]),
            -row["percent_vcf_name_length_match"],
            row["ref_contigs_not_in_vcf"],
        )
    )
    return near[0]


def format_terminal_summary(rows: list[dict[str, Any]], selected: dict[str, Any] | None) -> str:
    lines = [
        "Reference contig comparison (VCF header only; no variant records)",
        f"class_definitions={CLASS_DEFINITIONS}",
        f"n_comparisons={len(rows)}",
    ]
    for row in rows:
        lines.append(
            f"{row['candidate']} transform={row['transform']} class={row['classification']} "
            f"vcf={row['vcf_contig_count']} ref={row['ref_contig_count']} "
            f"name={row['exact_name_matches']} name_len={row['exact_name_and_length_matches']} "
            f"missing={row['vcf_contigs_missing_from_ref']} extra={row['ref_contigs_not_in_vcf']} "
            f"primary={row['primary_chrom_length_matches']}/25 "
            f"pct={row['percent_vcf_name_length_match']}"
        )
    if selected:
        lines.append(
            f"selected={selected['candidate']} transform={selected['transform']} "
            f"class={selected['classification']}"
        )
    else:
        lines.append("selected=NONE (no exact or near_exact candidate)")
    return "\n".join(lines) + "\n"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare VCF contig headers to public reference metadata."
    )
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--vcf", type=Path, default=None)
    parser.add_argument("--metadata-dir", type=Path, default=Path("D:/mva-reference/metadata"))
    parser.add_argument("--json-output", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    from mva_track1.private_outputs import data_root

    args = parse_args(argv)
    try:
        vcf_path = args.vcf or find_proband_vcf(args.data_dir or data_root())
        vcf = read_vcf_contig_header(vcf_path)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    candidates = default_candidates(args.metadata_dir)
    if not candidates:
        print("No reference metadata files found.", file=sys.stderr)
        return 2
    rows = compare_all(vcf, candidates)
    selected = select_reference(rows)
    report = {
        "vcf_contig_count": vcf.n,
        "class_definitions": CLASS_DEFINITIONS,
        "comparisons": rows,
        "selected": selected,
    }
    print(format_terminal_summary(rows, selected), end="")
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0 if selected else 3


if __name__ == "__main__":
    raise SystemExit(main())
