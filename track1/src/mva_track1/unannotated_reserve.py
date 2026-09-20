"""Aggregate rescue-lane stats for VEP-unannotated contigs.

Does not print contig names, coordinates, or alleles. Sites remain unranked.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

from mva_track1.annotation_profile import _open_text, parse_gt_zygosity, parse_info


def classify_contig(name: str) -> str:
    raw = name.casefold()
    if "decoy" in raw:
        return "decoy"
    if "random" in raw:
        return "random"
    if "_alt" in raw or "-alt" in raw:
        return "alternate"
    if "hs38d1" in raw or "ebv" in raw or raw.endswith("_eiv"):
        return "hs38d1_or_viral"
    if raw.startswith("un") or raw.startswith("chrun") or "unplaced" in raw:
        return "unplaced"
    if "hla" in raw:
        return "hla"
    return "other"


def _vcf_contigs(path: Path) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            counts[chrom] += 1
    return dict(counts)


def unannotated_contigs(input_vcf: Path, vep_vcf: Path) -> set[str]:
    input_counts = _vcf_contigs(input_vcf)
    vep_counts = _vcf_contigs(vep_vcf)
    skipped = {chrom for chrom, n in input_counts.items() if vep_counts.get(chrom, 0) == 0}
    return skipped


def profile_unannotated_reserve(
    annotation_input_vcf: Path,
    vep_vcf: Path,
    genotype_vcf: Path | None = None,
) -> dict[str, Any]:
    skipped = unannotated_contigs(annotation_input_vcf, vep_vcf)
    category_contigs: Counter[str] = Counter()
    category_variants: Counter[str] = Counter()
    n_variants = 0
    pass_n = 0
    het_n = 0
    hom_alt_n = 0
    dp_values: list[int] = []
    gq_values: list[int] = []
    af_present = 0
    input_on_skipped = 0

    with _open_text(annotation_input_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            if chrom not in skipped:
                continue
            input_on_skipped += 1
            category = classify_contig(chrom)
            category_variants[category] += 1
            n_variants += 1
            filt = parts[6]
            if filt in {"PASS", ".", ""}:
                pass_n += 1
            info = parse_info(parts[7])
            if any(key in info for key in ("MAX_AF", "AF", "CAF", "gnomADg_AF", "gnomADe_AF")):
                af_present += 1

    # Unique contig category counts.
    category_contigs = Counter(classify_contig(name) for name in skipped)

    if genotype_vcf is not None:
        with _open_text(genotype_vcf) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10:
                    continue
                chrom = parts[0]
                if chrom not in skipped:
                    continue
                fmt_map = dict(zip(parts[8].split(":"), parts[9].split(":"), strict=False))
                zyg = parse_gt_zygosity(fmt_map.get("GT", ""))
                if zyg == "het":
                    het_n += 1
                elif zyg == "hom_alt":
                    hom_alt_n += 1
                dp_raw = fmt_map.get("DP", "")
                gq_raw = fmt_map.get("GQ", "")
                try:
                    dp_values.append(int(float(dp_raw)))
                except (TypeError, ValueError):
                    pass
                try:
                    gq_values.append(int(float(gq_raw)))
                except (TypeError, ValueError):
                    pass

    def _dist(values: list[int]) -> dict[str, float | int]:
        if not values:
            return {"n": 0}
        ordered = sorted(values)
        n = len(ordered)
        return {
            "n": n,
            "min": ordered[0],
            "p25": ordered[n // 4],
            "median": ordered[n // 2],
            "p75": ordered[(3 * n) // 4],
            "max": ordered[-1],
        }

    return {
        "rescue_lane": True,
        "phenotype_ranking_applied": False,
        "unrecognized_contig_count": len(skipped),
        "variant_count": n_variants,
        "pass_count": pass_n,
        "heterozygous_count": het_n,
        "homozygous_alt_count": hom_alt_n,
        "dp_distribution": _dist(dp_values),
        "gq_distribution": _dist(gq_values),
        "af_info_present_from_vcf": af_present,
        "contig_categories": dict(category_contigs),
        "variants_by_contig_category": dict(category_variants),
        "note": (
            "These sites were skipped by Ensembl VEP because their contigs are "
            "absent from the GRCh38 cache. They remain in a rescue lane and are "
            "not phenotype-ranked in this iteration."
        ),
    }
