"""Parse VEP CSQ annotations and emit aggregate statistics only.

Never print variant records, coordinates, alleles, gene names, or sample IDs.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

CONSEQUENCE_GROUPS = {
    "transcript_ablation": ("transcript_ablation",),
    "splice_acceptor_variant": ("splice_acceptor_variant",),
    "splice_donor_variant": ("splice_donor_variant",),
    "stop_gained": ("stop_gained",),
    "frameshift_variant": ("frameshift_variant",),
    "stop_lost": ("stop_lost",),
    "start_lost": ("start_lost",),
    "transcript_amplification": ("transcript_amplification",),
    "inframe_insertion": ("inframe_insertion",),
    "inframe_deletion": ("inframe_deletion",),
    "missense_variant": ("missense_variant",),
    "protein_altering_variant": ("protein_altering_variant",),
    "splice_region_variant": (
        "splice_region_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    ),
    "synonymous_variant": ("synonymous_variant", "stop_retained_variant", "start_retained_variant"),
    "intron_variant": ("intron_variant",),
    "UTR": ("5_prime_UTR_variant", "3_prime_UTR_variant"),
    "upstream_downstream": ("upstream_gene_variant", "downstream_gene_variant"),
    "intergenic": ("intergenic_variant",),
}

SEVERITY_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "feature_elongation": 9,
    "feature_truncation": 10,
    "inframe_insertion": 11,
    "inframe_deletion": 12,
    "missense_variant": 13,
    "protein_altering_variant": 14,
    "splice_donor_5th_base_variant": 15,
    "splice_region_variant": 16,
    "splice_donor_region_variant": 17,
    "splice_polypyrimidine_tract_variant": 18,
    "incomplete_terminal_codon_variant": 19,
    "start_retained_variant": 20,
    "stop_retained_variant": 21,
    "synonymous_variant": 22,
    "coding_sequence_variant": 23,
    "mature_miRNA_variant": 24,
    "5_prime_UTR_variant": 25,
    "3_prime_UTR_variant": 26,
    "non_coding_transcript_exon_variant": 27,
    "intron_variant": 28,
    "NMD_transcript_variant": 29,
    "non_coding_transcript_variant": 30,
    "upstream_gene_variant": 31,
    "downstream_gene_variant": 32,
    "TFBS_ablation": 33,
    "TFBS_amplification": 34,
    "TF_binding_site_variant": 35,
    "regulatory_region_ablation": 36,
    "regulatory_region_amplification": 37,
    "regulatory_region_variant": 38,
    "intergenic_variant": 39,
}

CSQ_RE = re.compile(r"ID=CSQ,.+Format: ([^>\"]+)")


def parse_csq_format(header_line: str) -> list[str]:
    match = CSQ_RE.search(header_line)
    if not match:
        return []
    return [part.strip() for part in match.group(1).split("|")]


def most_severe_consequence(consequences: list[str]) -> str:
    best = "other"
    best_rank = 10_000
    for raw in consequences:
        for token in raw.split("&"):
            token = token.strip()
            rank = SEVERITY_RANK.get(token, 9_000)
            if rank < best_rank:
                best_rank = rank
                best = token
    return best


def consequence_group(term: str) -> str:
    for group, members in CONSEQUENCE_GROUPS.items():
        if term in members:
            return group
    return "other"


def parse_af(value: str) -> float | None:
    if not value or value == ".":
        return None
    token = value.split("&")[0]
    try:
        return float(token)
    except ValueError:
        return None


def af_bucket(value: float | None) -> str:
    """Keep missing distinct from observed-zero / rare."""
    if value is None:
        return "missing"
    if value == 0:
        return "eq_0"
    if value <= 0.0001:
        return "gt0_le_0.0001"
    if value <= 0.001:
        return "gt_0.0001_le_0.001"
    if value <= 0.01:
        return "gt_0.001_le_0.01"
    return "gt_0.01"


def clin_sig_category(raw: str) -> str:
    text = (raw or "").strip().lower().replace(" ", "_")
    if not text or text == ".":
        return "missing"
    tokens = {part.strip() for part in re.split(r"[&/,]", text) if part.strip()}
    if "conflicting_interpretations_of_pathogenicity" in tokens or "conflicting" in text:
        return "conflicting"
    if "pathogenic/likely_pathogenic" in text or (
        "pathogenic" in tokens and "likely_pathogenic" in tokens
    ):
        return "pathogenic/likely_pathogenic"
    if "pathogenic" in tokens:
        return "pathogenic"
    if "likely_pathogenic" in tokens:
        return "likely_pathogenic"
    if "uncertain_significance" in tokens or "vus" in tokens:
        return "VUS"
    if "likely_benign" in tokens:
        return "likely_benign"
    if "benign" in tokens:
        return "benign"
    return "other"


def parse_gt_zygosity(gt: str) -> str:
    raw = (gt or "").strip()
    if raw in {"", ".", "./.", ".|."}:
        return "missing"
    parts = [part for part in re.split(r"[/|]", raw) if part != "."]
    try:
        alleles = [int(part) for part in parts]
    except ValueError:
        return "other"
    if not alleles:
        return "missing"
    if len(alleles) == 1:
        return "hom_alt" if alleles[0] > 0 else "hom_ref"
    if all(a == 0 for a in alleles):
        return "hom_ref"
    if len(set(alleles)) > 1:
        return "het"
    if alleles[0] > 0:
        return "hom_alt"
    return "other"


def phase_pair_status(pid_a: str, pid_b: str, pgt_a: str, pgt_b: str) -> str:
    """Classify two heterozygous sites. Never claims trans without evidence."""
    pid_a = (pid_a or "").strip()
    pid_b = (pid_b or "").strip()
    if not pid_a or pid_a == "." or not pid_b or pid_b == ".":
        return "insufficient"
    if pid_a != pid_b:
        return "insufficient"
    hap_a = (pgt_a or "").replace("/", "|")
    hap_b = (pgt_b or "").replace("/", "|")
    if "|" not in hap_a or "|" not in hap_b:
        return "same_phase_set_unresolved"
    a_parts = hap_a.split("|")
    b_parts = hap_b.split("|")
    if len(a_parts) != 2 or len(b_parts) != 2:
        return "same_phase_set_unresolved"
    # Same haplotype carrying the alt in PGT (e.g. 0|1 with 0|1) => likely cis.
    a_alt_on_right = a_parts[1] not in {"0", "."}
    b_alt_on_right = b_parts[1] not in {"0", "."}
    a_alt_on_left = a_parts[0] not in {"0", "."}
    b_alt_on_left = b_parts[0] not in {"0", "."}
    if (a_alt_on_right and b_alt_on_right and not a_alt_on_left and not b_alt_on_left) or (
        a_alt_on_left and b_alt_on_left and not a_alt_on_right and not b_alt_on_right
    ):
        return "same_phase"
    if (a_alt_on_left and b_alt_on_right) or (a_alt_on_right and b_alt_on_left):
        return "opposite_phase"
    return "same_phase_set_unresolved"


def _open_text(path: Path) -> TextIO:
    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            out[key] = value
        else:
            out[part] = "1"
    return out


def parse_orig(raw: str) -> dict[str, Any] | None:
    """Parse bcftools --old-rec-tag ORIG=CHROM|POS|REF|ALTS|ALLELE_INDEX."""
    text = (raw or "").strip()
    if not text or text == ".":
        return None
    parts = text.split("|")
    if len(parts) < 4:
        return None
    try:
        pos = int(parts[1])
    except ValueError:
        return None
    return {
        "chrom": parts[0],
        "pos": pos,
        "ref": parts[2],
        "alts": parts[3],
        "allele_index": parts[4] if len(parts) > 4 else "",
    }


@dataclass
class AnnotationStats:
    input_variants: int = 0
    annotated: int = 0
    skipped_no_csq: int = 0
    pass_n: int = 0
    csq_entries: int = 0
    csq_group_entries: Counter[str] = field(default_factory=Counter)
    most_severe_group: Counter[str] = field(default_factory=Counter)
    max_af: Counter[str] = field(default_factory=Counter)
    gnomade: Counter[str] = field(default_factory=Counter)
    gnomadg: Counter[str] = field(default_factory=Counter)
    clin: Counter[str] = field(default_factory=Counter)
    has_sift: int = 0
    has_polyphen: int = 0
    has_mane: int = 0
    has_canonical: int = 0
    has_gene_pheno: int = 0
    waterfall: dict[str, int] = field(default_factory=dict)
    genes_ge2_rare_het: int = 0
    genes_ge2_rare_hm_het: int = 0
    variants_per_gene_hist: Counter[str] = field(default_factory=Counter)
    candidate_pairs: int = 0
    pairs_both_pass: int = 0
    pairs_any_clinvar_p: int = 0
    pairs_both_pid: int = 0
    pairs_same_phase: int = 0
    pairs_opposite_phase: int = 0
    pairs_insufficient_phase: int = 0
    missing_af_high_moderate: int = 0
    observed_rare_high_moderate: int = 0


def _info_af(info: dict[str, str], *keys: str) -> float | None:
    for key in keys:
        if key in info:
            return parse_af(info[key])
    return None


def stream_vep_vcf(
    path: Path,
    genotypes: dict[tuple[str, int, str, str], dict[str, Any]] | None = None,
    phase_by_locus: dict[tuple[str, int], dict[str, Any]] | None = None,
) -> AnnotationStats:
    stats = AnnotationStats()
    fields: list[str] = []
    # gene -> list of rare het variant keys (hashes of chrom/pos/ref/alt not printed)
    gene_rare_hets: dict[str, list[dict[str, Any]]] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                if "ID=CSQ" in line:
                    fields = parse_csq_format(line)
                continue
            if line.startswith("#CHROM"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, _vid, ref, alt, _qual, filt, info_s = parts[:8]
            try:
                pos = int(pos_s)
            except ValueError:
                pos = -1
            stats.input_variants += 1
            passed = filt in {"PASS", ".", ""}
            if passed:
                stats.pass_n += 1
            info = parse_info(info_s)
            csq_blob = info.get("CSQ", "")
            if not csq_blob:
                stats.skipped_no_csq += 1
                continue
            stats.annotated += 1
            entries = []
            for raw in csq_blob.split(","):
                values = raw.split("|")
                row = {
                    fields[i] if i < len(fields) else f"f{i}": values[i] if i < len(values) else ""
                    for i in range(max(len(fields), len(values)))
                }
                entries.append(row)
                stats.csq_entries += 1
                cons = row.get("Consequence", "")
                stats.csq_group_entries[consequence_group(most_severe_consequence([cons]))] += 1
            cons_terms = [row.get("Consequence", "") for row in entries]
            severe = most_severe_consequence(cons_terms)
            stats.most_severe_group[consequence_group(severe)] += 1
            max_af = None
            gnomade = None
            gnomadg = None
            clin_raw = ""
            sift = False
            poly = False
            mane = False
            canon = False
            pheno = False
            gene_ids: set[str] = set()
            impacts: set[str] = set()
            for row in entries:
                max_af = max_af if max_af is not None else parse_af(row.get("MAX_AF", ""))
                if max_af is None:
                    max_af = parse_af(row.get("AF", ""))
                gnomade = gnomade if gnomade is not None else parse_af(
                    row.get("gnomADe_AF", "") or row.get("gnomAD_e_AF", "")
                )
                gnomadg = gnomadg if gnomadg is not None else parse_af(
                    row.get("gnomADg_AF", "") or row.get("gnomAD_g_AF", "")
                )
                if not clin_raw:
                    clin_raw = row.get("CLIN_SIG", "") or row.get("ClinVar", "")
                if row.get("SIFT") and row.get("SIFT") not in {".", ""}:
                    sift = True
                if row.get("PolyPhen") and row.get("PolyPhen") not in {".", ""}:
                    poly = True
                if row.get("MANE_SELECT") and row.get("MANE_SELECT") not in {".", ""}:
                    mane = True
                if row.get("CANONICAL") == "YES":
                    canon = True
                if row.get("GENE_PHENO") and row.get("GENE_PHENO") not in {".", "0", ""}:
                    pheno = True
                symbol = row.get("Gene") or row.get("SYMBOL") or ""
                if symbol and symbol != ".":
                    gene_ids.add(symbol)
                if row.get("IMPACT"):
                    impacts.add(row.get("IMPACT", ""))
            if max_af is None:
                max_af = _info_af(info, "MAX_AF")
            stats.max_af[af_bucket(max_af)] += 1
            stats.gnomade[af_bucket(gnomade)] += 1
            stats.gnomadg[af_bucket(gnomadg)] += 1
            stats.clin[clin_sig_category(clin_raw)] += 1
            stats.has_sift += int(sift)
            stats.has_polyphen += int(poly)
            stats.has_mane += int(mane)
            stats.has_canonical += int(canon)
            stats.has_gene_pheno += int(pheno)
            high_mod = "HIGH" in impacts or "MODERATE" in impacts
            orig = parse_orig(info.get("ORIG", ""))
            gt_info = _lookup_genotype(genotypes, chrom, pos, ref, alt)
            if phase_by_locus:
                locus = (orig["chrom"], orig["pos"]) if orig else (chrom, pos)
                phase = phase_by_locus.get(locus) or phase_by_locus.get(
                    (locus[0][3:] if str(locus[0]).startswith("chr") else str(locus[0]), locus[1])
                )
                if phase:
                    merged = dict(gt_info or {})
                    merged["pid"] = phase.get("pid", merged.get("pid", ""))
                    merged["pgt"] = phase.get("pgt", merged.get("pgt", ""))
                    gt_info = merged
            _update_waterfall(stats, passed, max_af, high_mod, impacts, gt_info)
            rare_observed = max_af is not None and max_af <= 0.001
            if high_mod and max_af is None:
                stats.missing_af_high_moderate += 1
            if high_mod and rare_observed:
                stats.observed_rare_high_moderate += 1
            if genotypes is None or gt_info is None:
                continue
            if gt_info.get("zygosity") != "het":
                continue
            if not (rare_observed or max_af is None):
                continue
            for gene in gene_ids:
                gene_rare_hets.setdefault(gene, []).append(
                    {
                        "pass": passed,
                        "high_mod": high_mod,
                        "clin": clin_sig_category(clin_raw),
                        "pid": gt_info.get("pid", ""),
                        "pgt": gt_info.get("pgt", ""),
                        "af_missing": max_af is None,
                        "rare": rare_observed,
                    }
                )
    _finalize_pairs(stats, gene_rare_hets)
    return stats


def _update_waterfall(
    stats: AnnotationStats,
    passed: bool,
    max_af: float | None,
    high_mod: bool,
    impacts: set[str],
    gt_info: dict[str, Any] | None,
) -> None:
    wf = stats.waterfall
    wf["all_annotated"] = wf.get("all_annotated", 0) + 1
    if not passed:
        return
    wf["pass"] = wf.get("pass", 0) + 1
    missing = max_af is None
    if missing or max_af <= 0.01:
        wf["pass_maxaf_le_0.01_or_missing"] = wf.get("pass_maxaf_le_0.01_or_missing", 0) + 1
    if not missing and max_af <= 0.01:
        wf["pass_observed_maxaf_le_0.01"] = wf.get("pass_observed_maxaf_le_0.01", 0) + 1
    if missing or max_af <= 0.001:
        wf["pass_maxaf_le_0.001_or_missing"] = wf.get("pass_maxaf_le_0.001_or_missing", 0) + 1
    if not missing and max_af <= 0.001:
        wf["pass_observed_maxaf_le_0.001"] = wf.get("pass_observed_maxaf_le_0.001", 0) + 1
    if missing or max_af <= 0.0001:
        wf["pass_maxaf_le_0.0001_or_missing"] = wf.get("pass_maxaf_le_0.0001_or_missing", 0) + 1
    if not missing and max_af <= 0.0001:
        wf["pass_observed_maxaf_le_0.0001"] = wf.get("pass_observed_maxaf_le_0.0001", 0) + 1
    rare_or_missing = missing or (max_af is not None and max_af <= 0.001)
    if rare_or_missing and "HIGH" in impacts:
        wf["rare_or_missing_HIGH"] = wf.get("rare_or_missing_HIGH", 0) + 1
        if not missing:
            wf["observed_rare_HIGH"] = wf.get("observed_rare_HIGH", 0) + 1
    if rare_or_missing and "MODERATE" in impacts:
        wf["rare_or_missing_MODERATE"] = wf.get("rare_or_missing_MODERATE", 0) + 1
        if not missing:
            wf["observed_rare_MODERATE"] = wf.get("observed_rare_MODERATE", 0) + 1
    if not high_mod or gt_info is None:
        return
    if rare_or_missing and gt_info.get("zygosity") == "het":
        wf["rare_or_missing_HIGH_MODERATE_het"] = wf.get("rare_or_missing_HIGH_MODERATE_het", 0) + 1
        if not missing:
            wf["observed_rare_HIGH_MODERATE_het"] = wf.get("observed_rare_HIGH_MODERATE_het", 0) + 1
    if rare_or_missing and gt_info.get("zygosity") == "hom_alt":
        wf["rare_or_missing_HIGH_MODERATE_hom_alt"] = (
            wf.get("rare_or_missing_HIGH_MODERATE_hom_alt", 0) + 1
        )
        if not missing:
            wf["observed_rare_HIGH_MODERATE_hom_alt"] = wf.get(
                "observed_rare_HIGH_MODERATE_hom_alt", 0
            ) + 1


def _hist_bin(n: int) -> str:
    if n <= 1:
        return "1"
    if n == 2:
        return "2"
    if n <= 5:
        return "3-5"
    if n <= 10:
        return "6-10"
    return "11+"


def _finalize_pairs(
    stats: AnnotationStats, gene_rare_hets: dict[str, list[dict[str, Any]]]
) -> None:
    for _gene, rows in gene_rare_hets.items():
        rare_rows = [row for row in rows if row["rare"] or row["af_missing"]]
        hm = [row for row in rare_rows if row["high_mod"]]
        if len(rare_rows) >= 2:
            stats.genes_ge2_rare_het += 1
            stats.variants_per_gene_hist[_hist_bin(len(rare_rows))] += 1
        if len(hm) >= 2:
            stats.genes_ge2_rare_hm_het += 1
            for i, left in enumerate(hm):
                for right in hm[i + 1 :]:
                    stats.candidate_pairs += 1
                    if left["pass"] and right["pass"]:
                        stats.pairs_both_pass += 1
                    if left["clin"] in {
                        "pathogenic",
                        "likely_pathogenic",
                        "pathogenic/likely_pathogenic",
                    } or right["clin"] in {
                        "pathogenic",
                        "likely_pathogenic",
                        "pathogenic/likely_pathogenic",
                    }:
                        stats.pairs_any_clinvar_p += 1
                    pid_ok = (
                        bool(left["pid"])
                        and left["pid"] != "."
                        and bool(right["pid"])
                        and right["pid"] != "."
                    )
                    if pid_ok:
                        stats.pairs_both_pid += 1
                    status = phase_pair_status(left["pid"], right["pid"], left["pgt"], right["pgt"])
                    if status == "same_phase":
                        stats.pairs_same_phase += 1
                    elif status == "opposite_phase":
                        stats.pairs_opposite_phase += 1
                    else:
                        stats.pairs_insufficient_phase += 1


def _lookup_genotype(
    genotypes: dict[tuple[str, int, str, str], dict[str, Any]] | None,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> dict[str, Any] | None:
    if genotypes is None:
        return None
    for key in (
        (chrom, pos, ref, alt),
        (chrom[3:] if chrom.startswith("chr") else chrom, pos, ref, alt),
        (f"chr{chrom}" if not chrom.startswith("chr") else chrom, pos, ref, alt),
    ):
        if key in genotypes:
            return genotypes[key]
    return None


def load_original_genotypes(vcf_path: Path) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Load GT/PGT/PID from a private VCF. Keys are not printed."""
    table: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, _vid, ref, alt_s, _qual, filt, _info, fmt, sample = parts[:10]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            fmt_map = dict(zip(fmt.split(":"), sample.split(":"), strict=False))
            payload = {
                "zygosity": parse_gt_zygosity(fmt_map.get("GT", "")),
                "pid": fmt_map.get("PID", ""),
                "pgt": fmt_map.get("PGT", ""),
                "filter": filt,
            }
            for alt in alt_s.split(","):
                table[(chrom, pos, ref, alt)] = payload
    return table


def load_phase_by_locus(vcf_path: Path) -> dict[tuple[str, int], dict[str, Any]]:
    """Load original PGT/PID keyed by CHROM+POS. Never prints identifiers."""
    table: dict[tuple[str, int], dict[str, Any]] = {}
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, _vid, _ref, _alt, _qual, _filt, _info, fmt, sample = parts[:10]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            fmt_map = dict(zip(fmt.split(":"), sample.split(":"), strict=False))
            table[(chrom, pos)] = {
                "pid": fmt_map.get("PID", ""),
                "pgt": fmt_map.get("PGT", ""),
            }
    return table


def stats_to_report(stats: AnnotationStats) -> dict[str, Any]:
    return {
        "input_variants": stats.input_variants,
        "annotated": stats.annotated,
        "skipped_no_csq": stats.skipped_no_csq,
        "pass": stats.pass_n,
        "csq_entry_counts": dict(stats.csq_group_entries),
        "unique_variant_most_severe": dict(stats.most_severe_group),
        "max_af_buckets": dict(stats.max_af),
        "gnomade_af_buckets": dict(stats.gnomade),
        "gnomadg_af_buckets": dict(stats.gnomadg),
        "clin_sig": dict(stats.clin),
        "prediction_availability": {
            "sift": stats.has_sift,
            "polyphen": stats.has_polyphen,
            "mane_select": stats.has_mane,
            "canonical": stats.has_canonical,
            "gene_phenotype": stats.has_gene_pheno,
            "denominator_annotated": stats.annotated,
        },
        "waterfall": stats.waterfall,
        "missing_af_kept_distinct_from_rare": True,
        "compound_het": {
            "genes_ge2_rare_het": stats.genes_ge2_rare_het,
            "genes_ge2_rare_high_moderate_het": stats.genes_ge2_rare_hm_het,
            "variants_per_gene_histogram": dict(stats.variants_per_gene_hist),
            "candidate_pairs": stats.candidate_pairs,
            "pairs_both_pass": stats.pairs_both_pass,
            "pairs_any_clinvar_pathogenic": stats.pairs_any_clinvar_p,
            "pairs_both_have_pid": stats.pairs_both_pid,
            "pairs_phase_same": stats.pairs_same_phase,
            "pairs_phase_opposite": stats.pairs_opposite_phase,
            "pairs_phase_insufficient": stats.pairs_insufficient_phase,
            "does_not_claim_trans_without_evidence": True,
        },
    }


def format_terminal_summary(report: dict[str, Any]) -> str:
    lines = [
        "VEP annotation profile (aggregates only)",
        f"input_variants={report['input_variants']}",
        f"annotated={report['annotated']} skipped_no_csq={report['skipped_no_csq']}",
        f"pass={report['pass']}",
        f"csq_entries={report['csq_entry_counts']}",
        f"most_severe={report['unique_variant_most_severe']}",
        f"max_af={report['max_af_buckets']}",
        f"clin_sig={report['clin_sig']}",
        f"predictions={report['prediction_availability']}",
        f"waterfall={report['waterfall']}",
        f"compound_het={report['compound_het']}",
    ]
    text = "\n".join(lines) + "\n"
    return text


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate VEP annotation statistics; no variant dump."
    )
    parser.add_argument("--vep-vcf", type=Path, required=True)
    parser.add_argument("--original-vcf", type=Path, default=None)
    parser.add_argument(
        "--normalized-vcf",
        type=Path,
        default=None,
        help="Post-split VCF used for allele-specific GT. PGT/PID still come from --original-vcf.",
    )
    parser.add_argument("--json-output", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    genotypes = None
    phase = None
    if args.normalized_vcf is not None:
        genotypes = load_original_genotypes(args.normalized_vcf)
    elif args.original_vcf is not None:
        genotypes = load_original_genotypes(args.original_vcf)
    if args.original_vcf is not None:
        phase = load_phase_by_locus(args.original_vcf)
    stats = stream_vep_vcf(args.vep_vcf, genotypes, phase)
    report = stats_to_report(stats)
    print(format_terminal_summary(report), end="")
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 0
