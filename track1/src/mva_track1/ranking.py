"""Transparent compound-het ranking: genotype-only and phenotype-aware lanes.

Weights are fixed and documented. No model is trained on the patient.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from mva_track1.candidates import (
    AF_MISSING,
    PairCandidate,
    VariantCandidate,
    build_pairs,
    frequency_knowledge,
    pair_class_counts,
    public_pair_row,
)
from mva_track1.hpo_ontology import HpoOntology
from mva_track1.phenotype_score import PhenotypeScores, score_gene_phenotype

# Genotype component weights (sum to 1.0).
W_FREQUENCY = 0.35
W_FUNCTIONAL = 0.30
W_QUALITY = 0.15
W_CLINVAR = 0.15
W_PREDICTION = 0.05

# Lane B mix. Negative phenotype is reported but not used to drop candidates.
W_GENOTYPE_LANE_B = 0.55
W_PHENOTYPE_LANE_B = 0.45
W_NEGATIVE_PENALTY = 0.05

PHASE_CIS_MULTIPLIER = 0.15

SENSITIVITY_KEYS = ("S1", "S2", "S3", "S4", "S5", "S6", "S7")


def frequency_score(item: VariantCandidate) -> float:
    mapping = {
        "ZERO_REPORTED": 1.00,
        "OBSERVED_ULTRARARE": 0.90,
        "OBSERVED_RARE": 0.55,
        "AF_MISSING": 0.30,
        "COMMON": 0.00,
    }
    return mapping.get(item.frequency_class, 0.0)


def functional_score(item: VariantCandidate) -> float:
    cons = item.consequence
    if cons in {"splice_acceptor_variant", "splice_donor_variant", "transcript_ablation"}:
        return 1.00
    if cons in {"frameshift_variant", "stop_gained"}:
        return 0.95
    if cons in {"start_lost", "stop_lost"}:
        return 0.85
    if item.impact == "HIGH":
        return 0.80
    if cons == "missense_variant":
        return 0.55
    if cons in {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}:
        return 0.50
    if item.impact == "MODERATE":
        return 0.40
    return 0.10


def clinvar_score(item: VariantCandidate) -> float:
    mapping = {
        "pathogenic": 1.00,
        "likely_pathogenic": 0.90,
        "pathogenic/likely_pathogenic": 0.95,
        "conflicting": 0.15,
        "VUS": 0.05,
        "missing": 0.00,
        "other": 0.00,
        "likely_benign": -0.25,
        "benign": -0.40,
    }
    return mapping.get(item.clinvar, 0.0)


def prediction_score(item: VariantCandidate) -> float:
    if item.consequence != "missense_variant":
        return 0.0
    sift = (item.sift or "").casefold()
    poly = (item.polyphen or "").casefold()
    sift_del = sift.startswith("deleterious")
    poly_dam = poly.startswith("probably_damaging")
    poly_pos = poly.startswith("possibly_damaging")
    sift_tol = sift.startswith("tolerated")
    poly_ben = poly.startswith("benign")
    if sift_del and poly_dam:
        return 0.20
    if sift_del or poly_dam:
        return 0.10
    if poly_pos:
        return 0.05
    if sift_tol and poly_ben:
        return -0.05
    return 0.0


def quality_score(item: VariantCandidate) -> float:
    base = 1.0 if item.filt in {"PASS", ".", ""} else 0.35
    if item.dp is not None and item.dp < 10:
        base *= 0.5
    elif item.dp is not None and item.dp < 20:
        base *= 0.8
    if item.gq is not None and item.gq < 20:
        base *= 0.5
    elif item.gq is not None and item.gq < 30:
        base *= 0.85
    if item.ab is not None and item.zygosity == "het" and (item.ab < 0.15 or item.ab > 0.85):
        base *= 0.6
    return min(base, 1.0)


def variant_genotype_score(item: VariantCandidate) -> dict[str, float]:
    freq = frequency_score(item)
    func = functional_score(item)
    qual = quality_score(item)
    clin = clinvar_score(item)
    pred = prediction_score(item)
    total = (
        W_FREQUENCY * freq
        + W_FUNCTIONAL * func
        + W_QUALITY * qual
        + W_CLINVAR * clin
        + W_PREDICTION * pred
    )
    return {
        "frequency_score": round(freq, 4),
        "functional_score": round(func, 4),
        "quality_score": round(qual, 4),
        "clinvar_score": round(clin, 4),
        "prediction_score": round(pred, 4),
        "genotype_score": round(max(0.0, min(total, 1.0)), 4),
    }


def phase_multiplier(phase_class: str) -> float:
    if phase_class == "likely_cis":
        return PHASE_CIS_MULTIPLIER
    return 1.0


def pair_component_scores(pair: PairCandidate) -> dict[str, float]:
    left = variant_genotype_score(pair.left)
    right = variant_genotype_score(pair.right)
    phase_s = phase_multiplier(pair.phase_class)
    genotype = 0.5 * (left["genotype_score"] + right["genotype_score"]) * phase_s
    return {
        "frequency_score": round(0.5 * (left["frequency_score"] + right["frequency_score"]), 4),
        "functional_score": round(0.5 * (left["functional_score"] + right["functional_score"]), 4),
        "quality_score": round(0.5 * (left["quality_score"] + right["quality_score"]), 4),
        "clinvar_score": round(0.5 * (left["clinvar_score"] + right["clinvar_score"]), 4),
        "prediction_score": round(0.5 * (left["prediction_score"] + right["prediction_score"]), 4),
        "phase_score": round(phase_s, 4),
        "genotype_score": round(genotype, 4),
        "left_genotype_score": left["genotype_score"],
        "right_genotype_score": right["genotype_score"],
    }


def combine_lane_b(genotype_score: float, phenotype_score: float, negative: float) -> float:
    raw = W_GENOTYPE_LANE_B * genotype_score + W_PHENOTYPE_LANE_B * phenotype_score
    raw -= W_NEGATIVE_PENALTY * negative
    return round(max(0.0, min(raw, 1.0)), 4)


@dataclass
class RankedPair:
    pair: PairCandidate
    scores: dict[str, float]
    phenotype: PhenotypeScores | None
    overall: float
    lane: str


def _passes_settings(
    pair: PairCandidate,
    *,
    exclude_missing_af: bool = False,
    pass_only: bool = False,
) -> bool:
    variants = (pair.left, pair.right)
    if exclude_missing_af and any(item.frequency_class == AF_MISSING for item in variants):
        return False
    if pass_only and any(item.filt not in {"PASS", ".", ""} for item in variants):
        return False
    return True


def score_pairs(
    pairs: Iterable[PairCandidate],
    ontology: HpoOntology | None,
    present_hpo: list[str],
    absent_hpo: list[str],
    *,
    lane: str,
    use_phenotype: bool,
) -> list[RankedPair]:
    cache: dict[str, PhenotypeScores] = {}
    ranked: list[RankedPair] = []
    for pair in pairs:
        components = pair_component_scores(pair)
        pheno = None
        pheno_pos = 0.0
        pheno_neg = 0.0
        if use_phenotype and ontology is not None:
            if pair.gene not in cache:
                cache[pair.gene] = score_gene_phenotype(
                    pair.gene, ontology, present_hpo, absent_hpo
                )
            pheno = cache[pair.gene]
            pheno_pos = pheno.positive_phenotype_score
            pheno_neg = pheno.negative_contradiction_score
        if lane == "genotype_only" or not use_phenotype:
            overall = components["genotype_score"]
        else:
            overall = combine_lane_b(components["genotype_score"], pheno_pos, pheno_neg)
        scores = dict(components)
        scores["phenotype_score"] = round(pheno_pos, 4)
        scores["negative_contradiction_score"] = round(pheno_neg, 4)
        scores["overall_score"] = overall
        ranked.append(
            RankedPair(pair=pair, scores=scores, phenotype=pheno, overall=overall, lane=lane)
        )
    ranked.sort(key=lambda item: (-item.overall, item.pair.gene, item.pair.left.variant_id))
    return ranked


def best_pair_per_gene(ranked: list[RankedPair]) -> list[RankedPair]:
    best: dict[str, RankedPair] = {}
    for item in ranked:
        current = best.get(item.pair.gene)
        if current is None or item.overall > current.overall:
            best[item.pair.gene] = item
    ordered = list(best.values())
    ordered.sort(key=lambda item: (-item.overall, item.pair.gene))
    return ordered


def gene_rank_map(ranked_genes: list[RankedPair]) -> dict[str, int]:
    return {item.pair.gene: i + 1 for i, item in enumerate(ranked_genes)}


def sensitivity_analyses(
    pairs: list[PairCandidate],
    ontology: HpoOntology | None,
    present_high: list[str],
    present_medium: list[str],
    absent_hpo: list[str],
) -> dict[str, list[RankedPair]]:
    high_med = list(dict.fromkeys(present_high + present_medium))
    configs = {
        "S1": {
            "present": present_high,
            "use_phenotype": True,
            "exclude_missing_af": False,
            "pass_only": False,
        },
        "S2": {
            "present": high_med,
            "use_phenotype": True,
            "exclude_missing_af": False,
            "pass_only": False,
        },
        "S3": {
            "present": present_high,
            "use_phenotype": False,
            "exclude_missing_af": False,
            "pass_only": False,
        },
        "S4": {
            "present": present_high,
            "use_phenotype": True,
            "exclude_missing_af": True,
            "pass_only": False,
        },
        "S5": {
            "present": present_high,
            "use_phenotype": True,
            "exclude_missing_af": False,
            "pass_only": False,
        },
        "S6": {
            "present": present_high,
            "use_phenotype": True,
            "exclude_missing_af": False,
            "pass_only": True,
        },
        "S7": {
            "present": present_high,
            "use_phenotype": True,
            "exclude_missing_af": False,
            "pass_only": False,
        },
    }
    out: dict[str, list[RankedPair]] = {}
    for key, cfg in configs.items():
        filtered = [
            pair
            for pair in pairs
            if _passes_settings(
                pair,
                exclude_missing_af=bool(cfg["exclude_missing_af"]),
                pass_only=bool(cfg["pass_only"]),
            )
        ]
        ranked = score_pairs(
            filtered,
            ontology,
            list(cfg["present"]),
            absent_hpo,
            lane="phenotype_aware" if cfg["use_phenotype"] else "genotype_only",
            use_phenotype=bool(cfg["use_phenotype"]),
        )
        out[key] = best_pair_per_gene(ranked)
    return out


def rank_stability(
    analyses: dict[str, list[RankedPair]], genes: list[str]
) -> dict[str, dict[str, int | None]]:
    maps = {key: gene_rank_map(rows) for key, rows in analyses.items()}
    stability: dict[str, dict[str, int | None]] = {}
    for gene in genes:
        stability[gene] = {key: maps[key].get(gene) for key in SENSITIVITY_KEYS}
    return stability


def clinvar_pair_category(pair: PairCandidate) -> str:
    cats = {pair.left.clinvar, pair.right.clinvar}
    positive = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
    if cats & positive:
        if "conflicting" in cats:
            return "any_pathogenic_plus_conflicting"
        return "any_pathogenic"
    if "conflicting" in cats:
        return "conflicting"
    if "VUS" in cats:
        return "vus"
    if cats <= {"missing", "other"}:
        return "missing"
    if cats & {"benign", "likely_benign"}:
        return "benign_or_likely_benign"
    return "other"


def missing_af_audit(pair: PairCandidate) -> dict[str, Any]:
    rows = []
    for item in (pair.left, pair.right):
        rows.append(
            {
                "variant_id": item.variant_id,
                "frequency_class": item.frequency_class,
                "max_af_missing": item.max_af is None,
                "gnomadg_present": item.gnomadg_af is not None,
                "existing_variant_status": item.existing_status,
                "frequency_knowledge": frequency_knowledge(item.max_af, item.gnomadg_af),
            }
        )
    knowledge = {row["frequency_knowledge"] for row in rows}
    if knowledge == {"frequency_verified"}:
        pair_knowledge = "frequency_verified"
    elif "frequency_unknown" in knowledge and "frequency_verified" in knowledge:
        pair_knowledge = "frequency_partially_known"
    elif knowledge == {"frequency_unknown"}:
        pair_knowledge = "frequency_unknown"
    else:
        pair_knowledge = "frequency_partially_known"
    return {"pair_frequency_knowledge": pair_knowledge, "variants": rows}


def format_gene_rank_report(ranked_genes: list[RankedPair], limit: int = 10) -> str:
    lines = ["rank\tgene\toverall\tphenotype\tgenotype\tfreq_class\tphase\tclinvar"]
    for i, item in enumerate(ranked_genes[:limit], start=1):
        pair = item.pair
        lines.append(
            "\t".join(
                [
                    str(i),
                    pair.gene,
                    f"{item.overall:.4f}",
                    f"{item.scores.get('phenotype_score', 0):.4f}",
                    f"{item.scores.get('genotype_score', 0):.4f}",
                    pair.pair_frequency_class,
                    pair.phase_class,
                    clinvar_pair_category(pair),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def write_pair_tsv(path: Path, ranked: list[RankedPair], *, public: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for item in ranked:
        if public:
            row = public_pair_row(item.pair, item.scores)
            # Flatten a bit for TSV.
            flat = {
                "gene": row["gene"],
                "overall_score": item.overall,
                "genotype_score": item.scores.get("genotype_score"),
                "phenotype_score": item.scores.get("phenotype_score"),
                "frequency_score": item.scores.get("frequency_score"),
                "functional_score": item.scores.get("functional_score"),
                "quality_score": item.scores.get("quality_score"),
                "clinvar_score": item.scores.get("clinvar_score"),
                "phase_score": item.scores.get("phase_score"),
                "pair_frequency_class": item.pair.pair_frequency_class,
                "phase_class": item.pair.phase_class,
                "phase_status": item.pair.phase_status,
                "clinvar_category": clinvar_pair_category(item.pair),
                "left_variant_id": item.pair.left.variant_id,
                "right_variant_id": item.pair.right.variant_id,
                "left_consequence": item.pair.left.consequence,
                "right_consequence": item.pair.right.consequence,
                "left_frequency_class": item.pair.left.frequency_class,
                "right_frequency_class": item.pair.right.frequency_class,
                "left_filter": item.pair.left.filt,
                "right_filter": item.pair.right.filt,
            }
        else:
            from mva_track1.candidates import private_pair_record

            flat = private_pair_record(item.pair, item.scores)
            flat["overall_score"] = item.overall
        rows.append(flat)
    if not rows:
        path.write_text("gene\n", encoding="utf-8")
        return
    # Private dicts contain nested dicts; stringify those.
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value) if isinstance(value, (dict, list)) else value
                    for key, value in row.items()
                }
            )


def summarize_phase(pairs: Iterable[PairCandidate]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for pair in pairs:
        counts[pair.phase_class] += 1
        counts[f"status:{pair.phase_status}"] += 1
    return dict(counts)


def summarize_clinvar(pairs: Iterable[PairCandidate]) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for pair in pairs:
        counts[clinvar_pair_category(pair)] += 1
    return dict(counts)


def build_ranked_lanes(
    candidates: list[VariantCandidate],
    ontology: HpoOntology | None,
    present_high: list[str],
    present_medium: list[str],
    absent_hpo: list[str],
) -> dict[str, Any]:
    pairs = build_pairs(candidates)
    genotype_ranked = score_pairs(
        pairs, ontology, present_high, absent_hpo, lane="genotype_only", use_phenotype=False
    )
    pheno_ranked = score_pairs(
        pairs, ontology, present_high, absent_hpo, lane="phenotype_aware", use_phenotype=True
    )
    genotype_genes = best_pair_per_gene(genotype_ranked)
    pheno_genes = best_pair_per_gene(pheno_ranked)
    analyses = sensitivity_analyses(pairs, ontology, present_high, present_medium, absent_hpo)
    top_genes = [item.pair.gene for item in pheno_genes[:10]]
    return {
        "pairs": pairs,
        "pair_class_counts": pair_class_counts(pairs),
        "genotype_pairs": genotype_ranked,
        "phenotype_pairs": pheno_ranked,
        "genotype_genes": genotype_genes,
        "phenotype_genes": pheno_genes,
        "sensitivity": analyses,
        "stability": rank_stability(
            analyses,
            list(dict.fromkeys(top_genes + [item.pair.gene for item in genotype_genes[:10]])),
        ),
        "phase_summary": summarize_phase(pairs),
        "clinvar_summary": summarize_clinvar(pairs),
    }
