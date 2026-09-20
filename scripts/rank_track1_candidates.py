#!/usr/bin/env python3
"""Rank compound-het gene pairs from local VEP + HPO.

Stdout is gene-level only: no coordinates, alleles, or clinical source text.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _ensure_src_on_path() -> None:
    root = Path(__file__).resolve().parent.parent
    src = root / "track1" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Blind gene ranking; no variant dump.")
    parser.add_argument("--vep-vcf", type=Path, required=True)
    parser.add_argument("--normalized-vcf", type=Path, required=True)
    parser.add_argument("--original-vcf", type=Path, default=None)
    parser.add_argument("--annotation-input-vcf", type=Path, default=None)
    parser.add_argument("--hpo-dir", type=Path, required=True)
    parser.add_argument("--phenotype-json", type=Path, required=True)
    parser.add_argument("--genotype-tsv", type=Path, required=True)
    parser.add_argument("--phenotype-tsv", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, default=None)
    return parser.parse_args(argv)


def _hpo_lists(rows: list[dict]) -> tuple[list[str], list[str], list[str]]:
    high_present: list[str] = []
    medium_present: list[str] = []
    absent: list[str] = []
    for row in rows:
        hpo_id = row.get("hpo_id", "")
        status = row.get("status", "")
        conf = row.get("confidence", "")
        method = row.get("mapping_method", "")
        if method == "lexical_candidate":
            continue
        if row.get("scope") == "family_history":
            continue
        if status == "absent" and conf == "high":
            absent.append(hpo_id)
        elif status == "present" and conf == "high":
            high_present.append(hpo_id)
        elif status == "present" and conf == "medium":
            medium_present.append(hpo_id)
    return high_present, medium_present, absent


def main(argv: list[str] | None = None) -> int:
    _ensure_src_on_path()
    from mva_track1.candidates import (
        AF_MISSING,
        frequency_class_counts,
        heterozygous_noncommon,
        load_phase_by_locus,
        load_sample_fields,
        stream_high_moderate_candidates,
    )
    from mva_track1.hpo_ontology import (
        load_gene_disease_map,
        load_gene_phenotype_map,
        load_phenotype_hpoa,
        parse_obo,
    )
    from mva_track1.ranking import (
        build_ranked_lanes,
        clinvar_pair_category,
        format_gene_rank_report,
        missing_af_audit,
        write_pair_tsv,
    )
    from mva_track1.unannotated_reserve import profile_unannotated_reserve

    args = parse_args(argv)
    pheno_rows = json.loads(args.phenotype_json.read_text(encoding="utf-8"))
    present_high, present_medium, absent = _hpo_lists(pheno_rows)

    ontology = parse_obo(args.hpo_dir / "hp.obo")
    load_phenotype_hpoa(ontology, args.hpo_dir / "phenotype.hpoa")
    load_gene_disease_map(ontology, args.hpo_dir / "genes_to_disease.txt")
    load_gene_phenotype_map(ontology, args.hpo_dir / "genes_to_phenotype.txt")

    genotypes = load_sample_fields(args.normalized_vcf, zygosities={"het", "hom_alt"})
    phase = load_phase_by_locus(args.original_vcf) if args.original_vcf else None
    candidates = stream_high_moderate_candidates(args.vep_vcf, genotypes, phase)
    het = heterozygous_noncommon(candidates)
    result = build_ranked_lanes(candidates, ontology, present_high, present_medium, absent)

    write_pair_tsv(args.genotype_tsv, result["genotype_pairs"], public=True)
    write_pair_tsv(args.phenotype_tsv, result["phenotype_pairs"], public=True)

    reserve = None
    if args.annotation_input_vcf is not None:
        reserve = profile_unannotated_reserve(
            args.annotation_input_vcf, args.vep_vcf, args.normalized_vcf
        )

    top_pheno = result["phenotype_genes"][:10]
    missing_audits = []
    for item in top_pheno:
        if (
            item.pair.left.frequency_class == AF_MISSING
            or item.pair.right.frequency_class == AF_MISSING
        ):
            missing_audits.append({"gene": item.pair.gene, **missing_af_audit(item.pair)})

    summary = {
        "hpo_release": args.hpo_dir.name,
        "high_confidence_present_hpo": present_high,
        "n_high_moderate": len(candidates),
        "n_het_noncommon": len(het),
        "frequency_class_counts_het_noncommon": frequency_class_counts(het),
        "frequency_class_counts_high_moderate": frequency_class_counts(candidates),
        "pair_class_counts": result["pair_class_counts"],
        "n_pairs": len(result["pairs"]),
        "phase_summary": result["phase_summary"],
        "clinvar_summary": result["clinvar_summary"],
        "genotype_only_top10": [
            {
                "rank": i,
                "gene": item.pair.gene,
                "overall_score": item.overall,
                **{k: item.scores[k] for k in item.scores},
                "pair_frequency_class": item.pair.pair_frequency_class,
                "phase_class": item.pair.phase_class,
            }
            for i, item in enumerate(result["genotype_genes"][:10], start=1)
        ],
        "phenotype_aware_top10": [
            {
                "rank": i,
                "gene": item.pair.gene,
                "overall_score": item.overall,
                **{k: item.scores[k] for k in item.scores},
                "pair_frequency_class": item.pair.pair_frequency_class,
                "phase_class": item.pair.phase_class,
                "clinvar_category": clinvar_pair_category(item.pair),
            }
            for i, item in enumerate(result["phenotype_genes"][:10], start=1)
        ],
        "overlap_top10": sorted(
            {item.pair.gene for item in result["genotype_genes"][:10]}
            & {item.pair.gene for item in result["phenotype_genes"][:10]}
        ),
        "rank_stability": result["stability"],
        "missing_af_audit_top_phenotype": missing_audits,
        "unannotated_reserve": reserve,
        "ontology_counts": {
            "n_terms": len(ontology.terms),
            "n_diseases": ontology.n_diseases,
            "n_disease_annotations": ontology.n_disease_annotations,
            "n_gene_disease": len(ontology.gene_to_diseases),
        },
    }

    print("frequency_class_counts_het_noncommon", summary["frequency_class_counts_het_noncommon"])
    print("pair_class_counts", summary["pair_class_counts"])
    print("phase_summary", summary["phase_summary"])
    print("clinvar_summary", summary["clinvar_summary"])
    print("genotype_only_top10")
    print(format_gene_rank_report(result["genotype_genes"]))
    print("phenotype_aware_top10")
    print(format_gene_rank_report(result["phenotype_genes"]))
    print("overlap_top10", summary["overlap_top10"])
    print("rank_stability", json.dumps(summary["rank_stability"]))
    if reserve:
        print("unannotated_reserve", json.dumps(reserve))
    if args.summary_json is not None:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
