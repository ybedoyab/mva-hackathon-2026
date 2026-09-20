from __future__ import annotations

from pathlib import Path

from hpo_fixtures import write_mini_hpoa, write_mini_ontology
from mva_track1.candidates import (
    VariantCandidate,
    build_pairs,
    frequency_class,
    pair_frequency_class,
    public_variant_view,
)
from mva_track1.hpo_ontology import load_gene_disease_map, load_phenotype_hpoa, parse_obo
from mva_track1.phenotype_score import score_gene_phenotype
from mva_track1.ranking import (
    build_ranked_lanes,
    format_gene_rank_report,
    functional_score,
    pair_component_scores,
    phase_multiplier,
)

SECRET_POS = "888777666"
SECRET_REF = "SECRETA"
SECRET_ALT = "SECRETG"
SECRET_CHROM = "SECRETCHR"


def _var(**kwargs) -> VariantCandidate:
    defaults = dict(
        variant_id="x",
        gene="GENEA",
        transcript="NM_1",
        consequence="missense_variant",
        impact="MODERATE",
        max_af=0.00005,
        gnomadg_af=0.00005,
        frequency_class="OBSERVED_ULTRARARE",
        clinvar="missing",
        sift="",
        polyphen="",
        mane="NM_1",
        canonical="YES",
        dp=40,
        gq=50,
        ab=0.48,
        filt="PASS",
        pid="",
        pgt="",
        existing_variant="",
        existing_status="missing",
        zygosity="het",
        chrom="1",
        pos=100,
        ref="A",
        alt="G",
    )
    defaults.update(kwargs)
    return VariantCandidate(**defaults)  # type: ignore[arg-type]


def test_frequency_and_pair_classes() -> None:
    assert frequency_class(None) == "AF_MISSING"
    assert frequency_class(0.0) == "ZERO_REPORTED"
    assert frequency_class(0.00005) == "OBSERVED_ULTRARARE"
    assert frequency_class(0.0005) == "OBSERVED_RARE"
    assert frequency_class(0.01) == "COMMON"
    assert pair_frequency_class("OBSERVED_ULTRARARE", "OBSERVED_RARE") == "observed/observed"
    assert pair_frequency_class("ZERO_REPORTED", "OBSERVED_RARE") == "observed/zero"
    assert pair_frequency_class("ZERO_REPORTED", "ZERO_REPORTED") == "zero/zero"
    assert pair_frequency_class("OBSERVED_RARE", "AF_MISSING") == "observed/missing"
    assert pair_frequency_class("ZERO_REPORTED", "AF_MISSING") == "zero/missing"
    assert pair_frequency_class("AF_MISSING", "AF_MISSING") == "missing/missing"


def test_phase_cis_versus_unknown() -> None:
    left = _var(variant_id="a", pid="P1", pgt="0|1", pos=1)
    right = _var(variant_id="b", pid="P1", pgt="0|1", pos=2)
    cis = build_pairs([left, right])[0]
    assert cis.phase_class == "likely_cis"
    assert cis.phase_status == "same_phase"
    assert phase_multiplier(cis.phase_class) < 0.2
    unknown = build_pairs(
        [_var(variant_id="c", gene="GENEB", pos=3), _var(variant_id="d", gene="GENEB", pos=4)]
    )[0]
    assert unknown.phase_class == "phase_unknown"
    assert unknown.phase_status == "insufficient"


def test_functional_score_ordering() -> None:
    splice = functional_score(_var(consequence="splice_acceptor_variant", impact="HIGH"))
    frameshift = functional_score(_var(consequence="frameshift_variant", impact="HIGH"))
    missense = functional_score(_var(consequence="missense_variant", impact="MODERATE"))
    inframe = functional_score(_var(consequence="inframe_deletion", impact="MODERATE"))
    assert splice > frameshift > missense >= inframe


def test_phenotype_score_and_negative_separate(tmp_path: Path) -> None:
    ontology = parse_obo(write_mini_ontology(tmp_path))
    load_phenotype_hpoa(ontology, write_mini_hpoa(tmp_path))
    g2d = tmp_path / "genes_to_disease.txt"
    g2d.write_text(
        "ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n"
        "1\tGENEA\tmendelian\tOMIM:1\tHPO\n"
        "2\tGENEB\tmendelian\tOMIM:4\tHPO\n",
        encoding="utf-8",
    )
    load_gene_disease_map(ontology, g2d)
    scored = score_gene_phenotype("GENEA", ontology, ["HP:0001250"], ["HP:0000256"])
    other = score_gene_phenotype("GENEB", ontology, ["HP:0001250"], ["HP:0000256"])
    assert scored.positive_phenotype_score > other.positive_phenotype_score
    assert scored.negative_contradiction_score >= 0
    assert "positive_phenotype_score" in scored.__dataclass_fields__
    assert "negative_contradiction_score" in scored.__dataclass_fields__


def test_rank_reproducibility_and_sensitivity(tmp_path: Path) -> None:
    ontology = parse_obo(write_mini_ontology(tmp_path))
    load_phenotype_hpoa(ontology, write_mini_hpoa(tmp_path))
    g2d = tmp_path / "genes_to_disease.txt"
    g2d.write_text(
        "ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n"
        "1\tGENEA\tmendelian\tOMIM:1\tHPO\n"
        "2\tGENEB\tmendelian\tOMIM:4\tHPO\n",
        encoding="utf-8",
    )
    load_gene_disease_map(ontology, g2d)
    candidates = [
        _var(variant_id="a1", gene="GENEA", consequence="stop_gained", impact="HIGH", max_af=0.0,
             frequency_class="ZERO_REPORTED"),
        _var(variant_id="a2", gene="GENEA", consequence="missense_variant", max_af=0.0002,
             frequency_class="OBSERVED_RARE", pos=101),
        _var(variant_id="b1", gene="GENEB", consequence="missense_variant", max_af=None,
             frequency_class="AF_MISSING", filt="LowQual", dp=8, gq=10),
        _var(variant_id="b2", gene="GENEB", consequence="missense_variant", max_af=0.0004,
             frequency_class="OBSERVED_RARE", pos=201),
    ]
    first = build_ranked_lanes(candidates, ontology, ["HP:0001250"], [], [])
    second = build_ranked_lanes(candidates, ontology, ["HP:0001250"], [], [])
    assert [item.pair.gene for item in first["phenotype_genes"]] == [
        item.pair.gene for item in second["phenotype_genes"]
    ]
    assert first["sensitivity"]["S1"][0].pair.gene == first["phenotype_genes"][0].pair.gene
    assert "S4" in first["stability"][first["phenotype_genes"][0].pair.gene]
    report = format_gene_rank_report(first["phenotype_genes"])
    assert SECRET_POS not in report
    assert SECRET_REF not in report
    assert SECRET_ALT not in report
    planted = _var(
        variant_id="z1",
        gene="GENEA",
        chrom=SECRET_CHROM,
        pos=int(SECRET_POS),
        ref=SECRET_REF,
        alt=SECRET_ALT,
    )
    planted2 = _var(
        variant_id="z2",
        gene="GENEA",
        chrom=SECRET_CHROM,
        pos=int(SECRET_POS) + 1,
        ref=SECRET_REF,
        alt=SECRET_ALT,
        consequence="frameshift_variant",
        impact="HIGH",
    )
    ranked = build_ranked_lanes([planted, planted2], ontology, ["HP:0001250"], [], [])
    text = format_gene_rank_report(ranked["genotype_genes"])
    assert SECRET_CHROM not in text
    assert SECRET_POS not in text
    assert SECRET_REF not in text
    assert SECRET_ALT not in text
    view = public_variant_view(planted)
    assert "chrom" not in view
    assert "ref" not in view
    assert "alt" not in view
    scores = pair_component_scores(first["pairs"][0])
    assert set(scores) >= {
        "frequency_score",
        "functional_score",
        "quality_score",
        "clinvar_score",
        "phase_score",
        "genotype_score",
    }
