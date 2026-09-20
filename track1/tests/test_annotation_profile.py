from __future__ import annotations

import gzip
import json
from pathlib import Path

from mva_track1.annotation_profile import (
    af_bucket,
    clin_sig_category,
    format_terminal_summary,
    load_original_genotypes,
    load_phase_by_locus,
    most_severe_consequence,
    parse_csq_format,
    parse_gt_zygosity,
    parse_orig,
    phase_pair_status,
    stats_to_report,
    stream_vep_vcf,
)
from mva_track1.vcf_profile import parse_gt

SECRET_GENE = "SECRETGENEABC"
SECRET_SAMPLE = "SECRETSAMPLEXYZ"
SECRET_POS = "888777666"
SECRET_PID = "SECRET_PID_VALUE"


def _vep_vcf() -> str:
    csq_header = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP CSQ. '
        'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|MAX_AF|gnomADe_AF|gnomADg_AF|CLIN_SIG|SIFT|'
        'PolyPhen|MANE_SELECT|CANONICAL|GENE_PHENO">\n'
    )
    header = (
        "##fileformat=VCFv4.2\n"
        + csq_header
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
    )
    # Two rare hets in same gene, plus a common missense, plus missing AF.
    rows = [
        (
            f"1\t{SECRET_POS}\t.\tA\tG\t.\tPASS\t"
            f"CSQ=G|stop_gained|HIGH|{SECRET_GENE}|ENSG00000000001|0.00005|0.00005|.|"
            f"pathogenic|.|.|NM_1|YES|1\tGT:PGT:PID\t0/1:0|1:{SECRET_PID}"
        ),
        (
            f"1\t{int(SECRET_POS)+10}\t.\tC\tT\t.\tPASS\t"
            f"CSQ=T|missense_variant|MODERATE|{SECRET_GENE}|ENSG00000000001|0.0002|0.0002|.|"
            f".|deleterious(0.01)|probably_damaging(0.9)|.|YES|\tGT:PGT:PID\t0/1:1|0:{SECRET_PID}"
        ),
        (
            "1\t100\t.\tG\tA\t.\tPASS\t"
            f"CSQ=A|intron_variant|MODIFIER|{SECRET_GENE}|ENSG00000000001|0.2|0.2|0.2|"
            f".|.|.|.|.|\tGT\t0/1"
        ),
        (
            "1\t200\t.\tT\tA\t.\tLowQual\t"
            f"CSQ=A|missense_variant|MODERATE|{SECRET_GENE}|ENSG00000000001|.|.|.|"
            f"uncertain_significance|.|.|.|.|\tGT\t1/1"
        ),
    ]
    return header + "\n".join(rows) + "\n"


def test_csq_and_severity_and_buckets() -> None:
    assert parse_csq_format(
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="x Format: Allele|Consequence">'
    ) == ["Allele", "Consequence"]
    assert most_severe_consequence(["intron_variant&missense_variant", "synonymous_variant"]) == (
        "missense_variant"
    )
    assert af_bucket(None) == "missing"
    assert af_bucket(0.0) == "eq_0"
    assert af_bucket(0.00005) == "gt0_le_0.0001"
    assert af_bucket(0.0005) == "gt_0.0001_le_0.001"
    assert af_bucket(0.005) == "gt_0.001_le_0.01"
    assert af_bucket(0.2) == "gt_0.01"
    assert clin_sig_category("Pathogenic") == "pathogenic"
    assert clin_sig_category("likely_pathogenic") == "likely_pathogenic"
    assert clin_sig_category("pathogenic/likely_pathogenic") == "pathogenic/likely_pathogenic"
    assert clin_sig_category("uncertain_significance") == "VUS"
    assert clin_sig_category("conflicting_interpretations_of_pathogenicity") == "conflicting"
    assert clin_sig_category("benign") == "benign"
    assert clin_sig_category("likely_benign") == "likely_benign"
    assert clin_sig_category("") == "missing"


def test_phase_and_gt() -> None:
    assert parse_gt_zygosity("0/1") == "het"
    assert parse_gt_zygosity("1/1") == "hom_alt"
    assert parse_gt_zygosity("0/0") == "hom_ref"
    assert parse_gt("1/2")["het"] is True
    assert phase_pair_status("set1", "set1", "0|1", "0|1") == "same_phase"
    assert phase_pair_status("set1", "set1", "0|1", "1|0") == "opposite_phase"
    assert phase_pair_status("set1", "set2", "0|1", "0|1") == "insufficient"
    assert phase_pair_status(".", "set1", "0|1", "0|1") == "insufficient"


def test_stream_aggregates_without_leaks(tmp_path: Path) -> None:
    vep = tmp_path / "vep.vcf.gz"
    orig = tmp_path / "orig.vcf.gz"
    body = _vep_vcf()
    with gzip.open(vep, "wt", encoding="utf-8") as handle:
        handle.write(body)
    with gzip.open(orig, "wt", encoding="utf-8") as handle:
        handle.write(body)
    stats = stream_vep_vcf(vep, load_original_genotypes(orig), load_phase_by_locus(orig))
    report = stats_to_report(stats)
    assert report["input_variants"] == 4
    assert report["annotated"] == 4
    assert report["unique_variant_most_severe"]["stop_gained"] == 1
    assert report["unique_variant_most_severe"]["missense_variant"] == 2
    assert report["max_af_buckets"]["missing"] == 1
    assert report["max_af_buckets"]["gt0_le_0.0001"] == 1
    assert report["clin_sig"]["pathogenic"] == 1
    assert report["clin_sig"]["VUS"] == 1
    assert report["compound_het"]["genes_ge2_rare_het"] == 1
    assert report["compound_het"]["candidate_pairs"] == 1
    assert report["compound_het"]["pairs_phase_opposite"] == 1
    text = format_terminal_summary(report) + json.dumps(report)
    for secret in (SECRET_GENE, SECRET_SAMPLE, SECRET_POS, SECRET_PID):
        assert secret not in text
    assert "888777666" not in text
    assert report["csq_entry_counts"]["stop_gained"] == 1
    assert report["waterfall"]["pass_observed_maxaf_le_0.001"] == 2
    assert report["waterfall"]["pass_maxaf_le_0.001_or_missing"] == 2
    assert report["missing_af_kept_distinct_from_rare"] is True


def test_orig_and_unique_vs_transcript(tmp_path: Path) -> None:
    assert parse_orig("1|21|N|A,C|1") == {
        "chrom": "1",
        "pos": 21,
        "ref": "N",
        "alts": "A,C",
        "allele_index": "1",
    }
    csq_header = (
        '##INFO=<ID=CSQ,Number=.,Type=String,Description="x '
        'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|MAX_AF|gnomADe_AF|gnomADg_AF|'
        'CLIN_SIG|SIFT|PolyPhen|MANE_SELECT|CANONICAL|GENE_PHENO">\n'
    )
    chrom_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    header = "##fileformat=VCFv4.2\n" + csq_header + chrom_line
    # One variant, two transcript consequences: entries=2, unique most-severe=1.
    body = (
        header
        + "1\t50\t.\tA\tG\t.\tPASS\t"
        + "CSQ=G|missense_variant&intron_variant|MODERATE|GENEA|ENSG1|0.00001|.|.|.|"
        + "deleterious(0.01)|probably_damaging(0.9)|NM_1|YES|1,"
        + "G|intron_variant|MODIFIER|GENEA|ENSG1|0.00001|.|.|.|.|.|.|.|\n"
        + "1\t21\t.\tN\tC\t.\tPASS\tORIG=1|21|N|A,C|2;"
        + "CSQ=C|missense_variant|MODERATE|GENEB|ENSG2|0.0002|.|.|.|.|.|.|.|\n"
    )
    vep = tmp_path / "vep.vcf"
    vep.write_text(body, encoding="utf-8")
    orig = tmp_path / "orig.vcf"
    orig.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
        "1\t21\t.\tN\tA,C\t.\tPASS\t.\tGT:PGT:PID\t0/1:0|1:ORIGPID\n",
        encoding="utf-8",
    )
    norm = tmp_path / "norm.vcf"
    norm.write_text(
        "##fileformat=VCFv4.2\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
        "1\t21\t.\tN\tA\t.\tPASS\tORIG=1|21|N|A,C|1\tGT:PGT:PID\t0/1:0|1:ORIGPID\n"
        "1\t21\t.\tN\tC\t.\tPASS\tORIG=1|21|N|A,C|2\tGT:PGT:PID\t0/0:.:.\n",
        encoding="utf-8",
    )
    stats = stream_vep_vcf(vep, load_original_genotypes(norm), load_phase_by_locus(orig))
    report = stats_to_report(stats)
    assert stats.csq_entries == 3
    assert report["unique_variant_most_severe"]["missense_variant"] == 2
    assert report["waterfall"].get("observed_rare_HIGH_MODERATE_het", 0) == 0
    # The split unused ALT is 0/0, so it must not enter compound-het hets.
    assert report["compound_het"]["candidate_pairs"] == 0
    text = format_terminal_summary(report)
    assert SECRET_SAMPLE not in text
    assert "ORIGPID" not in text
    assert "GENEA" not in text
    assert "GENEB" not in text
