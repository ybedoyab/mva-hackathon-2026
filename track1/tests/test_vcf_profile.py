from __future__ import annotations

import gzip
import json
from pathlib import Path

from mva_track1.private_outputs import ensure_private_layout, future_private_table_path
from mva_track1.vcf_profile import (
    allele_kind,
    classify_phase_availability,
    format_terminal_summary,
    parse_gt,
    profile_vcf,
    record_kind,
    redact_text,
    render_markdown,
    text_contains_locus_like_token,
)

SECRET_SAMPLE = "SECRET_SAMPLE_XYZ"
SECRET_PID = "SECRET_PID_COORD"
SECRET_PID_2 = "SECRET_PID_TWO"
SECRET_PATH_USER = "secretuser"
SECRET_POS = "246801357"
SECRET_REF_LONG = "GGGGCCCC"
SECRET_ALT_INS = "ATTTG"


def _header() -> str:
    return (
        "##fileformat=VCFv4.2\n"
        f"##reference=file:///home/{SECRET_PATH_USER}/GRCh38.fa\n"
        "##source=HaplotypeCaller\n"
        "##GATKCommandLine=<ID=HaplotypeCaller,Version=4.2.6.1,"
        f'CommandLine="HaplotypeCaller -R /home/{SECRET_PATH_USER}/GRCh38.fa '
        f'-I /data/{SECRET_SAMPLE}.bam">\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
        "##contig=<ID=1,length=100000000>\n"
        "##contig=<ID=2,length=100000000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
        '##FORMAT=<ID=PGT,Number=1,Type=String,Description="PGT">\n'
        '##FORMAT=<ID=PID,Number=1,Type=String,Description="PID">\n'
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
    )


def _records_sorted() -> str:
    fmt = "GT:AD:DP:GQ:PGT:PID:PL"
    rows = [
        f"1\t{SECRET_POS}\t.\tA\tG\t80\tPASS\t.\t{fmt}\t"
        f"0/1:10,10:20:30:0|1:{SECRET_PID}:30,0,40",
        "1\t246801457\t.\tA\tC\t70\tPASS\t.\tGT:AD:DP:GQ:PL\t1/0:8,8:12:25:0,30,40",
        f"1\t246801557\t.\tA\tG\t90\tPASS\t.\t{fmt}\t"
        f"0|1:20,20:40:90:0|1:{SECRET_PID}:0,90,200",
        f"1\t246801657\t.\tC\tT\t50\tPASS\t.\tGT:AD:DP:GQ:PID:PL\t"
        f"1|0:5,15:15:40:{SECRET_PID_2}:10,0,80",
        "1\t246801757\t.\tA\tT\t99\tPASS\t.\tGT:AD:DP:GQ:PL\t1/1:0,30:30:99:99,20,0",
        "1\t246801857\t.\tA\tC,T\t40\tPASS\t.\tGT:AD:DP:GQ:PL\t1/2:6,6,6:18:22:0,10,20",
        f"1\t246801957\t.\tA\t{SECRET_ALT_INS}\t30\tPASS\t.\tGT:AD:DP:GQ:PL\t0/1:5,6:11:21:20,0,30",
        f"1\t246802057\t.\t{SECRET_REF_LONG}\tG\t20\tPASS\t.\tGT:AD:DP:GQ:PL\t0/1:4,5:9:19:10,0,20",
        "1\t246802157\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:GQ\t./.:.:.:.",
        "1\t246802257\t.\tA\tG\t10\tPASS\t.\tGT:DP:GQ\t1:8:12",
        "1\t246802357\t.\tA\tG\t5\tPASS\t.\tGT:AD:DP:GQ\t0/0:10,0:10:20",
        "1\t246802457\t.\tA\tC\t1\tLowQual\t.\tGT:AD:DP:GQ\t0/1:2,2:4:5",
        "1\t246802557\t.\tA\t*\t20\tPASS\t.\tGT:AD:DP:GQ\t0/1:5,5:10:20",
        "1\t246802657\t.\tA\t<DEL>\t20\tPASS\t.\tGT:DP:GQ\t0/1:10:20",
        "2\t100\t.\tG\tA\t40\tPASS\t.\tGT:AD:DP:GQ\t0/1:12,12:24:50",
    ]
    return "\n".join(rows) + "\n"


def _write_vcf(path: Path, body: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(body)


def _assert_no_leaks(text: str) -> None:
    secrets = (
        SECRET_SAMPLE,
        SECRET_PID,
        SECRET_PID_2,
        SECRET_PATH_USER,
        SECRET_POS,
        SECRET_REF_LONG,
        SECRET_ALT_INS,
        "246801457",
        "0/1:10,10",
        "HaplotypeCaller -R",
    )
    for secret in secrets:
        assert secret not in text, secret
    assert not text_contains_locus_like_token(text)


def test_allele_and_record_kinds() -> None:
    assert allele_kind("A", "G") == "snp"
    assert allele_kind("A", "AT") == "ins"
    assert allele_kind("AT", "A") == "del"
    assert allele_kind("AG", "TC") == "mnv"
    assert allele_kind("A", "*") == "star"
    assert allele_kind("A", "<DEL>") == "symbolic"
    assert record_kind("A", ["G"]) == "snp"
    assert record_kind("A", ["AT"]) == "ins"
    assert record_kind("AT", ["A"]) == "del"
    assert record_kind("A", ["C", "T"]) == "snp"
    assert record_kind("A", ["G", "AT"]) == "complex"


def test_parse_gt_cases() -> None:
    assert parse_gt("0/1")["het"] is True
    assert parse_gt("0/1")["label"] == "0/1"
    assert parse_gt("1/0")["het"] is True
    assert parse_gt("1/0")["label"] == "0/1"
    assert parse_gt("0|1")["phased"] is True
    assert parse_gt("0|1")["het"] is True
    assert parse_gt("1|0")["phased"] is True
    assert parse_gt("1/1")["hom_alt"] is True
    assert parse_gt("1/1")["label"] == "1/1"
    assert parse_gt("1/2")["het"] is True
    assert parse_gt("1/2")["label"] == "other_diploid"
    assert parse_gt("./.")["missing"] is True
    assert parse_gt(".")["missing"] is True
    assert parse_gt("1")["haploid"] is True
    assert parse_gt("1")["hom_alt"] is True
    assert parse_gt("0/0")["hom_ref"] is True
    assert parse_gt("0/0")["nonref"] is False


def test_phase_classification_thresholds() -> None:
    assert classify_phase_availability(0, 0) == "none"
    assert classify_phase_availability(100, 0) == "none"
    assert classify_phase_availability(1000, 5) == "very_sparse"
    assert classify_phase_availability(100, 5) == "sparse"
    assert classify_phase_availability(100, 20) == "moderate"
    assert classify_phase_availability(100, 50) == "widespread"


def test_redact_paths() -> None:
    unix = redact_text("reference=file:///home/secretuser/GRCh38.fa")
    assert "secretuser" not in unix
    assert "[REDACTED_PATH]" in unix
    assert unix.startswith("reference=file:")
    windows = redact_text("reference=file:///C:/Users/secretuser/GRCh38.fa")
    assert "secretuser" not in windows
    assert "[REDACTED_PATH]" in windows


def test_profile_synthetic_vcf_aggregates(tmp_path: Path) -> None:
    vcf = tmp_path / "synth.vcf.gz"
    _write_vcf(vcf, _header() + _records_sorted())
    (tmp_path / "synth.vcf.gz.tbi").write_bytes(b"fake")
    out = tmp_path / "profile"
    report = profile_vcf(vcf, out)
    rec = report["records"]
    gt = report["genotype"]
    ph = report["phasing"]

    assert rec["total"] == 15
    assert rec["pass"] == 14
    assert rec["filtered"] == 1
    assert rec["biallelic"] == 14
    assert rec["multiallelic"] == 1
    assert rec["snp"] == 11
    assert rec["insertion"] == 1
    assert rec["deletion"] == 1
    assert rec["mnv_or_complex"] == 2
    assert rec["symbolic_allele_records"] == 1
    assert rec["spanning_star_alleles"] == 1
    assert rec["records_per_chrom_bucket"]["1"] == 14
    assert rec["records_per_chrom_bucket"]["2"] == 1
    assert rec["autosomal_records"] == 15

    assert gt["0/1"] == 10
    assert gt["1/1"] == 1
    assert gt["other_diploid"] == 1
    assert gt["haploid"] == 1
    assert gt["missing"] == 1
    assert gt["0/0"] == 1
    assert gt["heterozygous"] == 11
    assert gt["homozygous_alt"] == 2
    assert gt["non_reference"] == 13
    assert gt["phased"] == 2
    assert gt["unphased"] == 12

    assert ph["heterozygous_phased_gt"] == 2
    assert ph["heterozygous_with_pgt"] == 2
    assert ph["heterozygous_with_pid"] == 3
    assert ph["records_with_pgt"] == 2
    assert ph["records_with_pid"] == 3
    assert ph["records_with_pgt_and_pid"] == 2
    assert ph["phase_sets"]["n_phase_sets"] == 2
    assert ph["phase_sets"]["max_variants_in_one_set"] == 2
    assert ph["availability"] == "moderate"
    assert ph["not_proof_of_cis_trans"] is True

    assert report["dp"]["n"] == 14
    assert report["dp"]["missing"] == 1
    assert report["dp"]["ge_10"] == 11
    assert report["gq"]["ge_20"] == 11
    ab = report["allele_balance"]["combined"]
    assert ab["n"] == 8
    assert abs(ab["median"] - 0.5) < 1e-9
    assert ab["band_0.25_0.75"] == 8

    assert report["titv"]["all_snps"]["ti"] == 7
    assert report["titv"]["all_snps"]["tv"] == 3
    assert report["titv"]["pass_snps"]["ti"] == 7
    assert report["titv"]["pass_snps"]["tv"] == 2
    assert report["titv"]["het_pass_snps"]["ti"] == 4
    assert report["titv"]["het_pass_snps"]["tv"] == 1
    assert report["titv"]["hom_alt_pass_snps"]["ti"] == 1
    assert report["titv"]["hom_alt_pass_snps"]["tv"] == 1

    cand = report["candidate_search_size"]
    assert cand["A_all_nonref"] == 13
    assert cand["B_pass_nonref"] == 12
    assert cand["heterozygous"]["A"] == 11
    assert cand["heterozygous"]["B"] == 10
    assert cand["heterozygous"]["C"] == 10
    assert cand["heterozygous"]["D"] == 9
    assert cand["heterozygous"]["E"] == 9
    assert cand["heterozygous"]["F"] == 6
    assert cand["homozygous_alt"]["A"] == 2
    assert cand["homozygous_alt"]["B"] == 2
    assert cand["homozygous_alt"]["D"] == 1
    assert cand["homozygous_alt"]["F"] == 1

    assert report["normalization_readiness"]["records_appear_sorted"] is True
    assert report["normalization_readiness"]["duplicate_chrom_pos_ref_alt"] == 0
    assert report["header"]["vcf_version"] == "4.2"
    assert report["header"]["caller"] == "GATK HaplotypeCaller"
    assert report["header"]["caller_version"] == "4.2.6.1"
    assert "HaplotypeCaller" in report["header"]["gatk_tool_ids"]
    assert report["input"]["tbi_present"] is True
    assert report["input"]["sample_count"] == 1
    assert report["input"]["filename_omitted"] is True

    dumped = json.dumps(report)
    md = (out / "vcf_profile.md").read_text(encoding="utf-8")
    term = format_terminal_summary(report)
    _assert_no_leaks(dumped)
    _assert_no_leaks(md)
    _assert_no_leaks(term)
    _assert_no_leaks(render_markdown(report))


def test_duplicate_and_unsorted_detection(tmp_path: Path) -> None:
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=100>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
    )
    dup_body = header + "1\t10\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n" * 2
    unsorted_body = header + (
        "1\t50\t.\tA\tG\t1\tPASS\t.\tGT\t0/1\n1\t10\t.\tA\tC\t1\tPASS\t.\tGT\t0/1\n"
    )
    dup = tmp_path / "dup.vcf.gz"
    unsorted = tmp_path / "unsorted.vcf.gz"
    _write_vcf(dup, dup_body)
    _write_vcf(unsorted, unsorted_body)
    dup_report = profile_vcf(dup)
    unsorted_report = profile_vcf(unsorted)
    assert dup_report["normalization_readiness"]["duplicate_chrom_pos_ref_alt"] == 1
    assert unsorted_report["normalization_readiness"]["records_appear_sorted"] is False
    _assert_no_leaks(json.dumps(dup_report))
    _assert_no_leaks(json.dumps(unsorted_report))
    _assert_no_leaks(format_terminal_summary(dup_report))


def test_stream_does_not_keep_sample_name(tmp_path: Path) -> None:
    vcf = tmp_path / "synth.vcf.gz"
    _write_vcf(vcf, _header() + _records_sorted())
    report = profile_vcf(vcf)
    dumped = json.dumps(report) + format_terminal_summary(report) + render_markdown(report)
    assert SECRET_SAMPLE not in dumped
    assert SECRET_PID not in dumped
    assert SECRET_PID_2 not in dumped


def test_private_layout_does_not_write_candidate_tables(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("MVA_OUTPUT_DIR", str(tmp_path / "work"))
    root = ensure_private_layout()
    csv_path = future_private_table_path("candidates")
    assert csv_path.parent == root
    assert not csv_path.exists()
    assert list(root.glob("*.csv")) == []
    assert (root / "README.txt").is_file()
