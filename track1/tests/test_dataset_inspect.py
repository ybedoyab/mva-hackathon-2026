from __future__ import annotations

import gzip
import json
from pathlib import Path

from mva_track1.dataset_inspect import format_summary, inspect_dataset
from mva_track1.dataset_inspect import main as inspect_main

SECRET_SAMPLE = "SECRET_SAMPLE_NAME"
SECRET_POS = "999999"
SECRET_REF = "AAA"
SECRET_ALT = "TTT"
SECRET_GT = "1|0"
SECRET_CLINICAL = "SECRET_CLINICAL_TEXT"
SECRET_SEQ = "ACGTSECRETSEQ"


def _write_vcf(path: Path) -> None:
    body = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=100>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SECRET_SAMPLE}\n"
        f"chr1\t{SECRET_POS}\t.\t{SECRET_REF}\t{SECRET_ALT}\t50\tPASS\tDP=9\tGT\t{SECRET_GT}\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(body)


def _write_fastq(path: Path) -> None:
    path.write_bytes(f"@synth\n{SECRET_SEQ}\n+\nIIIIIIIIIIIIII\n".encode("ascii"))


def _synthetic_dataset(tmp_path: Path) -> Path:
    data_dir = tmp_path / "dataset"
    data_dir.mkdir()
    _write_vcf(data_dir / "synth.vcf.gz")
    (data_dir / "synth.vcf.gz.tbi").write_bytes(b"fake-index")
    (data_dir / "notes.docx").write_bytes(SECRET_CLINICAL.encode("utf-8"))
    for lane in (1, 2, 3, 4):
        _write_fastq(data_dir / f"synth_L00{lane}_R1_001.fastq.gz")
        _write_fastq(data_dir / f"synth_L00{lane}_R2_001.fastq.gz")
    return data_dir


def _assert_no_leaks(text: str) -> None:
    secrets = (
        SECRET_SAMPLE,
        SECRET_POS,
        SECRET_REF,
        SECRET_ALT,
        SECRET_GT,
        SECRET_CLINICAL,
        SECRET_SEQ,
    )
    for secret in secrets:
        assert secret not in text


def test_safe_inspector_summary(tmp_path: Path) -> None:
    data_dir = _synthetic_dataset(tmp_path)
    manifest = inspect_dataset(data_dir)
    assert manifest.vcf_present
    assert manifest.vcf_index_present
    assert manifest.phenotype_document_present
    assert manifest.fastq_pair_count == 4
    assert manifest.vcf_header.readable
    assert manifest.vcf_header.sample_count == 1
    assert manifest.vcf_header.contig_style == "chr"
    assert "DP" in manifest.vcf_header.info_ids
    assert "GT" in manifest.vcf_header.format_ids
    assert manifest.genome_build_expected == "GRCh38"
    summary = format_summary(data_dir, manifest)
    _assert_no_leaks(summary)
    assert "VCF sample count: 1" in summary


def test_json_manifest_redacts_sample_and_records(tmp_path: Path) -> None:
    data_dir = _synthetic_dataset(tmp_path)
    json_path = tmp_path / "safe.json"
    assert inspect_main(["--data-dir", str(data_dir), "--json-output", str(json_path)]) == 0
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    dumped = json.dumps(payload)
    _assert_no_leaks(dumped)
    assert payload["vcf_present"] is True
    assert payload["vcf_header"]["sample_count"] == 1
    assert "sample_names" not in payload["vcf_header"]


def test_refuses_to_emit_vcf_body_fields(tmp_path: Path) -> None:
    data_dir = _synthetic_dataset(tmp_path)
    json_path = tmp_path / "out.json"
    inspect_main(["--data-dir", str(data_dir), "--json-output", str(json_path)])
    text = json_path.read_text(encoding="utf-8")
    assert SECRET_POS not in text
    assert SECRET_REF not in text
    assert SECRET_ALT not in text
    assert SECRET_GT not in text
