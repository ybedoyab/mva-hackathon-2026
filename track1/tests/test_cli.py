from __future__ import annotations

import gzip
from pathlib import Path

import pytest
from typer.testing import CliRunner

from mva_track1.cli import app


def test_help_exits_zero() -> None:
    result = CliRunner().invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "preflight" in result.output.lower()
    assert "info" in result.output.lower()
    assert "profile-vcf" in result.output.lower()


def test_info_exits_zero() -> None:
    result = CliRunner().invoke(app, ["info"])
    assert result.exit_code == 0
    assert "mva-track1" in result.output
    compact = " ".join(result.output.lower().split())
    assert "mva-track1" in compact
    assert "scientific pipeline" in compact
    assert "not implemented" in compact


def test_preflight_without_data(monkeypatch: pytest.MonkeyPatch) -> None:
    secret = "hf_" + ("q" * 32)
    monkeypatch.setenv("HF_TOKEN", secret)
    result = CliRunner().invoke(app, ["preflight"])
    assert result.exit_code == 0
    assert "preflight passed" in result.output.lower()
    assert "GRCh38" in result.output
    assert secret not in result.output
    assert "HF_TOKEN" in result.output


def test_preflight_with_synthetic_data_dir(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    vcf = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHIDDENNAME\n"
        "1\t12345\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n"
    )
    with gzip.open(data_dir / "x.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write(vcf)
    result = CliRunner().invoke(app, ["preflight", "--data-dir", str(data_dir)])
    assert result.exit_code == 0
    assert "VCF: present" in result.output
    assert "HIDDENNAME" not in result.output
    assert "12345" not in result.output


def test_profile_vcf_cli_aggregates_only(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    vcf_body = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHIDDENNAME\n"
        "1\t12345\t.\tA\tC\t40\tPASS\t.\tGT:AD:DP:GQ\t0/1:5,5:10:30\n"
    )
    with gzip.open(data_dir / "x.vcf.gz", "wt", encoding="utf-8") as handle:
        handle.write(vcf_body)
    out = tmp_path / "work" / "vcf_profile"
    result = CliRunner().invoke(
        app,
        ["profile-vcf", "--data-dir", str(data_dir), "--output-dir", str(out)],
    )
    assert result.exit_code == 0, result.output
    assert "total_records=1" in result.output
    assert "HIDDENNAME" not in result.output
    assert "12345" not in result.output
    assert (out / "vcf_profile.json").is_file()
    payload = (out / "vcf_profile.json").read_text(encoding="utf-8")
    assert "HIDDENNAME" not in payload
    assert "12345" not in payload
