from __future__ import annotations

import csv
from pathlib import Path

import pytest

from mva_track1.submission import (
    SubmissionError,
    empty_template_csv,
    load_submission,
    parse_raw_rows,
    write_empty_template,
)


def test_empty_template_has_official_columns_and_no_predictions() -> None:
    text = empty_template_csv()
    header = text.strip().split(",")
    assert header[:10] == [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
    ]
    assert "finding_type" in header
    assert text.count("\n") == 1


def test_write_empty_template(tmp_path: Path) -> None:
    path = tmp_path / "template.csv"
    write_empty_template(path)
    assert path.read_text(encoding="utf-8") == empty_template_csv()


def test_valid_compound_het_and_optional_finding_type() -> None:
    rows = [
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr1",
            "pos_1": "1000",
            "ref_1": "a",
            "alt_1": "g",
            "chrom_2": "chr1",
            "pos_2": "2000",
            "ref_2": "c",
            "alt_2": "t",
            "epcr": "0.9",
            "finding_type": "primary",
        }
    ]
    parsed, was_sorted = parse_raw_rows(rows)
    assert was_sorted
    assert parsed[0].rank == 1
    assert parsed[0].variants == frozenset({("chr1", 1000, "A", "G"), ("chr1", 2000, "C", "T")})


def test_optional_second_variant_may_be_blank() -> None:
    rows = [
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr2",
            "pos_1": "10",
            "ref_1": "T",
            "alt_1": "A",
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": "0.4",
        }
    ]
    parsed, _ = parse_raw_rows(rows)
    assert parsed[0].variants == frozenset({("chr2", 10, "T", "A")})


def test_partial_second_variant_is_rejected() -> None:
    rows = [
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr2",
            "pos_1": "10",
            "ref_1": "T",
            "alt_1": "A",
            "chrom_2": "chr3",
            "pos_2": "",
            "ref_2": "G",
            "alt_2": "C",
            "epcr": "0.4",
        }
    ]
    with pytest.raises(SubmissionError, match="variant 2"):
        parse_raw_rows(rows)


def test_epcr_bounds_and_sort() -> None:
    with pytest.raises(SubmissionError, match="out of range"):
        parse_raw_rows(
            [
                {
                    "proband_id": "PROBAND01",
                    "chrom_1": "chr1",
                    "pos_1": "1",
                    "ref_1": "A",
                    "alt_1": "T",
                    "chrom_2": "",
                    "pos_2": "",
                    "ref_2": "",
                    "alt_2": "",
                    "epcr": "0",
                }
            ]
        )
    rows = [
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr1",
            "pos_1": "1",
            "ref_1": "A",
            "alt_1": "T",
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": "0.2",
        },
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr1",
            "pos_1": "2",
            "ref_1": "A",
            "alt_1": "C",
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": "0.8",
        },
    ]
    parsed, was_sorted = parse_raw_rows(rows)
    assert was_sorted is False
    assert parsed[0].epcr == 0.8
    assert parsed[0].rank == 1
    assert parsed[1].rank == 2


def test_max_ten_rows() -> None:
    rows = [
        {
            "proband_id": "PROBAND01",
            "chrom_1": "chr1",
            "pos_1": str(index + 1),
            "ref_1": "A",
            "alt_1": "T",
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": "0.1",
        }
        for index in range(11)
    ]
    with pytest.raises(SubmissionError, match="max is 10"):
        parse_raw_rows(rows)


def test_load_submission_csv(tmp_path: Path) -> None:
    path = tmp_path / "pred.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "proband_id",
                "chrom_1",
                "pos_1",
                "ref_1",
                "alt_1",
                "chrom_2",
                "pos_2",
                "ref_2",
                "alt_2",
                "epcr",
                "finding_type",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "proband_id": "PROBAND01",
                "chrom_1": "chr1",
                "pos_1": "5",
                "ref_1": "G",
                "alt_1": "A",
                "chrom_2": "",
                "pos_2": "",
                "ref_2": "",
                "alt_2": "",
                "epcr": "1",
                "finding_type": "secondary",
            }
        )
    parsed, was_sorted = load_submission(path)
    assert was_sorted
    assert parsed[0].finding_type == "secondary"
    assert parsed[0].epcr == 1.0
