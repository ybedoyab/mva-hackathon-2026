from __future__ import annotations

from mva_track1.reference_match import (
    CLASS_EXACT,
    CLASS_INCOMPATIBLE,
    CLASS_NEAR_EXACT,
    ContigDict,
    compare_contig_dicts,
    format_terminal_summary,
    select_reference,
    strip_chr,
)


def _dict(name: str, lengths: dict[str, int]) -> ContigDict:
    return ContigDict(name=name, lengths=lengths, source="synthetic")


def test_exact_match() -> None:
    vcf = _dict("vcf", {"1": 100, "2": 200})
    ref = _dict("ref", {"1": 100, "2": 200})
    row = compare_contig_dicts(vcf, ref)
    assert row["classification"] == CLASS_EXACT
    assert row["percent_vcf_name_length_match"] == 100.0
    assert row["vcf_contigs_missing_from_ref"] == 0


def test_chr_prefix_is_near_exact_not_exact() -> None:
    vcf = _dict("vcf", {str(i): 100 + i for i in range(1, 23)} | {"X": 10, "Y": 11, "M": 16569})
    ref = _dict(
        "ref",
        {f"chr{i}": 100 + i for i in range(1, 23)} | {"chrX": 10, "chrY": 11, "chrM": 16569},
    )
    identity = compare_contig_dicts(vcf, ref)
    stripped = compare_contig_dicts(vcf, ref, transform=strip_chr, transform_name="strip_chr")
    assert identity["classification"] == CLASS_INCOMPATIBLE
    assert stripped["classification"] == CLASS_NEAR_EXACT
    assert stripped["equal_contig_count"] is True
    selected = select_reference([identity, stripped])
    assert selected is not None
    assert selected["transform"] == "strip_chr"


def test_length_mismatch_is_incompatible() -> None:
    vcf = _dict("vcf", {"1": 248956422, "2": 10})
    ref = _dict("ref", {"1": 248956422, "2": 99})
    row = compare_contig_dicts(vcf, ref)
    assert row["classification"] == CLASS_INCOMPATIBLE
    assert row["exact_name_and_length_matches"] == 1


def test_summary_does_not_list_contigs() -> None:
    vcf = _dict("vcf", {"SECRETCONTIG": 12345678})
    ref = _dict("ref", {"chrSECRETCONTIG": 12345678})
    row = compare_contig_dicts(vcf, ref, transform=strip_chr, transform_name="strip_chr")
    text = format_terminal_summary([row], row)
    assert "SECRETCONTIG" not in text
    assert "12345678" not in text
    assert "classification" in text or "class=" in text
