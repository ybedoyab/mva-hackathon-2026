from __future__ import annotations

from pathlib import Path

from mva_track1.unannotated_reserve import classify_contig, profile_unannotated_reserve


def _vcf(path: Path, rows: list[str], samples: bool = False) -> None:
    header = "##fileformat=VCFv4.2\n"
    if samples:
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    else:
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")


def test_unannotated_reserve_aggregates(tmp_path: Path) -> None:
    assert classify_contig("chrUn_decoy") == "decoy"
    assert classify_contig("chr1_KI270762v1_alt") == "alternate"
    input_vcf = tmp_path / "in.vcf"
    vep_vcf = tmp_path / "vep.vcf"
    gt_vcf = tmp_path / "gt.vcf"
    _vcf(
        input_vcf,
        [
            "1\t10\t.\tA\tG\t.\tPASS\t.",
            "chrUn_decoy\t20\t.\tA\tC\t.\tPASS\tAF=0.1",
            "chrUn_decoy\t30\t.\tG\tT\t.\tLowQual\t.",
        ],
    )
    _vcf(vep_vcf, ["1\t10\t.\tA\tG\t.\tPASS\tCSQ=."])
    _vcf(
        gt_vcf,
        [
            "1\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:40",
            "chrUn_decoy\t20\t.\tA\tC\t.\tPASS\t.\tGT:DP:GQ\t0/1:12:21",
            "chrUn_decoy\t30\t.\tG\tT\t.\tLowQual\t.\tGT:DP:GQ\t1/1:8:10",
        ],
        samples=True,
    )
    report = profile_unannotated_reserve(input_vcf, vep_vcf, gt_vcf)
    assert report["rescue_lane"] is True
    assert report["phenotype_ranking_applied"] is False
    assert report["unrecognized_contig_count"] == 1
    assert report["variant_count"] == 2
    assert report["pass_count"] == 1
    assert report["heterozygous_count"] == 1
    assert report["homozygous_alt_count"] == 1
    assert report["contig_categories"]["decoy"] == 1
    assert "chrUn_decoy" not in str(report)
