from __future__ import annotations

from pathlib import Path

from mva_track1.vep_run import summarize_vep_warnings, vep_args


def test_vep_args_are_offline_and_not_everything() -> None:
    args = vep_args(
        input_vcf=Path("D:/mva-hackathon-2026-work/annotation/in.vcf.gz"),
        output_vcf=Path("D:/mva-hackathon-2026-work/annotation/out.vcf.gz"),
        cache_dir=Path("D:/mva-vep-data"),
        warning_file=Path("D:/mva-hackathon-2026-work/logs/warn.txt"),
        stats_file=Path("D:/mva-hackathon-2026-work/logs/stats.html"),
        forks=4,
        synonyms=Path("D:/mva-vep-data/homo_sapiens/116_GRCh38/chr_synonyms.txt"),
    )
    joined = " ".join(args)
    assert "--offline" in args
    assert "--cache_version" in args
    assert "116" in args
    assert "--flag_pick_allele_gene" in args
    assert "--everything" not in args
    assert "--hgvs" not in joined
    assert "--database" not in args
    assert "--pick" not in args or args[args.index("--flag_pick_allele_gene") :]


def test_warning_summary_hides_contig_names(tmp_path: Path) -> None:
    secret = "SECRETCONTIGXYZ"
    warn = tmp_path / "warn.txt"
    err = tmp_path / "err.txt"
    warn.write_text(
        f"WARNING: Chromosome '{secret}' not found in cache\n"
        f"WARNING: Chromosome '{secret}' not found in cache\n"
        "WARNING: something deprecated happened\n",
        encoding="utf-8",
    )
    err.write_text("command_argc=12\n", encoding="utf-8")
    summary = summarize_vep_warnings(warn, err)
    text = str(summary)
    assert secret not in text
    assert summary["unrecognized_contig_count"] == 1
    assert summary["warning_categories"]["chrom_not_found"] == 2
    assert summary["does_not_list_contig_names"] is True
