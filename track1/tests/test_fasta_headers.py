from __future__ import annotations

from pathlib import Path

from mva_track1.fasta_headers import strip_chr_fasta_headers


def test_strip_chr_headers_preserves_sequence(tmp_path: Path) -> None:
    src = tmp_path / "ref.fa"
    src.write_text(">chr1 AC:xx\nACGT\n>chrM\nN\n>2\nTT\n", encoding="utf-8")
    dest = tmp_path / "ref.stripped.fa"
    stats = strip_chr_fasta_headers(src, dest)
    text = dest.read_text(encoding="utf-8")
    assert stats["records"] == 3
    assert stats["headers_renamed"] == 2
    assert ">1 AC:xx\nACGT\n" in text
    assert ">M\nN\n" in text
    assert ">2\nTT\n" in text
    assert "chr" not in text
