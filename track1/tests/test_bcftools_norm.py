from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from mva_track1.vcf_norm import docker_available, faidx, run_bcftools

SYNTH_FASTA = """>1
NNNNNNNNNNACGTACGTACNNNNNNNNNN
>2
NNNNNNNNNNTTTTTTTTTTNNNNNNNNNN
"""

SYNTH_VCF = """##fileformat=VCFv4.2
##contig=<ID=1,length=32>
##contig=<ID=2,length=32>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="PGT">
##FORMAT=<ID=PID,Number=1,Type=String,Description="PID">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SYNTH
1	11	.	A	C	50	PASS	.	GT:AD:PGT:PID	0/1:10,10:0|1:SYNTHPID
1	15	.	ACGT	A	50	PASS	.	GT:AD:PGT:PID	0/1:8,8:0|1:SYNTHPID
1	16	.	C	CA	50	PASS	.	GT:AD:PGT:PID	0/1:7,7:.:.
1	21	.	N	A,C	50	PASS	.	GT:AD:PGT:PID	0/1:5,4,1:.:.
1	22	.	N	A,T	50	PASS	.	GT:AD:PGT:PID	0/2:5,1,4:.:.
1	23	.	N	G,C	50	PASS	.	GT:AD:PGT:PID	1/2:1,6,5:.:.
"""


@pytest.mark.skipif(not docker_available(), reason="docker daemon not available")
def test_bcftools_norm_split_preserves_or_transforms_format(tmp_path: Path) -> None:
    fasta = tmp_path / "ref.fa"
    vcf = tmp_path / "in.vcf.gz"
    out = tmp_path / "out.vcf.gz"
    fasta.write_text(SYNTH_FASTA, encoding="utf-8")
    faidx(fasta, tmp_path / "faidx.log")
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write(SYNTH_VCF)
    log = tmp_path / "norm.log"
    proc = run_bcftools(
        [
            "bcftools",
            "norm",
            "-f",
            str(fasta),
            "-m",
            "-any",
            "--old-rec-tag",
            "ORIG",
            "--keep-sum",
            "AD",
            "-Oz",
            "-o",
            str(out),
            str(vcf),
        ],
        mounts=[tmp_path],
        log_path=log,
        check=False,
    )
    assert proc.returncode == 0, log.read_text(encoding="utf-8", errors="replace")
    assert out.is_file()
    text = gzip.open(out, "rt", encoding="utf-8").read()
    records = [line.split("\t") for line in text.splitlines() if line and not line.startswith("#")]
    assert len(records) == 9
    # Biallelic SNP/indel keep GT/AD/PGT/PID; ORIG is added only when the record changes.
    snp = next(row for row in records if row[1] == "11")
    assert snp[9].startswith("0/1:10,10:0|1:SYNTHPID")
    assert "ORIG=" not in snp[7]
    # 0/1 split: carried ALT stays 0/1; unused ALT becomes 0/0. PGT/PID are copied.
    split_01 = [row for row in records if row[1] == "21"]
    gts = {row[4]: row[9].split(":")[0] for row in split_01}
    assert gts == {"A": "0/1", "C": "0/0"}
    assert all("ORIG=1|21|N|A,C|" in row[7] for row in split_01)
    # 0/2 split: ALT2 is the carried allele.
    split_02 = [row for row in records if row[1] == "22"]
    gts02 = {row[4]: row[9].split(":")[0] for row in split_02}
    assert gts02 == {"A": "0/0", "T": "0/1"}
    # 1/2 split: both alts remain non-reference, but PGT is copied verbatim and is
    # not scientifically allele-specific after decomposition.
    split_12 = [row for row in records if row[1] == "23"]
    gts12 = {row[4]: row[9].split(":")[0] for row in split_12}
    assert gts12 == {"G": "1/0", "C": "0/1"}
    assert all("ORIG=1|23|N|G,C|" in row[7] for row in split_12)
    # PGT/PID may be copied or dropped after 1/2 split; they are not allele-specific.
    assert "PGT" in text
    assert "PID" in text
