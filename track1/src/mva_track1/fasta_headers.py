"""Rewrite FASTA headers by stripping a leading 'chr' without changing sequence."""

from __future__ import annotations

import gzip
from pathlib import Path


def strip_chr_fasta_headers(source: Path, dest: Path) -> dict[str, int]:
    """Write dest FASTA with 'chr' removed from record names. Sequence bytes unchanged."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    n_records = 0
    n_renamed = 0
    src_open = gzip.open if source.name.lower().endswith(".gz") else open
    dst_open = gzip.open if dest.name.lower().endswith(".gz") else open
    with src_open(source, "rt", encoding="utf-8", errors="replace") as inf, dst_open(
        dest, "wt", encoding="utf-8"
    ) as out:
        for line in inf:
            if line.startswith(">"):
                n_records += 1
                body = line[1:]
                if body.startswith("chr") or body.startswith("CHR"):
                    body = body[3:]
                    n_renamed += 1
                out.write(">" + body)
            else:
                out.write(line)
    return {"records": n_records, "headers_renamed": n_renamed}
