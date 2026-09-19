"""Stream a local VCF and compute aggregate QC statistics.

Never print individual variant records, coordinates, alleles, genotypes at a
locus, sample identifiers, or PID strings.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, TextIO

import numpy as np

NT = frozenset("ACGTNacgtn")
TRANSITIONS = frozenset({("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")})
AUTOSOMES = frozenset(str(i) for i in range(1, 23))
SEX_MT = {"X": "X", "Y": "Y", "M": "MT", "MT": "MT", "MITO": "MT"}
PATH_RE = re.compile(
    r"(?:(?<![A-Za-z])[A-Za-z]:[\\/]|\\\\|/(?:home|Users|usr|data|opt|mnt)/)[^\s,\"']+"
)
HARD_FILTER_ID_RE = re.compile(
    r"^(QD|FS|MQ|SOR|ReadPosRankSum|MQRankSum|RPRS)",
    re.IGNORECASE,
)
FILEFMT_RE = re.compile(r"##fileformat=VCFv([0-9.]+)", re.IGNORECASE)
CONTIG_RE = re.compile(r"##contig=<ID=([^,>]+)")
FILTER_RE = re.compile(r"##FILTER=<ID=([^,>]+)")
REFERENCE_RE = re.compile(r"##reference=([^\s]+)", re.IGNORECASE)


@dataclass
class NumericSummary:
    n: int = 0
    missing: int = 0
    minimum: float | None = None
    p5: float | None = None
    p25: float | None = None
    median: float | None = None
    mean: float | None = None
    p75: float | None = None
    p95: float | None = None
    maximum: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "missing": self.missing,
            "min": self.minimum,
            "p5": self.p5,
            "p25": self.p25,
            "median": self.median,
            "mean": self.mean,
            "p75": self.p75,
            "p95": self.p95,
            "max": self.maximum,
        }


@dataclass
class StreamStats:
    header_lines: list[str] = field(default_factory=list)
    chrom_header: str = ""
    n_samples: int = 0
    total: int = 0
    pass_n: int = 0
    filtered_n: int = 0
    biallelic: int = 0
    multiallelic: int = 0
    snp: int = 0
    insertion: int = 0
    deletion: int = 0
    mnv_complex: int = 0
    symbolic: int = 0
    spanning_star: int = 0
    chrom_counts: Counter[str] = field(default_factory=Counter)
    chrom_bucket_counts: Counter[str] = field(default_factory=Counter)
    gt_00: int = 0
    gt_01: int = 0
    gt_11: int = 0
    gt_other_diploid: int = 0
    gt_haploid: int = 0
    gt_missing: int = 0
    gt_phased: int = 0
    gt_unphased: int = 0
    het: int = 0
    hom_alt: int = 0
    nonref: int = 0
    het_phased_gt: int = 0
    het_with_pid: int = 0
    het_with_pgt: int = 0
    het_phased_or_pid: int = 0
    has_pgt: int = 0
    has_pid: int = 0
    has_pgt_and_pid: int = 0
    pid_sizes: Counter[str] = field(default_factory=Counter)
    dup_keys: int = 0
    seen_keys: set[tuple[str, int, str, str]] = field(default_factory=set)
    sorted_ok: bool = True
    prev_chrom_idx: int = -1
    prev_pos: int = -1
    contig_order: list[str] = field(default_factory=list)
    contig_index: dict[str, int] = field(default_factory=dict)
    unknown_contig_records: int = 0
    ref_blocks: int = 0
    qual_vals: list[float] = field(default_factory=list)
    qual_missing: int = 0
    dp_vals: list[float] = field(default_factory=list)
    dp_missing: int = 0
    gq_vals: list[float] = field(default_factory=list)
    gq_missing: int = 0
    ab_snv: list[float] = field(default_factory=list)
    ab_indel: list[float] = field(default_factory=list)
    ti_all: int = 0
    tv_all: int = 0
    ti_pass: int = 0
    tv_pass: int = 0
    ti_het_pass: int = 0
    tv_het_pass: int = 0
    ti_hom_pass: int = 0
    tv_hom_pass: int = 0
    nr_dp: int = 0
    nr_gq: int = 0
    nr_ad: int = 0
    nr_pl: int = 0
    nr_pgt: int = 0
    nr_pid: int = 0
    subset_a_nonref: int = 0
    subset_b_pass_nonref: int = 0
    subsets_het: dict[str, int] = field(
        default_factory=lambda: {k: 0 for k in "ABCDEF"}
    )
    subsets_hom: dict[str, int] = field(
        default_factory=lambda: {k: 0 for k in "ABCDEF"}
    )
    prev_chrom: str = ""


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def redact_text(text: str) -> str:
    return PATH_RE.sub("[REDACTED_PATH]", text)


def normalize_chrom(chrom: str) -> str:
    value = chrom.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def chrom_bucket(chrom: str) -> str:
    core = normalize_chrom(chrom).upper()
    if core in SEX_MT:
        return SEX_MT[core]
    if core in AUTOSOMES:
        return core
    if core.isdigit():
        return core
    return "other"


def allele_kind(ref: str, alt: str) -> str:
    if alt == "*":
        return "star"
    if alt.startswith("<") and alt.endswith(">"):
        return "symbolic"
    if not ref or not alt or any(ch not in NT for ch in ref + alt):
        return "complex"
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    if len(alt) > len(ref):
        return "ins"
    if len(alt) < len(ref):
        return "del"
    if len(ref) > 1:
        return "mnv"
    return "complex"


def record_kind(ref: str, alts: list[str]) -> str:
    kinds = {allele_kind(ref, alt) for alt in alts}
    if kinds <= {"snp"}:
        return "snp"
    if kinds <= {"ins"}:
        return "ins"
    if kinds <= {"del"}:
        return "del"
    if kinds <= {"mnv"}:
        return "mnv"
    return "complex"


def is_transition(ref: str, alt: str) -> bool:
    return (ref.upper(), alt.upper()) in TRANSITIONS


def parse_gt(gt: str) -> dict[str, Any]:
    raw = (gt or "").strip()
    if raw == "" or raw == "." or raw in {"./.", ".|."}:
        return {
            "missing": True,
            "phased": "|" in raw,
            "haploid": False,
            "alleles": [],
            "het": False,
            "hom_ref": False,
            "hom_alt": False,
            "nonref": False,
            "label": "missing",
        }
    phased = "|" in raw
    haploid = "/" not in raw and "|" not in raw
    parts = re.split(r"[/|]", raw)
    missing = any(part == "." for part in parts)
    alleles = parts
    numeric = []
    for part in parts:
        if part == ".":
            continue
        try:
            numeric.append(int(part))
        except ValueError:
            numeric.append(-1)
    het = (not missing) and (not haploid) and len(set(numeric)) > 1
    hom_ref = (not missing) and numeric and all(a == 0 for a in numeric)
    hom_alt = (not missing) and numeric and all(a == numeric[0] for a in numeric) and numeric[0] > 0
    if haploid and not missing:
        hom_ref = numeric == [0]
        hom_alt = len(numeric) == 1 and numeric[0] > 0
        het = False
    nonref = (not missing) and (not hom_ref)
    if missing:
        label = "missing"
    elif haploid:
        label = "haploid"
    elif hom_ref:
        label = "0/0"
    elif set(numeric) == {0, 1} and len(numeric) == 2:
        label = "0/1"
    elif numeric == [1, 1]:
        label = "1/1"
    else:
        label = "other_diploid"
    return {
        "missing": missing,
        "phased": phased,
        "haploid": haploid,
        "alleles": alleles,
        "het": het,
        "hom_ref": hom_ref,
        "hom_alt": hom_alt,
        "nonref": nonref,
        "label": label,
    }


def parse_ad(ad: str) -> list[int] | None:
    if not ad or ad == ".":
        return None
    values = []
    for part in ad.split(","):
        if part in {"", "."}:
            return None
        try:
            values.append(int(part))
        except ValueError:
            return None
    return values


def allele_balance(ad: list[int], alleles: list[str]) -> float | None:
    if len(ad) < 2 or len(alleles) != 2:
        return None
    try:
        ref_i = int(alleles[0])
        alt_i = int(alleles[1])
    except ValueError:
        return None
    if ref_i < 0 or alt_i < 0 or ref_i >= len(ad) or alt_i >= len(ad):
        return None
    ref_ad = ad[ref_i]
    alt_ad = ad[alt_i]
    denom = ref_ad + alt_ad
    if denom <= 0:
        return None
    return alt_ad / denom


def numeric_summary(values: list[float], missing: int) -> NumericSummary:
    if not values:
        return NumericSummary(n=0, missing=missing)
    arr = np.asarray(values, dtype=np.float64)
    return NumericSummary(
        n=int(arr.size),
        missing=missing,
        minimum=float(arr.min()),
        p5=float(np.percentile(arr, 5)),
        p25=float(np.percentile(arr, 25)),
        median=float(np.median(arr)),
        mean=float(arr.mean()),
        p75=float(np.percentile(arr, 75)),
        p95=float(np.percentile(arr, 95)),
        maximum=float(arr.max()),
    )


def count_ge(values: list[float], threshold: float) -> int:
    return sum(1 for value in values if value >= threshold)


def ab_band(values: list[float], low: float, high: float) -> int:
    return sum(1 for value in values if low <= value <= high)


def classify_phase_availability(n_het: int, n_het_phased_or_pid: int) -> str:
    """Classify phase-field availability among heterozygous records.

    none: no heterozygous records, or none with phased GT or PID
    very_sparse: 0% < fraction < 1%
    sparse: 1% <= fraction < 10%
    moderate: 10% <= fraction < 50%
    widespread: fraction >= 50%
    """
    if n_het <= 0 or n_het_phased_or_pid <= 0:
        return "none"
    frac = n_het_phased_or_pid / n_het
    if frac < 0.01:
        return "very_sparse"
    if frac < 0.10:
        return "sparse"
    if frac < 0.50:
        return "moderate"
    return "widespread"


def classify_callset(stats: StreamStats) -> tuple[str, list[str]]:
    evidence: list[str] = []
    auto = sum(stats.chrom_bucket_counts[c] for c in AUTOSOMES)
    n_auto_present = sum(1 for c in AUTOSOMES if stats.chrom_bucket_counts[c] > 0)
    evidence.append(f"total_records={stats.total}")
    evidence.append(f"autosomes_with_records={n_auto_present}/22")
    evidence.append(f"autosomal_records={auto}")
    evidence.append(f"nonref={stats.nonref}")
    evidence.append(f"hom_ref_gt_00={stats.gt_00}")
    evidence.append(f"possible_gvcf_ref_blocks={stats.ref_blocks}")
    if stats.nonref:
        evidence.append(f"nonref_to_homref_ratio={stats.nonref / max(stats.gt_00, 1):.4g}")
    if stats.ref_blocks > 0 or stats.gt_00 > stats.nonref:
        evidence.append("hom-ref or END-like records suggest a gVCF-style or all-sites callset")
    # Non-reference count is the WGS/WES discriminator; gVCF 0/0 blocks inflate totals.
    if stats.nonref >= 1_000_000 and n_auto_present >= 20:
        return "likely WGS callset", evidence
    if stats.nonref >= 10_000 and stats.nonref < 500_000 and n_auto_present >= 15:
        return "likely WES callset", evidence
    if stats.nonref < 10_000 or n_auto_present <= 5:
        return "targeted callset", evidence
    return "indeterminate", evidence


def parse_header_metadata(header_lines: list[str]) -> dict[str, Any]:
    fileformat = None
    reference = None
    filters: list[str] = []
    contig_ids: list[str] = []
    sources: list[str] = []
    gatk_lines = 0
    gatk_tool_ids: list[str] = []
    vqsr = False
    hard_filter_hint = False
    joint = False
    gvcf_header_hint = False
    caller = None
    caller_version = None
    commandline_headers_present = False
    for raw in header_lines:
        line = raw.rstrip("\n")
        m = FILEFMT_RE.match(line)
        if m:
            fileformat = m.group(1)
        m = REFERENCE_RE.match(line)
        if m:
            reference = redact_text(m.group(1))
        m = CONTIG_RE.match(line)
        if m:
            contig_ids.append(m.group(1))
        m = FILTER_RE.match(line)
        if m:
            fid = m.group(1)
            filters.append(fid)
            if HARD_FILTER_ID_RE.match(fid):
                hard_filter_hint = True
        if line.startswith("##source="):
            sources.append(redact_text(line.split("=", 1)[1]))
        if "GATKCommandLine" in line or line.startswith("##GATK"):
            gatk_lines += 1
            commandline_headers_present = True
            tool_m = re.search(r"ID=([^,>]+)", line)
            if tool_m:
                tool_id = tool_m.group(1)
                if tool_id not in gatk_tool_ids:
                    gatk_tool_ids.append(tool_id)
        if "VariantRecalibrator" in line or "ApplyVQSR" in line or "VQSR" in line:
            vqsr = True
        if "hard-filter" in line.lower() or "HardFiltered" in line or "QD_filter" in line:
            hard_filter_hint = True
        if any(token in line for token in ("QD_filter", "FS_filter", "MQ_filter", "SOR_filter")):
            hard_filter_hint = True
        if "GenomicsDB" in line or "GenotypeGVCFs" in line or "CombineGVCFs" in line:
            joint = True
        if "GVCFBlock" in line:
            gvcf_header_hint = True
        if "HaplotypeCaller" in line and caller is None:
            caller = "GATK HaplotypeCaller"
        if "Version=" in line and "GATK" in line and caller_version is None:
            vm = re.search(r"Version=([^,\s>]+)", line)
            if vm:
                caller_version = vm.group(1).strip().strip('"').strip("'")
    if caller is None and (gatk_lines or gatk_tool_ids):
        caller = "GATK"
    contig_style = None
    if contig_ids:
        contig_style = "chr" if contig_ids[0].startswith("chr") else "no_chr"
    return {
        "vcf_version": fileformat,
        "reference_declaration": reference,
        "caller": caller,
        "caller_version": caller_version,
        "gatk_commandline_headers": gatk_lines,
        "gatk_commandline_headers_present": commandline_headers_present,
        "gatk_tool_ids": gatk_tool_ids,
        "vqsr_metadata_present": vqsr,
        "hard_filter_metadata_present": hard_filter_hint,
        "joint_genotyping_metadata_present": joint,
        "gvcf_block_metadata_present": gvcf_header_hint,
        "filter_ids": filters,
        "source_headers": sources,
        "contig_count": len(contig_ids),
        "contig_style": contig_style,
        "contig_ids_header": contig_ids,
    }


def _open_text(path: Path) -> TextIO:
    name = path.name.lower()
    if name.endswith(".gz") or name.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def stream_vcf(path: Path) -> tuple[StreamStats, dict[str, Any]]:
    stats = StreamStats()
    header_meta: dict[str, Any] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                stats.header_lines.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                stats.chrom_header = line.rstrip("\n")
                cols = stats.chrom_header.split("\t")
                try:
                    fmt_i = cols.index("FORMAT")
                    stats.n_samples = max(0, len(cols) - fmt_i - 1)
                except ValueError:
                    stats.n_samples = 0
                header_meta = parse_header_metadata(stats.header_lines)
                stats.contig_order = list(header_meta.get("contig_ids_header") or [])
                stats.contig_index = {name: i for i, name in enumerate(stats.contig_order)}
                continue
            if not line.strip() or line.startswith("#"):
                continue
            _consume_record(stats, line)
    if not header_meta:
        header_meta = parse_header_metadata(stats.header_lines)
    return stats, header_meta


def _consume_record(stats: StreamStats, line: str) -> None:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return
    stats.total += 1
    chrom, pos_s, _vid, ref, alt_s, qual_s, filt, info = fields[:8]
    fmt = fields[8] if len(fields) > 8 else ""
    sample = fields[9] if len(fields) > 9 else ""
    try:
        pos = int(pos_s)
    except ValueError:
        pos = -1
    alts = [a for a in alt_s.split(",") if a]
    if len(alts) <= 1:
        stats.biallelic += 1
    else:
        stats.multiallelic += 1
    kind = record_kind(ref, alts)
    if kind == "snp":
        stats.snp += 1
    elif kind == "ins":
        stats.insertion += 1
    elif kind == "del":
        stats.deletion += 1
    else:
        stats.mnv_complex += 1
    if any(allele_kind(ref, a) == "symbolic" for a in alts):
        stats.symbolic += 1
    if any(a == "*" for a in alts):
        stats.spanning_star += 1
    if "END=" in info:
        stats.ref_blocks += 1
    passed = filt in {"PASS", ".", ""}
    if passed:
        stats.pass_n += 1
    else:
        stats.filtered_n += 1
    stats.chrom_counts[chrom] += 1
    stats.chrom_bucket_counts[chrom_bucket(chrom)] += 1
    if stats.contig_index:
        idx = stats.contig_index.get(chrom)
        if idx is None:
            stats.unknown_contig_records += 1
            idx = 10_000 + stats.unknown_contig_records
        if stats.total > 1:
            if idx < stats.prev_chrom_idx or (
                idx == stats.prev_chrom_idx and pos < stats.prev_pos
            ):
                stats.sorted_ok = False
        stats.prev_chrom_idx = idx
    elif stats.total > 1:
        if chrom < stats.prev_chrom or (chrom == stats.prev_chrom and pos < stats.prev_pos):
            stats.sorted_ok = False
    stats.prev_chrom = chrom
    stats.prev_pos = pos
    key = (chrom, pos, ref, alt_s)
    if key in stats.seen_keys:
        stats.dup_keys += 1
    else:
        stats.seen_keys.add(key)
    if qual_s in {"", "."}:
        stats.qual_missing += 1
    else:
        try:
            stats.qual_vals.append(float(qual_s))
        except ValueError:
            stats.qual_missing += 1
    fmt_keys = fmt.split(":") if fmt else []
    sample_vals = sample.split(":") if sample else []
    fmt_map = dict(zip(fmt_keys, sample_vals, strict=False))
    gt_info = parse_gt(fmt_map.get("GT", ""))
    if gt_info["missing"]:
        stats.gt_missing += 1
    elif gt_info["label"] == "0/0":
        stats.gt_00 += 1
    elif gt_info["label"] == "0/1":
        stats.gt_01 += 1
    elif gt_info["label"] == "1/1":
        stats.gt_11 += 1
    elif gt_info["haploid"]:
        stats.gt_haploid += 1
    else:
        stats.gt_other_diploid += 1
    if gt_info["phased"]:
        stats.gt_phased += 1
    elif not gt_info["missing"]:
        stats.gt_unphased += 1
    pgt = fmt_map.get("PGT", "")
    pid = fmt_map.get("PID", "")
    pgt_ok = bool(pgt) and pgt != "."
    pid_ok = bool(pid) and pid != "."
    if gt_info["het"]:
        stats.het += 1
        if gt_info["phased"]:
            stats.het_phased_gt += 1
        if pgt_ok:
            stats.het_with_pgt += 1
        if pid_ok:
            stats.het_with_pid += 1
        if gt_info["phased"] or pid_ok:
            stats.het_phased_or_pid += 1
    if gt_info["hom_alt"]:
        stats.hom_alt += 1
    if gt_info["nonref"]:
        stats.nonref += 1
    if pgt_ok:
        stats.has_pgt += 1
    if pid_ok:
        stats.has_pid += 1
        stats.pid_sizes[pid] += 1
    if pgt_ok and pid_ok:
        stats.has_pgt_and_pid += 1
    dp_raw = fmt_map.get("DP", "")
    gq_raw = fmt_map.get("GQ", "")
    dp_val = None
    gq_val = None
    if dp_raw in {"", "."}:
        stats.dp_missing += 1
    else:
        try:
            dp_val = float(dp_raw)
            stats.dp_vals.append(dp_val)
        except ValueError:
            stats.dp_missing += 1
    if gq_raw in {"", "."}:
        stats.gq_missing += 1
    else:
        try:
            gq_val = float(gq_raw)
            stats.gq_vals.append(gq_val)
        except ValueError:
            stats.gq_missing += 1
    ad = parse_ad(fmt_map.get("AD", ""))
    biallelic = len(alts) == 1
    ab = None
    if (
        gt_info["het"]
        and biallelic
        and ad is not None
        and kind in {"snp", "ins", "del"}
    ):
        ab = allele_balance(ad, gt_info["alleles"])
        if ab is not None:
            if kind == "snp":
                stats.ab_snv.append(ab)
            else:
                stats.ab_indel.append(ab)
    if biallelic and kind == "snp":
        ti = is_transition(ref, alts[0])
        if ti:
            stats.ti_all += 1
        else:
            stats.tv_all += 1
        if passed:
            if ti:
                stats.ti_pass += 1
            else:
                stats.tv_pass += 1
            if gt_info["het"]:
                if ti:
                    stats.ti_het_pass += 1
                else:
                    stats.tv_het_pass += 1
            if gt_info["hom_alt"]:
                if ti:
                    stats.ti_hom_pass += 1
                else:
                    stats.tv_hom_pass += 1
    if gt_info["nonref"]:
        if dp_val is not None:
            stats.nr_dp += 1
        if gq_val is not None:
            stats.nr_gq += 1
        if ad is not None:
            stats.nr_ad += 1
        if fmt_map.get("PL") not in {None, "", "."}:
            stats.nr_pl += 1
        if pgt_ok:
            stats.nr_pgt += 1
        if pid_ok:
            stats.nr_pid += 1
        stats.subset_a_nonref += 1
        if passed:
            stats.subset_b_pass_nonref += 1
        _update_subsets(stats, passed, gt_info, dp_val, gq_val, ab)


def _update_subsets(
    stats: StreamStats,
    passed: bool,
    gt_info: dict[str, Any],
    dp_val: float | None,
    gq_val: float | None,
    ab: float | None,
) -> None:
    target = None
    if gt_info["het"]:
        target = stats.subsets_het
    elif gt_info["hom_alt"]:
        target = stats.subsets_hom
    if target is None:
        return
    target["A"] += 1
    if not passed:
        return
    target["B"] += 1
    if gt_info["het"]:
        target["C"] += 1
        if dp_val is not None and dp_val >= 10:
            target["D"] += 1
            if gq_val is not None and gq_val >= 20:
                target["E"] += 1
                if ab is not None and 0.25 <= ab <= 0.75:
                    target["F"] += 1
    elif gt_info["hom_alt"]:
        target["C"] += 1
        if dp_val is not None and dp_val >= 10:
            target["D"] += 1
            if gq_val is not None and gq_val >= 20:
                target["E"] += 1
                target["F"] += 1


def _ratio(num: int, den: int) -> float | None:
    if den <= 0:
        return None
    return num / den


def _titv(ti: int, tv: int) -> float | None:
    if tv <= 0:
        return None
    return ti / tv


def phase_set_size_distribution(pid_sizes: Counter[str]) -> dict[str, Any]:
    sizes = list(pid_sizes.values())
    if not sizes:
        return {
            "n_phase_sets": 0,
            "max_variants_in_one_set": 0,
            "median_set_size": None,
            "size_histogram": {},
        }
    hist: Counter[str] = Counter()
    for n in sizes:
        if n == 1:
            hist["1"] += 1
        elif n == 2:
            hist["2"] += 1
        elif n <= 5:
            hist["3-5"] += 1
        elif n <= 10:
            hist["6-10"] += 1
        else:
            hist["11+"] += 1
    return {
        "n_phase_sets": len(sizes),
        "max_variants_in_one_set": max(sizes),
        "median_set_size": float(np.median(np.asarray(sizes, dtype=np.float64))),
        "size_histogram": dict(hist),
    }


def scientific_recommendations(
    stats: StreamStats,
    phase_class: str,
    callset: str,
    header_meta: dict[str, Any],
) -> list[str]:
    n_nr = max(stats.nonref, 1)
    dp_frac = stats.nr_dp / n_nr
    gq_frac = stats.nr_gq / n_nr
    ad_frac = stats.nr_ad / n_nr
    recs: list[str] = []
    if dp_frac >= 0.90 and gq_frac >= 0.90:
        recs.append(
            "The VCF is technically suitable to begin annotation/prioritization from the "
            "provided callset; FORMAT DP and GQ are sufficiently populated for QC filters."
        )
    else:
        recs.append(
            "FORMAT DP/GQ completeness is incomplete; treat depth/quality filters as optional "
            "and inspect missingness before discarding sites."
        )
    recs.append(
        "Allele balance is usable for heterozygous QC."
        if ad_frac >= 0.80 and (stats.ab_snv or stats.ab_indel)
        else "Allele-balance filters may be weakly powered because AD is sparse or unused."
    )
    recs.append(
        f"Local read-backed phase fields (phased GT / PGT / PID) are classified as "
        f"{phase_class}; this is not proof of cis/trans and should not replace later "
        f"compound-het logic."
    )
    if stats.multiallelic or stats.insertion or stats.deletion:
        recs.append(
            "Later left-alignment/splitting (`bcftools norm -m -any`) is recommended, but "
            "only against a checksum-verified FASTA that matches this VCF's GRCh38 contig set."
        )
    else:
        recs.append(
            "Multiallelic/indel burden is low; normalization can still wait for a matching FASTA."
        )
    recs.append(
        f"Callset classification from aggregate VCF properties is {callset}; "
        "FASTQs were not inspected."
    )
    recs.append(
        "Do not re-call from FASTQ in the next iteration unless subsequent annotation shows "
        "the supplied VCF is unusable; start from this VCF."
    )
    recs.append("Do not apply population AF filters until an annotation source is chosen.")
    if header_meta.get("caller") or header_meta.get("gatk_tool_ids"):
        recs.append(
            "Upstream caller identity is inferred from header metadata only; command-line "
            "paths were redacted and sample identifiers were not retained."
        )
    return recs


def build_report(
    path: Path,
    stats: StreamStats,
    header_meta: dict[str, Any],
    sha256: str,
    tbi_present: bool,
) -> dict[str, Any]:
    stats.seen_keys.clear()
    n = max(stats.total, 1)
    n_het = stats.het
    n_het_phased_or_pid = stats.het_phased_or_pid
    phase_class = classify_phase_availability(n_het, n_het_phased_or_pid)
    phase_sets = phase_set_size_distribution(stats.pid_sizes)
    stats.pid_sizes.clear()
    callset, evidence = classify_callset(stats)
    ab_all = stats.ab_snv + stats.ab_indel
    auto = sum(stats.chrom_bucket_counts[c] for c in AUTOSOMES)
    header_meta.pop("contig_ids_header", None)
    return {
        "input": {
            "size_bytes": path.stat().st_size,
            "sha256": sha256,
            "tbi_present": tbi_present,
            "sample_count": stats.n_samples,
            "parser": "python_gzip_stream",
            "filename_omitted": True,
        },
        "header": header_meta,
        "records": {
            "total": stats.total,
            "pass": stats.pass_n,
            "filtered": stats.filtered_n,
            "biallelic": stats.biallelic,
            "multiallelic": stats.multiallelic,
            "snp": stats.snp,
            "insertion": stats.insertion,
            "deletion": stats.deletion,
            "mnv_or_complex": stats.mnv_complex,
            "symbolic_allele_records": stats.symbolic,
            "spanning_star_alleles": stats.spanning_star,
            "possible_reference_blocks": stats.ref_blocks,
            "records_per_chrom_bucket": dict(stats.chrom_bucket_counts),
            "autosomal_records": auto,
            "x_records": stats.chrom_bucket_counts.get("X", 0),
            "y_records": stats.chrom_bucket_counts.get("Y", 0),
            "mt_records": stats.chrom_bucket_counts.get("MT", 0),
            "other_contig_records": stats.chrom_bucket_counts.get("other", 0),
        },
        "qual": numeric_summary(stats.qual_vals, stats.qual_missing).to_dict(),
        "genotype": {
            "0/0": stats.gt_00,
            "0/1": stats.gt_01,
            "1/1": stats.gt_11,
            "other_diploid": stats.gt_other_diploid,
            "haploid": stats.gt_haploid,
            "missing": stats.gt_missing,
            "phased": stats.gt_phased,
            "unphased": stats.gt_unphased,
            "heterozygous": stats.het,
            "homozygous_alt": stats.hom_alt,
            "non_reference": stats.nonref,
            "pct_0/0": 100 * stats.gt_00 / n,
            "pct_0/1": 100 * stats.gt_01 / n,
            "pct_1/1": 100 * stats.gt_11 / n,
            "pct_other_diploid": 100 * stats.gt_other_diploid / n,
            "pct_haploid": 100 * stats.gt_haploid / n,
            "pct_missing": 100 * stats.gt_missing / n,
        },
        "dp": {
            **numeric_summary(stats.dp_vals, stats.dp_missing).to_dict(),
            "ge_5": count_ge(stats.dp_vals, 5),
            "ge_10": count_ge(stats.dp_vals, 10),
            "ge_20": count_ge(stats.dp_vals, 20),
            "ge_30": count_ge(stats.dp_vals, 30),
        },
        "gq": {
            **numeric_summary(stats.gq_vals, stats.gq_missing).to_dict(),
            "ge_10": count_ge(stats.gq_vals, 10),
            "ge_20": count_ge(stats.gq_vals, 20),
            "ge_30": count_ge(stats.gq_vals, 30),
            "ge_60": count_ge(stats.gq_vals, 60),
            "ge_90": count_ge(stats.gq_vals, 90),
        },
        "allele_balance": {
            "snv": {
                **numeric_summary(stats.ab_snv, 0).to_dict(),
                "band_0.20_0.80": ab_band(stats.ab_snv, 0.20, 0.80),
                "band_0.25_0.75": ab_band(stats.ab_snv, 0.25, 0.75),
                "band_0.30_0.70": ab_band(stats.ab_snv, 0.30, 0.70),
                "band_0.35_0.65": ab_band(stats.ab_snv, 0.35, 0.65),
            },
            "indel": {
                **numeric_summary(stats.ab_indel, 0).to_dict(),
                "band_0.20_0.80": ab_band(stats.ab_indel, 0.20, 0.80),
                "band_0.25_0.75": ab_band(stats.ab_indel, 0.25, 0.75),
                "band_0.30_0.70": ab_band(stats.ab_indel, 0.30, 0.70),
                "band_0.35_0.65": ab_band(stats.ab_indel, 0.35, 0.65),
            },
            "combined": {
                **numeric_summary(ab_all, 0).to_dict(),
                "band_0.20_0.80": ab_band(ab_all, 0.20, 0.80),
                "band_0.25_0.75": ab_band(ab_all, 0.25, 0.75),
                "band_0.30_0.70": ab_band(ab_all, 0.30, 0.70),
                "band_0.35_0.65": ab_band(ab_all, 0.35, 0.65),
            },
        },
        "phasing": {
            "heterozygous_records": stats.het,
            "heterozygous_phased_gt": stats.het_phased_gt,
            "heterozygous_with_pgt": stats.het_with_pgt,
            "heterozygous_with_pid": stats.het_with_pid,
            "heterozygous_with_phased_gt_or_pid": n_het_phased_or_pid,
            "records_with_pgt": stats.has_pgt,
            "records_with_pid": stats.has_pid,
            "records_with_pgt_and_pid": stats.has_pgt_and_pid,
            "availability": phase_class,
            "availability_definition": (
                "fraction of heterozygous records with phased GT (|) or PID: "
                "none=0; very_sparse<1%; sparse<10%; moderate<50%; widespread>=50%. "
                "PID strings are never emitted."
            ),
            "phase_sets": phase_sets,
            "not_proof_of_cis_trans": True,
        },
        "titv": {
            "all_snps": {
                "ti": stats.ti_all,
                "tv": stats.tv_all,
                "ratio": _titv(stats.ti_all, stats.tv_all),
            },
            "pass_snps": {
                "ti": stats.ti_pass,
                "tv": stats.tv_pass,
                "ratio": _titv(stats.ti_pass, stats.tv_pass),
            },
            "het_pass_snps": {
                "ti": stats.ti_het_pass,
                "tv": stats.tv_het_pass,
                "ratio": _titv(stats.ti_het_pass, stats.tv_het_pass),
            },
            "hom_alt_pass_snps": {
                "ti": stats.ti_hom_pass,
                "tv": stats.tv_hom_pass,
                "ratio": _titv(stats.ti_hom_pass, stats.tv_hom_pass),
            },
        },
        "field_availability_nonref": {
            "n_nonref": stats.nonref,
            "dp": _ratio(stats.nr_dp, stats.nonref),
            "gq": _ratio(stats.nr_gq, stats.nonref),
            "ad": _ratio(stats.nr_ad, stats.nonref),
            "pl": _ratio(stats.nr_pl, stats.nonref),
            "pgt": _ratio(stats.nr_pgt, stats.nonref),
            "pid": _ratio(stats.nr_pid, stats.nonref),
        },
        "candidate_search_size": {
            "A_all_nonref": stats.subset_a_nonref,
            "B_pass_nonref": stats.subset_b_pass_nonref,
            "heterozygous": stats.subsets_het,
            "homozygous_alt": stats.subsets_hom,
            "definitions": {
                "A_global": "all non-reference variants",
                "B_global": "PASS non-reference variants",
                "A": "all non-reference of this zygosity",
                "B": "PASS of this zygosity",
                "C": "PASS het or PASS hom-alt",
                "D": "C and DP>=10",
                "E": "D and GQ>=20",
                "F_het": "E and AB 0.25-0.75",
                "F_hom_alt": "E (AB band not applied to hom-alt)",
            },
        },
        "normalization_readiness": {
            "multiallelic_records": stats.multiallelic,
            "symbolic_allele_records": stats.symbolic,
            "spanning_star_alleles": stats.spanning_star,
            "indel_records": stats.insertion + stats.deletion,
            "duplicate_chrom_pos_ref_alt": stats.dup_keys,
            "records_appear_sorted": stats.sorted_ok,
            "unknown_contig_records": stats.unknown_contig_records,
            "ref_validation_requires_external_fasta": True,
            "recommend_bcftools_norm": bool(
                stats.multiallelic or stats.insertion or stats.deletion or stats.spanning_star
            ),
            "do_not_normalize_without_matching_fasta": True,
        },
        "callset_type": {
            "classification": callset,
            "evidence": evidence,
            "noncoding_not_assessed": True,
        },
        "recommendations": scientific_recommendations(
            stats, phase_class, callset, header_meta
        ),
    }


FORBIDDEN_STDOUT_PATTERNS = (
    re.compile(r"\bchr\d+\s+\d{4,}"),
    re.compile(r"\b[1-9XYM]{1,2}\t\d{4,}"),
)


def text_contains_locus_like_token(text: str) -> bool:
    return any(pattern.search(text) for pattern in FORBIDDEN_STDOUT_PATTERNS)


def _fmt_num_summary(block: dict[str, Any]) -> str:
    parts = []
    for key in ("n", "missing", "min", "p5", "p25", "median", "mean", "p75", "p95", "max"):
        value = block.get(key)
        if isinstance(value, float):
            parts.append(f"{key}={value:.4g}")
        else:
            parts.append(f"{key}={value}")
    return " ".join(parts)


def format_terminal_summary(report: dict[str, Any]) -> str:
    rec = report["records"]
    gt = report["genotype"]
    ph = report["phasing"]
    cand = report["candidate_search_size"]
    ab = report["allele_balance"]["combined"]
    miss = report["field_availability_nonref"]
    lines = [
        "VCF profile (aggregates only; original VCF not modified)",
        f"compressed_size_bytes={report['input']['size_bytes']}",
        f"sha256={report['input']['sha256']}",
        f"tbi_present={report['input']['tbi_present']}",
        f"vcf_version={report['header'].get('vcf_version')}",
        f"sample_count={report['input']['sample_count']}",
        f"contig_style={report['header'].get('contig_style')}",
        f"contig_count={report['header'].get('contig_count')}",
        f"caller={report['header'].get('caller')}",
        f"caller_version={report['header'].get('caller_version')}",
        f"gatk_tool_ids={report['header'].get('gatk_tool_ids')}",
        f"vqsr_metadata={report['header'].get('vqsr_metadata_present')}",
        f"hard_filter_metadata={report['header'].get('hard_filter_metadata_present')}",
        f"joint_genotyping_metadata={report['header'].get('joint_genotyping_metadata_present')}",
        f"total_records={rec['total']}",
        f"pass={rec['pass']} filtered={rec['filtered']}",
        f"biallelic={rec['biallelic']} multiallelic={rec['multiallelic']}",
        f"snp={rec['snp']} ins={rec['insertion']} del={rec['deletion']} "
        f"mnv_complex={rec['mnv_or_complex']}",
        f"symbolic={rec['symbolic_allele_records']} spanning_star={rec['spanning_star_alleles']}",
        f"chrom_buckets={rec['records_per_chrom_bucket']}",
        f"autosomal={rec['autosomal_records']} X={rec['x_records']} Y={rec['y_records']} "
        f"MT={rec['mt_records']}",
        f"gt_0/0={gt['0/0']} gt_0/1={gt['0/1']} gt_1/1={gt['1/1']} "
        f"other_diploid={gt['other_diploid']} haploid={gt['haploid']} missing={gt['missing']}",
        f"gt_phased={gt['phased']} gt_unphased={gt['unphased']}",
        f"het={gt['heterozygous']} hom_alt={gt['homozygous_alt']} nonref={gt['non_reference']}",
        f"qual={_fmt_num_summary(report['qual'])}",
        (
            f"dp={_fmt_num_summary(report['dp'])} ge5={report['dp']['ge_5']} "
            f"ge10={report['dp']['ge_10']}"
        ),
        f"dp_ge20={report['dp']['ge_20']} dp_ge30={report['dp']['ge_30']}",
        (
            f"gq={_fmt_num_summary(report['gq'])} ge10={report['gq']['ge_10']} "
            f"ge20={report['gq']['ge_20']}"
        ),
        (
            f"gq_ge30={report['gq']['ge_30']} ge60={report['gq']['ge_60']} "
            f"ge90={report['gq']['ge_90']}"
        ),
        f"ab_n={ab['n']} ab_median={ab['median']} ab_p5={ab['p5']} ab_p25={ab['p25']} "
        f"ab_p75={ab['p75']} ab_p95={ab['p95']}",
        f"ab_bands 0.20-0.80={ab['band_0.20_0.80']} 0.25-0.75={ab['band_0.25_0.75']} "
        f"0.30-0.70={ab['band_0.30_0.70']} 0.35-0.65={ab['band_0.35_0.65']}",
        f"phase_availability={ph['availability']}",
        f"het_phased_gt={ph['heterozygous_phased_gt']} het_pgt={ph['heterozygous_with_pgt']} "
        f"het_pid={ph['heterozygous_with_pid']}",
        f"phase_sets={ph['phase_sets']}",
        (
            f"titv_all={report['titv']['all_snps']['ratio']} "
            f"titv_pass={report['titv']['pass_snps']['ratio']}"
        ),
        f"field_availability_nonref dp={miss['dp']} gq={miss['gq']} ad={miss['ad']} "
        f"pl={miss['pl']} pgt={miss['pgt']} pid={miss['pid']}",
        (
            f"candidate_A_all_nonref={cand['A_all_nonref']} "
            f"candidate_B_pass_nonref={cand['B_pass_nonref']}"
        ),
        f"candidate_het={cand['heterozygous']}",
        f"candidate_hom_alt={cand['homozygous_alt']}",
        f"callset_type={report['callset_type']['classification']}",
        f"records_appear_sorted={report['normalization_readiness']['records_appear_sorted']}",
        f"duplicate_keys={report['normalization_readiness']['duplicate_chrom_pos_ref_alt']}",
    ]
    return "\n".join(lines) + "\n"


def render_markdown(report: dict[str, Any]) -> str:
    rec = report["records"]
    header = report["header"]
    ph = report["phasing"]
    norm = report["normalization_readiness"]
    cs = report["callset_type"]
    cand = report["candidate_search_size"]
    recs = report.get("recommendations") or []
    return "\n".join(
        [
            "# Local VCF profile (aggregates only)",
            "",
            "The official VCF was not modified. No variant records, coordinates, alleles, "
            "sample identifiers, or PID strings are listed here.",
            "",
            "## Input",
            f"- compressed size bytes: {report['input']['size_bytes']}",
            f"- SHA-256: `{report['input']['sha256']}`",
            f"- TBI present: {report['input']['tbi_present']}",
            f"- sample count: {report['input']['sample_count']}",
            f"- parser: {report['input']['parser']}",
            "",
            "## Header / caller",
            f"- VCF version: {header.get('vcf_version')}",
            f"- reference declaration: {header.get('reference_declaration')}",
            f"- caller: {header.get('caller')}",
            f"- caller version: {header.get('caller_version')}",
            f"- GATK command-line headers: {header.get('gatk_commandline_headers')}",
            f"- GATK tool IDs: {header.get('gatk_tool_ids')}",
            f"- VQSR metadata: {header.get('vqsr_metadata_present')}",
            f"- hard-filter metadata: {header.get('hard_filter_metadata_present')}",
            f"- joint-genotyping metadata: {header.get('joint_genotyping_metadata_present')}",
            f"- gVCF-block metadata: {header.get('gvcf_block_metadata_present')}",
            f"- FILTER IDs: {', '.join(header.get('filter_ids') or [])}",
            f"- contig style: {header.get('contig_style')} (n={header.get('contig_count')})",
            "- Local filesystem paths in header text were redacted.",
            "- Command-line strings were not copied.",
            "",
            "## Records",
            f"- total {rec['total']}; PASS {rec['pass']}; filtered {rec['filtered']}",
            f"- biallelic {rec['biallelic']}; multiallelic {rec['multiallelic']}",
            f"- SNP {rec['snp']}; insertion {rec['insertion']}; deletion {rec['deletion']}; "
            f"MNV/complex {rec['mnv_or_complex']}",
            (
                f"- symbolic {rec['symbolic_allele_records']}; "
                f"spanning * {rec['spanning_star_alleles']}"
            ),
            f"- autosomal {rec['autosomal_records']}; X {rec['x_records']}; Y {rec['y_records']}; "
            f"MT {rec['mt_records']}",
            f"- per-chrom buckets: `{rec['records_per_chrom_bucket']}`",
            "",
            "## QUAL / DP / GQ",
            f"- QUAL: `{report['qual']}`",
            f"- DP: `{report['dp']}`",
            f"- GQ: `{report['gq']}`",
            "",
            "## Allele balance (het biallelic SNV+indel with usable AD)",
            f"- combined: `{report['allele_balance']['combined']}`",
            f"- SNV: `{report['allele_balance']['snv']}`",
            f"- indel: `{report['allele_balance']['indel']}`",
            "",
            "## Field availability among non-reference genotypes",
            f"- `{report['field_availability_nonref']}`",
            "",
            "## Callset type",
            f"- classification: **{cs['classification']}**",
            f"- evidence: {cs['evidence']}",
            "- noncoding content was not assessed (no annotation yet).",
            "",
            "## Phasing",
            f"- heterozygous records: {ph['heterozygous_records']}",
            f"- heterozygous with phased GT: {ph['heterozygous_phased_gt']}",
            (
                f"- heterozygous with PGT: {ph['heterozygous_with_pgt']}; "
                f"PID: {ph['heterozygous_with_pid']}"
            ),
            f"- records with PGT: {ph['records_with_pgt']}; PID: {ph['records_with_pid']}; both: "
            f"{ph['records_with_pgt_and_pid']}",
            f"- availability: **{ph['availability']}**",
            f"- definition: {ph['availability_definition']}",
            f"- phase-set size histogram: `{ph['phase_sets']}`",
            "- This is not proof of cis/trans configuration.",
            "",
            "## Normalization readiness",
            f"- duplicates: {norm['duplicate_chrom_pos_ref_alt']}",
            f"- appears sorted: {norm['records_appear_sorted']}",
            f"- unknown contig records vs header: {norm['unknown_contig_records']}",
            f"- recommend later `bcftools norm -m -any` against a **matching** GRCh38 FASTA: "
            f"{norm['recommend_bcftools_norm']}",
            "- Do not normalize against an unverified FASTA.",
            "",
            "## Candidate-search size (descriptive, no AF filter)",
            f"- A all nonref: {cand['A_all_nonref']}",
            f"- B PASS nonref: {cand['B_pass_nonref']}",
            f"- heterozygous A–F: `{cand['heterozygous']}`",
            f"- homozygous-alt A–F: `{cand['homozygous_alt']}`",
            "",
            "## Ti/Tv (biallelic SNPs)",
            f"- all: {report['titv']['all_snps']}",
            f"- PASS: {report['titv']['pass_snps']}",
            f"- het PASS: {report['titv']['het_pass_snps']}",
            f"- hom-alt PASS: {report['titv']['hom_alt_pass_snps']}",
            "",
            "## Recommendations for the next iteration",
            *[f"- {item}" for item in recs],
            "",
        ]
    )


def find_proband_vcf(data_dir: Path) -> Path:
    if not data_dir.is_dir():
        raise FileNotFoundError("dataset directory not found")
    candidates: list[Path] = []
    for path in data_dir.rglob("*"):
        if not path.is_file():
            continue
        rel_parts = path.relative_to(data_dir).parts
        if any(part.startswith(".") for part in rel_parts[:-1]):
            continue
        name = path.name.lower()
        if name.endswith(".vcf.gz") or name.endswith(".vcf.bgz"):
            candidates.append(path)
    if not candidates:
        raise FileNotFoundError("no compressed VCF found in the dataset directory")
    candidates.sort(key=lambda item: item.stat().st_size, reverse=True)
    return candidates[0]


def profile_vcf(vcf_path: Path, output_dir: Path | None = None) -> dict[str, Any]:
    tbi = Path(str(vcf_path) + ".tbi")
    if not tbi.is_file():
        tbi = vcf_path.with_suffix(vcf_path.suffix + ".tbi")
    sha = sha256_file(vcf_path)
    stats, header_meta = stream_vcf(vcf_path)
    report = build_report(vcf_path, stats, header_meta, sha, tbi.is_file())
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "vcf_profile.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        (output_dir / "vcf_profile.md").write_text(
            render_markdown(report) + "\n", encoding="utf-8"
        )
    return report


def parse_args(argv: list[str] | None = None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Profile a local Track 1 VCF and emit aggregate statistics only."
    )
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--vcf", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    import sys

    from mva_track1.download import is_inside_repo
    from mva_track1.private_outputs import (
        data_root,
        ensure_private_layout,
        vcf_profile_root,
    )

    args = parse_args(argv)
    ensure_private_layout()
    data_dir = args.data_dir or data_root()
    try:
        vcf_path = args.vcf or find_proband_vcf(data_dir)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if is_inside_repo(vcf_path):
        print("Refusing to profile a VCF inside the git repository.", file=sys.stderr)
        return 2
    output_dir = args.output_dir or vcf_profile_root()
    if is_inside_repo(output_dir):
        print("Refusing to write profile reports inside the git repository.", file=sys.stderr)
        return 2
    report = profile_vcf(vcf_path, output_dir)
    print(format_terminal_summary(report), end="")
    print("Wrote aggregate JSON/Markdown reports to the local work directory.")
    return 0
