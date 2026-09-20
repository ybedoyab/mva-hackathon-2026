"""Private variant candidate model and compound-het pair construction.

Patient coordinates and alleles may exist in memory and private local files.
They must not be printed. Tests must use synthetic variants only.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

from mva_track1.annotation_profile import (
    _lookup_genotype,
    _open_text,
    clin_sig_category,
    most_severe_consequence,
    parse_af,
    parse_csq_format,
    parse_gt_zygosity,
    parse_info,
    parse_orig,
    phase_pair_status,
)

AF_OBSERVED_ULTRARARE = "OBSERVED_ULTRARARE"
AF_OBSERVED_RARE = "OBSERVED_RARE"
AF_ZERO_REPORTED = "ZERO_REPORTED"
AF_MISSING = "AF_MISSING"
AF_COMMON = "COMMON"

OBSERVED_CLASSES = {AF_OBSERVED_ULTRARARE, AF_OBSERVED_RARE}
PAIR_CLASS_ORDER = (
    "observed/observed",
    "observed/zero",
    "zero/zero",
    "observed/missing",
    "zero/missing",
    "missing/missing",
)

HIGH_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
}


def frequency_class(max_af: float | None) -> str:
    if max_af is None:
        return AF_MISSING
    if max_af == 0:
        return AF_ZERO_REPORTED
    if max_af <= 0.0001:
        return AF_OBSERVED_ULTRARARE
    if max_af <= 0.001:
        return AF_OBSERVED_RARE
    return AF_COMMON


def coarse_frequency(freq_class: str) -> str:
    if freq_class in OBSERVED_CLASSES:
        return "observed"
    if freq_class == AF_ZERO_REPORTED:
        return "zero"
    if freq_class == AF_MISSING:
        return "missing"
    return "common"


def pair_frequency_class(class_a: str, class_b: str) -> str:
    order = {"observed": 0, "zero": 1, "missing": 2, "common": 3}
    left, right = coarse_frequency(class_a), coarse_frequency(class_b)
    if order.get(left, 9) > order.get(right, 9):
        left, right = right, left
    return f"{left}/{right}"


def provenance_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    raw = f"{chrom}|{pos}|{ref}|{alt}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:16]


def parse_ab(ad: str, gt: str) -> float | None:
    if not ad or ad == ".":
        return None
    parts = []
    for token in ad.split(","):
        try:
            parts.append(int(token))
        except ValueError:
            return None
    if not parts:
        return None
    total = sum(parts)
    if total <= 0:
        return None
    alleles = [token for token in re.split(r"[/|]", gt or "") if token not in {"", "."}]
    alt_indexes = []
    for token in alleles:
        try:
            idx = int(token)
        except ValueError:
            continue
        if idx > 0:
            alt_indexes.append(idx)
    if not alt_indexes:
        return None
    alt_depth = sum(parts[i] for i in alt_indexes if i < len(parts))
    return alt_depth / total


def frequency_knowledge(max_af: float | None, gnomadg_af: float | None) -> str:
    if max_af is not None:
        return "frequency_verified"
    if gnomadg_af is not None:
        return "frequency_partially_known"
    return "frequency_unknown"


@dataclass
class VariantCandidate:
    variant_id: str
    gene: str
    transcript: str
    consequence: str
    impact: str
    max_af: float | None
    gnomadg_af: float | None
    frequency_class: str
    clinvar: str
    sift: str
    polyphen: str
    mane: str
    canonical: str
    dp: int | None
    gq: int | None
    ab: float | None
    filt: str
    pid: str
    pgt: str
    existing_variant: str
    existing_status: str
    zygosity: str
    pick: bool = False
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: str = ""


@dataclass
class PairCandidate:
    gene: str
    left: VariantCandidate
    right: VariantCandidate
    pair_frequency_class: str
    phase_status: str
    phase_class: str


def phase_class_from_status(status: str) -> str:
    if status == "same_phase":
        return "likely_cis"
    return "phase_unknown"


def _int_or_none(raw: str) -> int | None:
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return None


def _pick_csq_row(entries: list[dict[str, str]]) -> dict[str, str]:
    def _score(row: dict[str, str]) -> tuple[int, int, int, int]:
        impact = (row.get("IMPACT") or "").upper()
        impact_rank = {"HIGH": 3, "MODERATE": 2, "LOW": 1}.get(impact, 0)
        return (
            1 if row.get("PICK") == "1" else 0,
            1 if row.get("MANE_SELECT") not in {"", "."} else 0,
            1 if row.get("CANONICAL") == "YES" else 0,
            impact_rank,
        )

    return max(entries, key=_score)


def _candidate_from_row(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    filt: str,
    row: dict[str, str],
    gt_info: dict[str, Any] | None,
) -> VariantCandidate | None:
    impact = (row.get("IMPACT") or "").upper()
    if impact not in {"HIGH", "MODERATE"}:
        return None
    gene = (row.get("SYMBOL") or row.get("Gene") or "").strip()
    if not gene or gene == ".":
        return None
    max_af = parse_af(row.get("MAX_AF", ""))
    if max_af is None:
        max_af = parse_af(row.get("AF", ""))
    gnomadg = parse_af(row.get("gnomADg_AF", "") or row.get("gnomAD_g_AF", ""))
    gt = (gt_info or {}).get("gt", "")
    zyg = (gt_info or {}).get("zygosity") or parse_gt_zygosity(gt)
    existing = (row.get("Existing_variation") or "").strip()
    existing_status = "present" if existing not in {"", "."} else "missing"
    return VariantCandidate(
        variant_id=provenance_id(chrom, pos, ref, alt),
        gene=gene,
        transcript=(row.get("Feature") or row.get("MANE_SELECT") or "").strip(),
        consequence=most_severe_consequence([row.get("Consequence", "")]),
        impact=impact,
        max_af=max_af,
        gnomadg_af=gnomadg,
        frequency_class=frequency_class(max_af),
        clinvar=clin_sig_category(row.get("CLIN_SIG", "") or row.get("ClinVar", "")),
        sift=(row.get("SIFT") or "").strip(),
        polyphen=(row.get("PolyPhen") or "").strip(),
        mane=(row.get("MANE_SELECT") or "").strip(),
        canonical=(row.get("CANONICAL") or "").strip(),
        dp=_int_or_none((gt_info or {}).get("dp", "")),
        gq=_int_or_none((gt_info or {}).get("gq", "")),
        ab=parse_ab((gt_info or {}).get("ad", ""), gt),
        filt=filt,
        pid=(gt_info or {}).get("pid", "") or "",
        pgt=(gt_info or {}).get("pgt", "") or "",
        existing_variant=existing,
        existing_status=existing_status,
        zygosity=zyg,
        pick=row.get("PICK") == "1",
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
    )


def load_sample_fields(
    vcf_path: Path,
    zygosities: set[str] | None = None,
) -> dict[tuple[str, int, str, str], dict[str, Any]]:
    """Load GT/DP/GQ/AD/PGT/PID. Keys stay in memory and are not printed."""
    table: dict[tuple[str, int, str, str], dict[str, Any]] = {}
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, _vid, ref, alt_s, _qual, _filt, _info, fmt, sample = parts[:10]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            fmt_map = dict(zip(fmt.split(":"), sample.split(":"), strict=False))
            payload = {
                "gt": fmt_map.get("GT", ""),
                "zygosity": parse_gt_zygosity(fmt_map.get("GT", "")),
                "dp": fmt_map.get("DP", ""),
                "gq": fmt_map.get("GQ", ""),
                "ad": fmt_map.get("AD", ""),
                "pid": fmt_map.get("PID", ""),
                "pgt": fmt_map.get("PGT", ""),
            }
            if zygosities is not None and payload["zygosity"] not in zygosities:
                continue
            for alt in alt_s.split(","):
                table[(chrom, pos, ref, alt)] = payload
    return table


def load_phase_by_locus(vcf_path: Path) -> dict[tuple[str, int], dict[str, str]]:
    table: dict[tuple[str, int], dict[str, str]] = {}
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos_s, _vid, _ref, _alt, _qual, _filt, _info, fmt, sample = parts[:10]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            fmt_map = dict(zip(fmt.split(":"), sample.split(":"), strict=False))
            table[(chrom, pos)] = {
                "pid": fmt_map.get("PID", ""),
                "pgt": fmt_map.get("PGT", ""),
            }
    return table


def stream_high_moderate_candidates(
    vep_vcf: Path,
    genotypes: dict[tuple[str, int, str, str], dict[str, Any]] | None = None,
    phase_by_locus: dict[tuple[str, int], dict[str, str]] | None = None,
) -> list[VariantCandidate]:
    fields: list[str] = []
    out: list[VariantCandidate] = []
    with _open_text(vep_vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                if "ID=CSQ" in line:
                    fields = parse_csq_format(line)
                continue
            if line.startswith("#CHROM"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, _vid, ref, alt, _qual, filt, info_s = parts[:8]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            info = parse_info(info_s)
            csq_blob = info.get("CSQ", "")
            if not csq_blob or not fields:
                continue
            entries = []
            for raw in csq_blob.split(","):
                values = raw.split("|")
                entries.append(
                    {
                        fields[i] if i < len(fields) else f"f{i}": (
                            values[i] if i < len(values) else ""
                        )
                        for i in range(max(len(fields), len(values)))
                    }
                )
            row = _pick_csq_row(entries)
            orig = parse_orig(info.get("ORIG", ""))
            gt_info = _lookup_genotype(genotypes, chrom, pos, ref, alt)
            if orig and genotypes is not None and gt_info is None:
                gt_info = _lookup_genotype(
                    genotypes, orig["chrom"], orig["pos"], orig["ref"], alt
                )
            if phase_by_locus:
                locus = (orig["chrom"], orig["pos"]) if orig else (chrom, pos)
                phase = phase_by_locus.get(locus)
                if phase is None and str(locus[0]).startswith("chr"):
                    phase = phase_by_locus.get((str(locus[0])[3:], locus[1]))
                if phase is None and not str(locus[0]).startswith("chr"):
                    phase = phase_by_locus.get((f"chr{locus[0]}", locus[1]))
                if phase:
                    merged = dict(gt_info or {})
                    if phase.get("pid"):
                        merged["pid"] = phase["pid"]
                    if phase.get("pgt"):
                        merged["pgt"] = phase["pgt"]
                    gt_info = merged
            candidate = _candidate_from_row(chrom, pos, ref, alt, filt, row, gt_info)
            if candidate is not None:
                out.append(candidate)
    return out


def heterozygous_noncommon(candidates: Iterable[VariantCandidate]) -> list[VariantCandidate]:
    return [
        item
        for item in candidates
        if item.zygosity == "het" and item.frequency_class != AF_COMMON
    ]


def frequency_class_counts(candidates: Iterable[VariantCandidate]) -> dict[str, int]:
    counts = {
        AF_OBSERVED_ULTRARARE: 0,
        AF_OBSERVED_RARE: 0,
        AF_ZERO_REPORTED: 0,
        AF_MISSING: 0,
        AF_COMMON: 0,
    }
    for item in candidates:
        counts[item.frequency_class] = counts.get(item.frequency_class, 0) + 1
    return counts


def build_pairs(candidates: Iterable[VariantCandidate]) -> list[PairCandidate]:
    by_gene: dict[str, list[VariantCandidate]] = {}
    for item in heterozygous_noncommon(candidates):
        by_gene.setdefault(item.gene, []).append(item)
    pairs: list[PairCandidate] = []
    for gene, rows in by_gene.items():
        unique: dict[str, VariantCandidate] = {}
        for row in rows:
            unique[row.variant_id] = row
        items = list(unique.values())
        for i, left in enumerate(items):
            for right in items[i + 1 :]:
                status = phase_pair_status(left.pid, right.pid, left.pgt, right.pgt)
                pairs.append(
                    PairCandidate(
                        gene=gene,
                        left=left,
                        right=right,
                        pair_frequency_class=pair_frequency_class(
                            left.frequency_class, right.frequency_class
                        ),
                        phase_status=status,
                        phase_class=phase_class_from_status(status),
                    )
                )
    return pairs


def pair_class_counts(pairs: Iterable[PairCandidate]) -> dict[str, int]:
    counts = {key: 0 for key in PAIR_CLASS_ORDER}
    extra: dict[str, int] = {}
    for pair in pairs:
        key = pair.pair_frequency_class
        if key in counts:
            counts[key] += 1
        else:
            extra[key] = extra.get(key, 0) + 1
    counts.update(extra)
    return counts


def public_variant_view(item: VariantCandidate) -> dict[str, Any]:
    """Safe fields for reports: no coordinates or alleles."""
    return {
        "variant_id": item.variant_id,
        "gene": item.gene,
        "transcript": item.transcript,
        "consequence": item.consequence,
        "impact": item.impact,
        "max_af": item.max_af,
        "gnomadg_af": item.gnomadg_af,
        "frequency_class": item.frequency_class,
        "clinvar": item.clinvar,
        "sift": item.sift,
        "polyphen": item.polyphen,
        "mane": item.mane,
        "canonical": item.canonical,
        "dp": item.dp,
        "gq": item.gq,
        "ab": item.ab,
        "filter": item.filt,
        "has_pid": bool(item.pid) and item.pid != ".",
        "has_pgt": bool(item.pgt) and item.pgt != "." and "|" in item.pgt.replace("/", "|"),
        "existing_status": item.existing_status,
        "frequency_knowledge": frequency_knowledge(item.max_af, item.gnomadg_af),
        "zygosity": item.zygosity,
    }


def public_pair_row(pair: PairCandidate, scores: dict[str, Any] | None = None) -> dict[str, Any]:
    row = {
        "gene": pair.gene,
        "pair_frequency_class": pair.pair_frequency_class,
        "phase_status": pair.phase_status,
        "phase_class": pair.phase_class,
        "left": public_variant_view(pair.left),
        "right": public_variant_view(pair.right),
    }
    if scores:
        row.update(scores)
    return row


def private_pair_record(
    pair: PairCandidate, scores: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Local-only record. May include coordinates. Must not be printed or committed."""
    payload = {
        "gene": pair.gene,
        "pair_frequency_class": pair.pair_frequency_class,
        "phase_status": pair.phase_status,
        "phase_class": pair.phase_class,
        "left": asdict(pair.left),
        "right": asdict(pair.right),
    }
    if scores:
        payload.update(scores)
    return payload
