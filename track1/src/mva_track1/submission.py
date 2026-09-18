"""Track 1 submission CSV schema (public official format).

Validates submissions without using ground truth.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COLUMNS = (
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
)
OPTIONAL_COLUMNS = ("finding_type", "notes")
MAX_ROWS = 10
ALLOWED_FINDING_TYPES = frozenset({"primary", "secondary"})
TEMPLATE_COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS

Variant = tuple[str, int, str, str]


class SubmissionError(ValueError):
    """Invalid Track 1 submission CSV."""


@dataclass(frozen=True)
class SubmissionRow:
    """One ranked prediction row after EPCR sort."""

    proband_id: str
    variants: frozenset[Variant]
    epcr: float
    rank: int
    finding_type: str = "primary"
    notes: str = ""


def empty_template_csv() -> str:
    """Return a header-only official submission template. No predictions."""
    return ",".join(TEMPLATE_COLUMNS) + "\n"


def write_empty_template(path: Path) -> None:
    path.write_text(empty_template_csv(), encoding="utf-8")


def _cell(row: dict[str, str], name: str) -> str:
    value = row.get(name, "")
    if value is None:
        return ""
    return str(value).strip()


def _is_blank(value: str) -> bool:
    return value == ""


def _parse_position(raw: str, field: str) -> int:
    if _is_blank(raw):
        raise SubmissionError(f"{field} must be an integer position")
    try:
        pos = int(raw)
    except ValueError as exc:
        raise SubmissionError(f"{field} must be an integer position, got {raw!r}") from exc
    if pos <= 0:
        raise SubmissionError(f"{field} must be a positive integer, got {pos}")
    return pos


def _parse_epcr(raw: str) -> float:
    try:
        epcr = float(raw)
    except ValueError as exc:
        raise SubmissionError(f"epcr must be a float in (0, 1], got {raw!r}") from exc
    if not (0 < epcr <= 1):
        raise SubmissionError(f"epcr {epcr} out of range (0, 1]")
    return epcr


def _parse_variant(chrom: str, pos: str, ref: str, alt: str, *, label: str) -> Variant:
    missing = [
        name
        for name, value in (
            (f"chrom_{label}", chrom),
            (f"pos_{label}", pos),
            (f"ref_{label}", ref),
            (f"alt_{label}", alt),
        )
        if _is_blank(value)
    ]
    if missing:
        raise SubmissionError(f"variant {label} is incomplete; missing {', '.join(missing)}")
    return (chrom, _parse_position(pos, f"pos_{label}"), ref.upper(), alt.upper())


def _second_variant_started(chrom: str, pos: str, ref: str, alt: str) -> bool:
    return any(not _is_blank(value) for value in (chrom, pos, ref, alt))


def parse_raw_rows(rows: Sequence[dict[str, str]]) -> tuple[list[SubmissionRow], bool]:
    """Validate raw CSV dicts and return ranked rows plus whether EPCR was pre-sorted."""
    if not rows:
        raise SubmissionError("submission has no data rows")
    if len(rows) > MAX_ROWS:
        raise SubmissionError(f"submission has {len(rows)} rows; max is {MAX_ROWS}")

    parsed: list[tuple[int, str, frozenset[Variant], float, str, str]] = []
    epcrs: list[float] = []
    for index, row in enumerate(rows):
        pid = _cell(row, "proband_id")
        if _is_blank(pid):
            raise SubmissionError(f"row {index + 1}: proband_id is required")

        v1 = _parse_variant(
            _cell(row, "chrom_1"),
            _cell(row, "pos_1"),
            _cell(row, "ref_1"),
            _cell(row, "alt_1"),
            label="1",
        )
        chrom2, pos2, ref2, alt2 = (
            _cell(row, "chrom_2"),
            _cell(row, "pos_2"),
            _cell(row, "ref_2"),
            _cell(row, "alt_2"),
        )
        if _second_variant_started(chrom2, pos2, ref2, alt2):
            v2 = _parse_variant(chrom2, pos2, ref2, alt2, label="2")
            if v1 == v2:
                raise SubmissionError(f"row {index + 1}: second variant duplicates the first")
            variants = frozenset({v1, v2})
        else:
            variants = frozenset({v1})

        epcr = _parse_epcr(_cell(row, "epcr"))
        finding_type = (_cell(row, "finding_type") or "primary").lower()
        if finding_type not in ALLOWED_FINDING_TYPES:
            raise SubmissionError(
                f"row {index + 1}: finding_type must be "
                f"'primary' or 'secondary', got {finding_type!r}"
            )
        notes = _cell(row, "notes")
        parsed.append((index, pid, variants, epcr, finding_type, notes))
        epcrs.append(epcr)

    was_sorted = epcrs == sorted(epcrs, reverse=True)
    ranked_source = sorted(parsed, key=lambda item: (-item[3], item[0]))
    ranked: list[SubmissionRow] = []
    for rank, item in enumerate(ranked_source, start=1):
        _orig, pid, variants, epcr, finding_type, notes = item
        ranked.append(
            SubmissionRow(
                proband_id=pid,
                variants=variants,
                epcr=epcr,
                rank=rank,
                finding_type=finding_type,
                notes=notes,
            )
        )
    return ranked, was_sorted


def load_submission(csv_path: Path) -> tuple[list[SubmissionRow], bool]:
    """Load and validate a Track 1 CSV. Does not use ground truth."""
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SubmissionError("CSV is missing a header row")
        fields = [name.strip() for name in reader.fieldnames if name is not None]
        missing = [name for name in REQUIRED_COLUMNS if name not in fields]
        if missing:
            raise SubmissionError(f"CSV missing required columns: {', '.join(missing)}")
        rows = [{(k or "").strip(): (v or "") for k, v in row.items()} for row in reader]
    return parse_raw_rows(rows)
