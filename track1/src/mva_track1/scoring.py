"""Public Track 1 scoring logic (synthetic use only).

Mirrors SageBio evaluation.py. This module never loads the challenge answer key
and must only be tested with synthetic variants.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva_track1.submission import SubmissionRow, Variant

# Official rank-point scale from SageBio evaluation.py.
RANK_POINT_TIERS = [
    (1, 100),
    (3, 50),
    (5, 25),
    (10, 10),
]


@dataclass(frozen=True)
class ScoreResult:
    proband_id: str
    full_match_rank: int | None
    partial_match_rank: int | None
    rank_points: float
    f_max: float
    f_max_threshold: float | None
    n_predictions_at_f_max: int


def rank_to_points(rank: int) -> int:
    for max_rank, points in RANK_POINT_TIERS:
        if rank <= max_rank:
            return points
    return 0


def score_proband(
    proband_id: str,
    submission_rows: list[SubmissionRow],
    true_variants: frozenset[Variant],
) -> ScoreResult:
    """Score a submission against a caller-supplied true variant set.

    `true_variants` must be synthetic in tests. Do not pass patient variants.
    """
    is_compound_het = len(true_variants) == 2

    full_match_rank = None
    partial_match_rank = None

    for row in submission_rows:
        if row.variants == true_variants:
            full_match_rank = row.rank
            break

    if is_compound_het and full_match_rank is None:
        for row in submission_rows:
            if row.variants & true_variants:
                partial_match_rank = row.rank
                break

    if full_match_rank is not None:
        rank_points = float(rank_to_points(full_match_rank))
    elif partial_match_rank is not None:
        rank_points = 0.5 * rank_to_points(partial_match_rank)
    else:
        rank_points = 0.0

    thresholds = sorted({row.epcr for row in submission_rows}, reverse=True)
    best_f = 0.0
    best_threshold: float | None = None
    best_n = 0

    for threshold in thresholds:
        predicted_variants: set[Variant] = set()
        n_rows_at_threshold = 0
        for row in submission_rows:
            if row.epcr >= threshold:
                predicted_variants |= set(row.variants)
                n_rows_at_threshold += 1

        true_pos = len(predicted_variants & true_variants)
        false_pos = len(predicted_variants - true_variants)
        false_neg = len(true_variants - predicted_variants)

        precision = 0.0
        if true_pos + false_pos > 0:
            precision = true_pos / (true_pos + false_pos)
        recall = 0.0
        if true_pos + false_neg > 0:
            recall = true_pos / (true_pos + false_neg)
        if precision + recall > 0:
            f_measure = 2 * precision * recall / (precision + recall)
        else:
            f_measure = 0.0

        if f_measure > best_f:
            best_f = f_measure
            best_threshold = threshold
            best_n = n_rows_at_threshold

    return ScoreResult(
        proband_id=proband_id,
        full_match_rank=full_match_rank,
        partial_match_rank=partial_match_rank,
        rank_points=rank_points,
        f_max=best_f,
        f_max_threshold=best_threshold,
        n_predictions_at_f_max=best_n,
    )
