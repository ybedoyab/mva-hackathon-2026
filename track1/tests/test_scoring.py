from __future__ import annotations

from mva_track1.scoring import score_proband
from mva_track1.submission import SubmissionRow

TRUE_PAIR = frozenset(
    {
        ("chr1", 1000, "A", "G"),
        ("chr1", 2000, "C", "T"),
    }
)
DECOY = ("chr9", 50, "G", "A")
WRONG_SECOND = ("chr3", 9, "T", "A")


def _row(variants: frozenset, epcr: float, rank: int) -> SubmissionRow:
    return SubmissionRow(proband_id="SYNTH01", variants=variants, epcr=epcr, rank=rank)


def test_full_match_rank_one() -> None:
    rows = [_row(TRUE_PAIR, 0.99, 1)]
    result = score_proband("SYNTH01", rows, TRUE_PAIR)
    assert result.full_match_rank == 1
    assert result.partial_match_rank is None
    assert result.rank_points == 100
    assert result.f_max == 1.0


def test_full_match_rank_two() -> None:
    rows = [
        _row(frozenset({DECOY}), 0.95, 1),
        _row(TRUE_PAIR, 0.80, 2),
    ]
    result = score_proband("SYNTH01", rows, TRUE_PAIR)
    assert result.full_match_rank == 2
    assert result.rank_points == 50


def test_one_correct_allele_paired_with_wrong_allele() -> None:
    mixed = frozenset({("chr1", 1000, "A", "G"), WRONG_SECOND})
    rows = [_row(mixed, 1.0, 1)]
    result = score_proband("SYNTH01", rows, TRUE_PAIR)
    assert result.full_match_rank is None
    assert result.partial_match_rank == 1
    assert result.rank_points == 50.0
    assert 0 < result.f_max < 1


def test_no_correct_alleles() -> None:
    rows = [_row(frozenset({DECOY}), 0.7, 1)]
    result = score_proband("SYNTH01", rows, TRUE_PAIR)
    assert result.full_match_rank is None
    assert result.partial_match_rank is None
    assert result.rank_points == 0
    assert result.f_max == 0


def test_fmax_drops_when_extra_false_positives_share_threshold() -> None:
    rows = [
        _row(TRUE_PAIR, 0.9, 1),
        _row(frozenset({DECOY}), 0.9, 2),
    ]
    result = score_proband("SYNTH01", rows, TRUE_PAIR)
    assert result.full_match_rank == 1
    assert result.rank_points == 100
    # Predicted {true, true, decoy} at the only EPCR threshold: precision 2/3, recall 1, F=0.8
    assert result.f_max == 0.8
