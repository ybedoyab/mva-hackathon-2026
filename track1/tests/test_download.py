from __future__ import annotations

from pathlib import Path

import pytest

from mva_track1.download import DownloadError, format_plan, plan_download, run_download
from mva_track1.paths import repo_root


def test_refuses_repository_local_destination() -> None:
    root = repo_root()
    dest = root / "controlled_data"
    plan = plan_download(dest, root=root)
    assert plan.inside_repo is True
    with pytest.raises(DownloadError, match="git repository tree"):
        run_download(plan, confirm=True, dry_run=False)


def test_dry_run_does_not_call_snapshot(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    called = {"snapshot": False}

    def fake_snapshot_download(**kwargs):  # type: ignore[no-untyped-def]
        called["snapshot"] = True
        raise AssertionError("snapshot_download should not run during dry-run")

    monkeypatch.setattr("huggingface_hub.snapshot_download", fake_snapshot_download, raising=False)
    plan = plan_download(tmp_path / "outside", root=repo_root())
    assert plan.inside_repo is False
    run_download(plan, confirm=False, dry_run=True)
    captured = capsys.readouterr()
    assert "DRY RUN" in captured.out
    assert "SageBio/mva-hackathon-2026-data" in captured.out
    assert "must not be committed" in captured.out
    assert called["snapshot"] is False


def test_confirm_required_for_live_download(tmp_path: Path) -> None:
    plan = plan_download(tmp_path / "outside", root=repo_root())
    with pytest.raises(DownloadError, match="confirm-download"):
        run_download(plan, confirm=False, dry_run=False)


def test_plan_does_not_embed_token(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    secret = "hf_" + ("z" * 32)
    monkeypatch.setenv("HF_TOKEN", secret)
    plan = plan_download(tmp_path / "outside", root=repo_root())
    text = format_plan(plan, dry_run=True)
    assert plan.token_present is True
    assert secret not in text
    assert "present" in text
