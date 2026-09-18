from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from shutil import which

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CHECKER_PATH = REPO_ROOT / "scripts" / "check_repo_safety.py"
HAS_GIT = which("git") is not None


def load_checker():
    spec = importlib.util.spec_from_file_location("check_repo_safety", CHECKER_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


checker = load_checker()


def test_genomic_extension_is_flagged(tmp_path: Path) -> None:
    vcf = tmp_path / "proband.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    issues = checker.inspect_file(vcf, tmp_path)
    assert any(issue.reason == "genomic data file extension" for issue in issues)


def test_env_example_is_allowed_but_env_is_not(tmp_path: Path) -> None:
    example = tmp_path / ".env.example"
    example.write_text("HF_TOKEN=\n", encoding="utf-8")
    env = tmp_path / ".env"
    env.write_text("HF_TOKEN=\n", encoding="utf-8")

    example_issues = checker.inspect_file(example, tmp_path)
    env_issues = checker.inspect_file(env, tmp_path)

    assert not any("credential" in issue.reason for issue in example_issues)
    assert any("credential" in issue.reason for issue in env_issues)


def test_large_file_requires_allowlist(tmp_path: Path) -> None:
    blob = tmp_path / "big.bin"
    blob.write_bytes(b"0" * (checker.MAX_FILE_BYTES + 1))

    flagged = checker.inspect_file(blob, tmp_path)
    assert any("exceeds 25 MB" in issue.reason for issue in flagged)

    allowed = checker.inspect_file(
        blob,
        tmp_path,
        large_file_allowlist=frozenset({"big.bin"}),
    )
    assert not any("exceeds 25 MB" in issue.reason for issue in allowed)


def test_huggingface_token_is_flagged_without_being_printed(tmp_path: Path) -> None:
    token = "hf_" + ("a" * 32)
    notes = tmp_path / "notes.txt"
    notes.write_text(f"export TOKEN={token}\n", encoding="utf-8")

    issues = checker.inspect_file(notes, tmp_path)
    report = checker.format_report(issues, scanned=1)

    assert any("Hugging Face" in issue.reason for issue in issues)
    assert token not in report
    assert "a" * 32 not in report


def test_github_token_is_flagged_without_being_printed(tmp_path: Path) -> None:
    token = "ghp_" + ("b" * 36)
    notes = tmp_path / "notes.txt"
    notes.write_text(f"export TOKEN={token}\n", encoding="utf-8")

    issues = checker.inspect_file(notes, tmp_path)
    report = checker.format_report(issues, scanned=1)

    assert any("GitHub" in issue.reason for issue in issues)
    assert token not in report
    assert "b" * 36 not in report


def test_clean_text_file_has_no_issues(tmp_path: Path) -> None:
    readme = tmp_path / "README.md"
    readme.write_text("# example\n", encoding="utf-8")
    assert checker.inspect_file(readme, tmp_path) == []


def _git_init(directory: Path) -> None:
    subprocess.run(["git", "init"], cwd=directory, check=True, capture_output=True)


@pytest.mark.skipif(not HAS_GIT, reason="git is required for candidate-file tests")
def test_gitignored_vcf_is_not_a_candidate(tmp_path: Path) -> None:
    _git_init(tmp_path)
    (tmp_path / ".gitignore").write_text("*.vcf\n", encoding="utf-8")
    (tmp_path / "hidden.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    (tmp_path / "ok.txt").write_text("hello\n", encoding="utf-8")

    candidates = {path.name for path in checker.iter_candidate_files(tmp_path)}
    assert "ok.txt" in candidates
    assert "hidden.vcf" not in candidates


@pytest.mark.skipif(not HAS_GIT, reason="git is required for candidate-file tests")
def test_main_fails_on_untracked_vcf(tmp_path: Path) -> None:
    _git_init(tmp_path)
    (tmp_path / "sample.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    assert checker.main(["--root", str(tmp_path)]) == 1


def test_checker_script_module_entry(tmp_path: Path) -> None:
    clean = tmp_path / "ok.txt"
    clean.write_text("hello\n", encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(CHECKER_PATH), "--root", str(tmp_path)],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "OK" in result.stdout
