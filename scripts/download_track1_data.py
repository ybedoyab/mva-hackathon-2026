#!/usr/bin/env python3
"""CLI wrapper: prepare a safe download of SageBio/mva-hackathon-2026-data."""

from __future__ import annotations

import sys
from pathlib import Path


def _ensure_src_on_path() -> None:
    root = Path(__file__).resolve().parent.parent
    src = root / "track1" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def main() -> int:
    _ensure_src_on_path()
    from mva_track1.download import main as download_main

    return download_main()


if __name__ == "__main__":
    raise SystemExit(main())
