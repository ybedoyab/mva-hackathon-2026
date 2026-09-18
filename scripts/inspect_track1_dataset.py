#!/usr/bin/env python3
"""CLI wrapper: inspect a local Track 1 dataset without reading patient contents."""

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
    from mva_track1.dataset_inspect import main as inspect_main

    return inspect_main()


if __name__ == "__main__":
    raise SystemExit(main())
