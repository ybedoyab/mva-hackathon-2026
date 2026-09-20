#!/usr/bin/env python3
"""CLI wrapper: aggregate VEP annotation statistics without printing variants."""

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
    from mva_track1.annotation_profile import main as profile_main

    return profile_main()


if __name__ == "__main__":
    raise SystemExit(main())
