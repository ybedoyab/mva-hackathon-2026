#!/usr/bin/env python3
"""Extract a local HPO phenotype profile from the clinical DOCX.

Never prints raw clinical text. Safe stdout is structure + HPO ids/labels only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _ensure_src_on_path() -> None:
    root = Path(__file__).resolve().parent.parent
    src = root / "track1" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Local DOCX to HPO; no source text on stdout.")
    parser.add_argument("--docx", type=Path, default=None)
    parser.add_argument("--hpo-dir", type=Path, required=True)
    parser.add_argument("--private-tsv", type=Path, required=True)
    parser.add_argument("--safe-tsv", type=Path, required=True)
    parser.add_argument("--safe-json", type=Path, required=True)
    parser.add_argument("--structure-json", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    _ensure_src_on_path()
    from mva_track1.hpo_ontology import parse_obo
    from mva_track1.phenotype_extract import (
        extract_phenotype,
        find_clinical_docx,
        format_safe_phenotype_report,
        format_structure_report,
        structure_as_dict,
        write_private_candidates,
        write_safe_phenotype,
    )
    from mva_track1.private_outputs import data_root

    args = parse_args(argv)
    docx = args.docx or find_clinical_docx(data_root())
    ontology = parse_obo(args.hpo_dir / "hp.obo")
    structure, candidates = extract_phenotype(docx, ontology)
    write_private_candidates(args.private_tsv, candidates)
    write_safe_phenotype(args.safe_tsv, args.safe_json, candidates)
    print(format_structure_report(structure), end="")
    print(format_safe_phenotype_report(candidates), end="")
    if args.structure_json is not None:
        args.structure_json.parent.mkdir(parents=True, exist_ok=True)
        args.structure_json.write_text(
            json.dumps(structure_as_dict(structure), indent=2) + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
