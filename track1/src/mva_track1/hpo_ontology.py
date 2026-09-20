"""Local Human Phenotype Ontology loading and semantic similarity.

Information content is estimated from disease annotations in phenotype.hpoa.
Similarity uses Resnik IC of the most informative common ancestor (MICA) and
Best-Match Average (BMA) over phenotype sets.

Equations (natural log):

    n(t)  = number of diseases annotated with term t or any descendant
    N     = number of diseases with at least one phenotypic (aspect P) annotation
    p(t)  = n(t) / N
    IC(t) = -ln(p(t))   if n(t) > 0 else 0

    Resnik(a, b) = IC(MICA(a, b))
    Resnik_norm(a, b) = Resnik(a, b) / IC_max

    BMA(A, B) = 0.5 * (
        mean_{a in A} max_{b in B} Resnik_norm(a, b)
      + mean_{b in B} max_{a in A} Resnik_norm(a, b)
    )

Generic roots such as HP:0000001 and HP:0000118 are excluded from query sets so
they cannot dominate scores. Low-IC query terms may be dropped by the caller.
"""

from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

HPO_ID_RE = re.compile(r"HP:\d{7}")
GENERIC_HPO_IDS = frozenset(
    {
        "HP:0000001",  # All
        "HP:0000118",  # Phenotypic abnormality
    }
)
SYNONYM_RE = re.compile(
    r'^synonym:\s*"(.+)"\s+(EXACT|BROAD|NARROW|RELATED)\b',
)

NORMALIZE_PUNCT_RE = re.compile(r"[^a-z0-9+\-]+")


def normalize_label(text: str) -> str:
    lowered = (text or "").casefold().replace("'", "")
    return " ".join(NORMALIZE_PUNCT_RE.sub(" ", lowered).split())


@dataclass
class HpoTerm:
    hpo_id: str
    name: str
    synonyms: dict[str, str] = field(default_factory=dict)  # normalized -> scope
    parents: set[str] = field(default_factory=set)
    alt_ids: set[str] = field(default_factory=set)
    obsolete: bool = False
    replaced_by: str = ""


@dataclass
class HpoOntology:
    terms: dict[str, HpoTerm]
    children: dict[str, set[str]]
    label_index: dict[str, list[tuple[str, str]]]
    alt_to_primary: dict[str, str]
    ancestors_cache: dict[str, frozenset[str]] = field(default_factory=dict)
    ic: dict[str, float] = field(default_factory=dict)
    ic_max: float = 0.0
    n_diseases: int = 0
    n_disease_annotations: int = 0
    n_excluded_annotations: int = 0
    disease_phenotypes: dict[str, set[str]] = field(default_factory=dict)
    disease_excluded: dict[str, set[str]] = field(default_factory=dict)
    gene_to_diseases: dict[str, set[str]] = field(default_factory=dict)
    gene_to_phenotypes: dict[str, set[str]] = field(default_factory=dict)

    def primary_id(self, hpo_id: str) -> str | None:
        if hpo_id in self.terms:
            return hpo_id
        return self.alt_to_primary.get(hpo_id)

    def label(self, hpo_id: str) -> str:
        term = self.terms.get(hpo_id)
        return term.name if term else ""

    def is_obsolete(self, hpo_id: str) -> bool:
        term = self.terms.get(hpo_id)
        return bool(term and term.obsolete)

    def ancestors(self, hpo_id: str, include_self: bool = True) -> frozenset[str]:
        if hpo_id in self.ancestors_cache:
            cached = self.ancestors_cache[hpo_id]
        else:
            seen: set[str] = set()
            stack = [hpo_id]
            while stack:
                current = stack.pop()
                term = self.terms.get(current)
                if term is None:
                    continue
                for parent in term.parents:
                    if parent not in seen:
                        seen.add(parent)
                        stack.append(parent)
            cached = frozenset(seen)
            self.ancestors_cache[hpo_id] = cached
        if include_self:
            return cached | {hpo_id}
        return cached

    def descendants_include(self, hpo_id: str) -> set[str]:
        out = {hpo_id}
        stack = [hpo_id]
        while stack:
            current = stack.pop()
            for child in self.children.get(current, ()):
                if child not in out:
                    out.add(child)
                    stack.append(child)
        return out


def parse_obo(path: Path) -> HpoOntology:
    terms: dict[str, HpoTerm] = {}
    children: dict[str, set[str]] = defaultdict(set)
    alt_to_primary: dict[str, str] = {}
    current: dict[str, str | list[str] | set[str] | bool] | None = None

    def _flush() -> None:
        nonlocal current
        if not current:
            return
        hpo_id = str(current.get("id", ""))
        if not hpo_id.startswith("HP:"):
            current = None
            return
        term = HpoTerm(
            hpo_id=hpo_id,
            name=str(current.get("name", "")),
            synonyms=dict(current.get("synonyms", {})),  # type: ignore[arg-type]
            parents=set(current.get("parents", set())),  # type: ignore[arg-type]
            alt_ids=set(current.get("alt_ids", set())),  # type: ignore[arg-type]
            obsolete=bool(current.get("obsolete", False)),
            replaced_by=str(current.get("replaced_by", "")),
        )
        terms[hpo_id] = term
        for parent in term.parents:
            children[parent].add(hpo_id)
        for alt in term.alt_ids:
            alt_to_primary[alt] = hpo_id
        current = None

    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                _flush()
                current = {
                    "synonyms": {},
                    "parents": set(),
                    "alt_ids": set(),
                    "obsolete": False,
                    "replaced_by": "",
                }
                continue
            if current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:].strip()
            elif line.startswith("name: "):
                current["name"] = line[6:].strip()
            elif line.startswith("is_a: "):
                parent = line[6:].split("!", 1)[0].strip()
                if parent.startswith("HP:"):
                    current["parents"].add(parent)  # type: ignore[union-attr]
            elif line.startswith("alt_id: "):
                current["alt_ids"].add(line[8:].strip())  # type: ignore[union-attr]
            elif line.startswith("is_obsolete:"):
                current["obsolete"] = "true" in line.casefold()
            elif line.startswith("replaced_by: "):
                current["replaced_by"] = line[13:].strip()
            else:
                match = SYNONYM_RE.match(line)
                if match:
                    syn_map: dict[str, str] = current["synonyms"]  # type: ignore[assignment]
                    syn_map[normalize_label(match.group(1))] = match.group(2)

    _flush()

    label_index: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for hpo_id, term in terms.items():
        if term.obsolete:
            continue
        if term.name:
            label_index[normalize_label(term.name)].append((hpo_id, "EXACT_LABEL"))
        for syn, scope in term.synonyms.items():
            label_index[syn].append((hpo_id, scope))

    return HpoOntology(
        terms=terms,
        children=dict(children),
        label_index=dict(label_index),
        alt_to_primary=alt_to_primary,
    )


def load_phenotype_hpoa(ontology: HpoOntology, path: Path) -> None:
    """Fill disease phenotype sets and information content from phenotype.hpoa."""
    disease_pheno: dict[str, set[str]] = defaultdict(set)
    disease_excl: dict[str, set[str]] = defaultdict(set)
    header: list[str] | None = None
    n_ann = 0
    n_not = 0
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts, strict=False))
            aspect = (row.get("aspect") or "").strip()
            if aspect != "P":
                continue
            disease = (row.get("database_id") or "").strip()
            hpo_id = ontology.primary_id((row.get("hpo_id") or "").strip())
            if not disease or not hpo_id:
                continue
            qualifier = (row.get("qualifier") or "").strip().upper()
            if qualifier == "NOT":
                disease_excl[disease].add(hpo_id)
                n_not += 1
                continue
            disease_pheno[disease].add(hpo_id)
            n_ann += 1

    term_diseases: dict[str, set[str]] = defaultdict(set)
    for disease, terms in disease_pheno.items():
        for term in terms:
            for anc in ontology.ancestors(term, include_self=True):
                term_diseases[anc].add(disease)

    n_diseases = len(disease_pheno)
    ic: dict[str, float] = {}
    ic_max = 0.0
    if n_diseases:
        for hpo_id, diseases in term_diseases.items():
            p = len(diseases) / n_diseases
            value = -math.log(p) if p > 0 else 0.0
            ic[hpo_id] = value
            if value > ic_max:
                ic_max = value

    ontology.ic = ic
    ontology.ic_max = ic_max
    ontology.n_diseases = n_diseases
    ontology.n_disease_annotations = n_ann
    ontology.n_excluded_annotations = n_not
    ontology.disease_phenotypes = dict(disease_pheno)
    ontology.disease_excluded = dict(disease_excl)


def load_gene_disease_map(ontology: HpoOntology, path: Path) -> None:
    mapping: dict[str, set[str]] = defaultdict(set)
    header: list[str] | None = None
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [part.lstrip("#") for part in parts]
                continue
            row = dict(zip(header, parts, strict=False))
            gene = (row.get("gene_symbol") or "").strip()
            disease = (row.get("disease_id") or "").strip()
            if gene and disease:
                mapping[gene].add(disease)
    ontology.gene_to_diseases = dict(mapping)


def load_gene_phenotype_map(ontology: HpoOntology, path: Path) -> None:
    mapping: dict[str, set[str]] = defaultdict(set)
    header: list[str] | None = None
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [part.lstrip("#") for part in parts]
                continue
            row = dict(zip(header, parts, strict=False))
            gene = (row.get("gene_symbol") or "").strip()
            hpo_id = ontology.primary_id((row.get("hpo_id") or "").strip())
            if gene and hpo_id:
                mapping[gene].add(hpo_id)
    ontology.gene_to_phenotypes = dict(mapping)


def mica(ontology: HpoOntology, term_a: str, term_b: str) -> str | None:
    shared = ontology.ancestors(term_a, include_self=True) & ontology.ancestors(
        term_b, include_self=True
    )
    if not shared:
        return None
    return max(shared, key=lambda hid: ontology.ic.get(hid, 0.0))


def resnik(ontology: HpoOntology, term_a: str, term_b: str) -> float:
    common = mica(ontology, term_a, term_b)
    if common is None:
        return 0.0
    return ontology.ic.get(common, 0.0)


def resnik_normalized(ontology: HpoOntology, term_a: str, term_b: str) -> float:
    if ontology.ic_max <= 0:
        return 0.0
    return resnik(ontology, term_a, term_b) / ontology.ic_max


def best_match_average(
    ontology: HpoOntology,
    set_a: set[str],
    set_b: set[str],
) -> float:
    query = {hid for hid in set_a if hid not in GENERIC_HPO_IDS}
    target = {hid for hid in set_b if hid not in GENERIC_HPO_IDS}
    if not query or not target:
        return 0.0

    def _directed(src: set[str], dst: set[str]) -> float:
        scores = [
            max(resnik_normalized(ontology, a, b) for b in dst)
            for a in src
        ]
        return sum(scores) / len(scores)

    return 0.5 * (_directed(query, target) + _directed(target, query))


def lookup_normalized(ontology: HpoOntology, text: str) -> list[tuple[str, str]]:
    key = normalize_label(text)
    if not key:
        return []
    return list(ontology.label_index.get(key, []))
