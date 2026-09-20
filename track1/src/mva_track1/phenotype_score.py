"""Gene-level HPO semantic similarity using disease annotations.

Positive similarity is Resnik+BMA against the best-matching disease of a gene.
Negative/excluded phenotypes are scored separately and never fail a candidate
on their own.
"""

from __future__ import annotations

from dataclasses import dataclass

from mva_track1.hpo_ontology import (
    GENERIC_HPO_IDS,
    HpoOntology,
    best_match_average,
    resnik_normalized,
)

# Grouping terms whose labels start with this prefix are excluded from queries.
GENERIC_NAME_PREFIXES = ("abnormality of ",)


@dataclass
class PhenotypeScores:
    gene: str
    positive_phenotype_score: float
    negative_contradiction_score: float
    best_disease: str
    n_query_terms: int
    n_disease_terms: int
    used_gene_phenotype_fallback: bool


def query_term_set(hpo_ids: list[str], ontology: HpoOntology) -> set[str]:
    out: set[str] = set()
    for hid in hpo_ids:
        primary = ontology.primary_id(hid) or hid
        if primary in GENERIC_HPO_IDS:
            continue
        if ontology.is_obsolete(primary):
            continue
        label = ontology.label(primary).casefold()
        if any(label.startswith(prefix) for prefix in GENERIC_NAME_PREFIXES):
            continue
        if ontology.ic.get(primary, 0.0) <= 0:
            continue
        out.add(primary)
    return out


def _best_disease_score(
    ontology: HpoOntology,
    query: set[str],
    diseases: set[str],
) -> tuple[str, float, int]:
    best_id = ""
    best = 0.0
    best_n = 0
    for disease in diseases:
        terms = ontology.disease_phenotypes.get(disease, set())
        score = best_match_average(ontology, query, terms)
        if score > best:
            best = score
            best_id = disease
            best_n = len(terms)
    return best_id, best, best_n


def _contradiction_score(
    ontology: HpoOntology,
    absent: set[str],
    disease_present: set[str],
    disease_excluded: set[str],
) -> float:
    """Higher means more tension between patient-absent and disease-typical terms.

    Two weak signals, averaged:
    - patient absent vs disease-present typical features (expressivity caution)
    - patient present is handled elsewhere; here we only use absent terms
    Disease-excluded annotations that overlap patient-absent are *concordant*
    and do not add contradiction.
    """
    if not absent or not disease_present:
        return 0.0
    hits = []
    for term in absent:
        if term in disease_excluded:
            continue
        best = max(resnik_normalized(ontology, term, other) for other in disease_present)
        hits.append(best)
    if not hits:
        return 0.0
    return sum(hits) / len(hits)


def score_gene_phenotype(
    gene: str,
    ontology: HpoOntology,
    present_hpo: list[str],
    absent_hpo: list[str] | None = None,
) -> PhenotypeScores:
    query = query_term_set(present_hpo, ontology)
    absent = query_term_set(absent_hpo or [], ontology)
    diseases = ontology.gene_to_diseases.get(gene, set())
    fallback = False
    if diseases:
        best_disease, positive, n_terms = _best_disease_score(ontology, query, diseases)
        disease_present = ontology.disease_phenotypes.get(best_disease, set())
        disease_excl = ontology.disease_excluded.get(best_disease, set())
        negative = _contradiction_score(ontology, absent, disease_present, disease_excl)
    else:
        fallback = True
        gene_pheno = ontology.gene_to_phenotypes.get(gene, set())
        best_disease = ""
        n_terms = len(gene_pheno)
        positive = best_match_average(ontology, query, gene_pheno) if gene_pheno else 0.0
        negative = _contradiction_score(ontology, absent, gene_pheno, set())
    return PhenotypeScores(
        gene=gene,
        positive_phenotype_score=round(positive, 6),
        negative_contradiction_score=round(negative, 6),
        best_disease=best_disease,
        n_query_terms=len(query),
        n_disease_terms=n_terms,
        used_gene_phenotype_fallback=fallback,
    )
