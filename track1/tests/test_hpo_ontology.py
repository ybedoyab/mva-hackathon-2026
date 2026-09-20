from __future__ import annotations

from pathlib import Path

from hpo_fixtures import write_mini_hpoa, write_mini_ontology
from mva_track1.hpo_ontology import (
    best_match_average,
    load_phenotype_hpoa,
    mica,
    parse_obo,
    resnik,
    resnik_normalized,
)


def test_ancestors_and_specific_over_generic(tmp_path: Path) -> None:
    ontology = parse_obo(write_mini_ontology(tmp_path))
    ancs = ontology.ancestors("HP:0001250")
    assert "HP:0001250" in ancs
    assert "HP:0000707" in ancs
    assert "HP:0000118" in ancs
    assert "HP:0000001" in ancs
    assert "HP:0001250" not in ontology.ancestors("HP:0001250", include_self=False)


def test_information_content_and_resnik_bma(tmp_path: Path) -> None:
    ontology = parse_obo(write_mini_ontology(tmp_path))
    load_phenotype_hpoa(ontology, write_mini_hpoa(tmp_path))
    assert ontology.n_diseases == 23
    ic_sz = ontology.ic["HP:0001250"]
    ic_ns = ontology.ic["HP:0000707"]
    ic_root = ontology.ic["HP:0000118"]
    assert ic_sz > ic_ns > ic_root
    assert mica(ontology, "HP:0001250", "HP:0001249") == "HP:0000707"
    same = resnik(ontology, "HP:0001250", "HP:0001250")
    assert same == ic_sz
    assert 0 < resnik_normalized(ontology, "HP:0001250", "HP:0001249") < 1
    bma_same = best_match_average(ontology, {"HP:0001250"}, {"HP:0001250"})
    bma_diff = best_match_average(ontology, {"HP:0001250"}, {"HP:0001249"})
    assert bma_same > bma_diff
    generic = best_match_average(ontology, {"HP:0000118"}, {"HP:0000118"})
    assert generic == 0.0
