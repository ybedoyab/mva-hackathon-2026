from __future__ import annotations

from pathlib import Path


def write_mini_ontology(tmp_path: Path) -> Path:
    obo = tmp_path / "hp.obo"
    obo.write_text(
        "\n".join(
            [
                "format-version: 1.2",
                "",
                "[Term]",
                "id: HP:0000001",
                "name: All",
                "",
                "[Term]",
                "id: HP:0000118",
                "name: Phenotypic abnormality",
                "is_a: HP:0000001",
                "",
                "[Term]",
                "id: HP:0000707",
                "name: Abnormality of the nervous system",
                'synonym: "Neurological abnormality" EXACT []',
                "is_a: HP:0000118",
                "",
                "[Term]",
                "id: HP:0001250",
                "name: Seizure",
                'synonym: "Seizures" EXACT []',
                "is_a: HP:0000707",
                "",
                "[Term]",
                "id: HP:0001249",
                "name: Intellectual disability",
                "is_a: HP:0000707",
                "",
                "[Term]",
                "id: HP:0000256",
                "name: Macrocephaly",
                'synonym: "Large head" RELATED []',
                "is_a: HP:0000118",
                "",
                "[Term]",
                "id: HP:0001507",
                "name: Growth abnormality",
                "is_a: HP:0000118",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return obo


def write_mini_hpoa(tmp_path: Path) -> Path:
    header = (
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\t"
        "frequency\tsex\tmodifier\taspect\tbiocuration"
    )
    rows = [
        header,
        "OMIM:1\tD1\t\tHP:0001250\t.\tTAS\t\t\t\t\tP\t.",
        "OMIM:2\tD2\t\tHP:0001249\t.\tTAS\t\t\t\t\tP\t.",
        "OMIM:3\tD3\t\tHP:0001250\t.\tTAS\t\t\t\t\tP\t.",
        "OMIM:3\tD3\t\tHP:0001249\t.\tTAS\t\t\t\t\tP\t.",
        "OMIM:3\tD3\tNOT\tHP:0000256\t.\tTAS\t\t\t\t\tP\t.",
    ]
    for i in range(4, 24):
        rows.append(f"OMIM:{i}\tDx{i}\t\tHP:0001507\t.\tTAS\t\t\t\t\tP\t.")
    path = tmp_path / "phenotype.hpoa"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path
