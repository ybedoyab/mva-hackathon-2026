# Track 1 methodology

This document records implemented Track 1 methods. It does **not** contain patient
HPO terms, gene symbols from the case, variant coordinates, alleles, or clinical
wording. Nothing here is a diagnosis.

Challenge configuration (submission schema, scoring, file inventory) is recorded
in [track1_challenge_spec.md](track1_challenge_spec.md).

## Challenge-informed analysis strategy

The public evaluator describes a **compound-heterozygous** clinically validated
answer key. That is organizer scoring design, not a shortcut to a named gene.

The pipeline therefore ranks **same-gene heterozygous pairs** from HIGH/MODERATE
VEP consequences, while keeping frequency-missing sites distinct from observed
rare sites. A later iteration will compare this blind ranking to disease-specific
knowledge and independent tools. This iteration does **not** use a named-syndrome
gene list.

## 1. Data inventory

Implemented as a filename/size/header-only inspect. Controlled files stay outside
Git.

## 2. Data quality control

Aggregate VCF profiling only (FILTER, DP/GQ/AB distributions, phasing-tag
presence). No read-level FASTQ inspection in this iteration.

## 3. Variant normalization

Local `bcftools norm` (Docker `staphb/bcftools:1.21`) against a contig-matched
GRCh38 FASTA. Multi-allelics are split; original CHROM/POS/REF/ALT are retained
in an `ORIG` INFO tag. REF-mismatch count is recorded. Sites-only VCF is used as
VEP input.

## 4. Variant annotation

Offline Ensembl VEP 116.2 (`ensemblorg/ensembl-vep:release_116.2`) with the
official `homo_sapiens/116_GRCh38` cache, `--network none`. Unrecognized-contig
sites are counted and held in a rescue lane (see §19-equivalent below).

## 5. Population frequency filtering

`MAX_AF` from VEP CSQ is the primary frequency field. `gnomADg_AF` is retained
independently. Each heterozygous HIGH/MODERATE variant is assigned **exactly
one** class:

- `OBSERVED_ULTRARARE`: `MAX_AF > 0` and `<= 1e-4`
- `OBSERVED_RARE`: `MAX_AF > 1e-4` and `<= 1e-3`
- `ZERO_REPORTED`: `MAX_AF == 0`
- `AF_MISSING`: `MAX_AF` absent
- `COMMON`: `MAX_AF > 1e-3`

`AF_MISSING` is **not** treated as rare. Pair classes are:

`observed/observed`, `observed/zero`, `zero/zero`, `observed/missing`,
`zero/missing`, `missing/missing`.

COMMON variants are excluded from compound-het pairing.

## 6. Consequence/severity filtering

Pair construction uses VEP IMPACT `HIGH` or `MODERATE` on a picked CSQ row
(PICK, then MANE_SELECT, then CANONICAL, then impact). Functional scores
distinguish splice acceptor/donor, frameshift, stop-gained, start/stop-lost,
missense, and inframe events. Loss-of-function is **not** assumed pathogenic.

## 7. Inheritance modeling

Single-proband VCF; no parental BAM/VCF. Same-gene heterozygous pairs are
generated. Original `PGT`/`PID` are used only to flag **likely cis** when both
sites share a PID and the alt alleles sit on the same encoded haplotype.
Opposite-phase is reported as a raw status but is **not** claimed as trans.
Missing phase is `phase_unknown`, never trans.

## 8. Phenotype-driven prioritization

### Clinical document

The phenotype DOCX is parsed locally with `python-docx`. Raw wording is never
printed or committed. A private review TSV (source phrases) may exist only under
the local working directory. Safe artifacts contain HPO ids, labels, status, confidence, mapping method,
and scope only.

### HPO release

Official Human Phenotype Ontology **v2026-06-23**
(https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-06-23).
Files used locally (not in Git): `hp.obo`, `phenotype.hpoa`,
`genes_to_phenotype.txt`, `genes_to_disease.txt`, `phenotype_to_genes.txt`.
Published GitHub release SHA-256 digests are verified at download.

### Mapping rules

Each mapped item has `status` ∈ {present, absent, uncertain},
`mapping_method` ∈ {explicit_hpo, exact_label, synonym, lexical_candidate},
and `scope` ∈ {proband, family_history, uncertain_context}.

- Explicit `HP:NNNNNNN` tokens are validated against this release.
- Exact primary labels and OBO synonyms are matched after punctuation/case
  normalization. Leading negation tokens are stripped before label lookup so
  that status is assigned on the same span as the match.
- EXACT/NARROW synonyms are high confidence; BROAD/RELATED are medium.
- Whole-label containment may be proposed only as `lexical_candidate` /
  low confidence and is **excluded** from ranking.
- Obvious negation (`no`, `not`, `without`, `absent`, …) → absent **on that
  span only**. A negated specific finding is not generalized to a broad parent
  such as “Abnormality of the kidney”.
- Grouping-term absences (`Abnormality of …`) are not used for quantitative
  negative scoring unless the document explicitly codes that parent HPO.
- If a more specific present term exists, ancestor matches are dropped.
- Family-history / maternal / paternal / miscarriage context is stored as
  `scope=family_history` and is **excluded** from the proband ranking profile.
- Hedge words (`possible`, `suspected`, …) → uncertain.
- “Not mentioned” is not converted to absent.

Ranking currently uses **high-confidence present proband** terms only.
Medium-confidence present terms are a sensitivity lane (S2). Family-history
terms and low-confidence automatic mappings are never used to rank.

### Semantic similarity

Disease-level annotations come from `phenotype.hpoa` (aspect `P`). Genes are
linked to diseases via `genes_to_disease.txt`. If a gene has no disease link,
`genes_to_phenotype.txt` is a fallback pseudo-profile.

Information content is estimated from the HPOA disease corpus:

```
n(t)  = number of diseases annotated with t or a descendant
N     = number of diseases with ≥1 phenotypic annotation
IC(t) = -ln(n(t)/N)
Resnik(a,b) = IC(MICA(a,b))
Resnik_norm = Resnik / max_t IC(t)
BMA(A,B) = 0.5 * (mean_a max_b Resnik_norm(a,b)
                + mean_b max_a Resnik_norm(a,b))
```

Query terms exclude ontology roots (`HP:0000001`, `HP:0000118`) and labels
starting with “Abnormality of ”. Gene score = max BMA over the gene’s diseases.

### Negative phenotypes

Patient-absent terms produce a separate `negative_contradiction_score` (mean
best Resnik to the matched disease’s typical present terms). A single absent
feature cannot fail a candidate. Primary ranking uses positive evidence;
a small documented penalty (0.05 × negative score) is applied only in the
phenotype-aware lane.

## 9. Pathogenicity evidence

ClinVar (`CLIN_SIG` in the local VEP cache only): pathogenic / likely_pathogenic
support; VUS is not upgraded; conflicting stays conflicting. SIFT/PolyPhen are
supporting missense evidence with a 5% genotype-weight cap. No remote ClinVar,
OMIM, Orphanet, PanelApp, or gnomAD queries.

## 10. Read-level validation

Not implemented in this iteration (FASTQ not aligned).

## 11. Candidate ranking

Two lanes, both disease-agnostic with respect to the named challenge syndrome:

**Lane A (genotype-only).** No HPO. Pair genotype score is the mean of per-allele
component scores:

```
variant = 0.35*frequency + 0.30*functional + 0.15*quality
        + 0.15*clinvar + 0.05*prediction
pair_genotype = mean(variant_left, variant_right) * phase_multiplier
```

Frequency: zero-reported 1.00; observed ≤1e-4 0.90; observed ≤1e-3 0.55;
missing 0.30 (unknown, not ultrarare); common 0.00.
`likely_cis` multiplies the pair by 0.15. `phase_unknown` multiplies by 1.0.

**Lane B (phenotype-aware, still syndrome-agnostic).**

```
overall = 0.55*pair_genotype + 0.45*positive_phenotype
        - 0.05*negative_contradiction
```

clipped to [0, 1]. Gene rank uses the best pair per gene. Ties break on gene
symbol. Reports expose gene symbols and component scores only — not coordinates
or alleles.

## 12. Sensitivity/ablation analyses

- S1: high-confidence present HPO only (primary Lane B)
- S2: high + medium-confidence present HPO
- S3: phenotype score removed (Lane A)
- S4: drop any pair containing `AF_MISSING`
- S5: keep `AF_MISSING` with unknown-frequency score (same policy as S1)
- S6: PASS-only pairs
- S7: PASS plus non-PASS rescue (default pair set; non-PASS already quality-penalized)

Rank stability across S1–S7 is recorded per gene. A robust candidate should not
exist in only one threshold setting.

## Unannotated-contig rescue

Variants whose contigs are absent from the Ensembl GRCh38 cache are **not**
discarded. They are profiled in aggregate (PASS/het/hom-alt, DP/GQ, AF-in-VCF,
decoy/random/alt/unplaced categories) and left unranked by phenotype.

## 13. Final submission generation

Not implemented. No submission file is produced in this iteration.
