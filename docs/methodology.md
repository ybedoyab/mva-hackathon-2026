# Track 1 methodology (skeleton)

This document is a **planning skeleton** for Track 1. None of the scientific analysis stages below are implemented. No results are claimed.

Tool choices, thresholds, and annotation sources will be decided after a local data inventory. Do not treat this outline as a finished pipeline or as evidence for any variant.

Challenge configuration (submission schema, scoring, file inventory) is recorded in [track1_challenge_spec.md](track1_challenge_spec.md). That specification work is separate from variant prioritization.

## Challenge-informed analysis strategy

The public evaluator describes a **compound-heterozygous** clinically validated answer key. That is organizer scoring design, not a shortcut to a named gene.

Therefore, when scientific stages are later implemented, the pipeline should **prioritize robust identification and ranking of plausible compound-heterozygous pairs** (same-gene rare deleterious combinations, cis/trans evidence when available, and pair-level EPCR), while still retaining **broader candidate-generation** (including single-variant hypotheses) for scientific defensibility.

We do **not** want a brittle solution that merely searches one known MVA gene.

Future scientific analysis should independently evaluate:

- rarity;
- consequence;
- gene-level recessive pairing;
- pathogenicity evidence;
- phenotype consistency;
- phasing / cis-trans evidence where possible;
- read support;
- alternative candidate explanations.

The current public file listing is single-proband (VCF + FASTQ lanes + clinical DOCX), with no parental sequencing files. Recessive pairing must therefore be argued from the proband VCF, annotation, and (if produced locally) alignments — not from a trio. BAM/CRAM files are not provided; read support would require local alignment of the FASTQs, with genome-scale intermediates kept outside git and later deleted per challenge rules.

No candidate gene names or variants are asserted here.

## 1. Data inventory

Status: not implemented

Record which files were provided, their types, genome build (once known), and which individuals they belong to. No analysis beyond listing and documenting local inputs.

## 2. Data quality control

Status: not implemented

Inspect sequencing and variant-call quality so later filters are not applied blindly. Metrics and exclusion rules are not defined yet.

## 3. Variant normalization

Status: not implemented

Put variants into a consistent representation (build, allele trimming/left-alignment, and multi-allelic handling) before annotation and comparison.

## 4. Variant annotation

Status: not implemented

Attach consequence, gene, and other reference-database fields. Annotation sources and versions will be recorded when chosen; none are selected here.

## 5. Population frequency filtering

Status: not implemented

Down-weight or exclude variants that are too common in population catalogs relative to a rare-disease prior. Thresholds are unset.

## 6. Consequence/severity filtering

Status: not implemented

Use predicted functional impact as one filter among others. The impact scale and exceptions (for example, noncoding hypotheses) are not defined yet.

## 7. Inheritance modeling

Status: not implemented

The public challenge inventory is single-proband; parental VCF/BAM files are not provided. Recessive / compound-het pairing must be modeled without a sequenced trio unless authorized family data appear later. Do not assume de novo evidence from parents that we do not have.

## 8. Phenotype-driven prioritization

Status: not implemented

Match gene/variant hypotheses to the recorded phenotype. Phenotype vocabularies and similarity methods are not chosen yet.

## 9. Pathogenicity evidence

Status: not implemented

Gather independent evidence (in silico scores, conservation, published criteria, and similar) without treating any single score as decisive.

## 10. Read-level validation

Status: not implemented

Inspect alignments for candidate variants to flag likely artifacts. This step requires local BAM/CRAM access and will not write reads into the repository.

## 11. Candidate ranking

Status: not implemented

Produce an ordered shortlist for submission. The ranking function, features, and score combination are not implemented.

## 12. Sensitivity/ablation analyses

Status: not implemented

Test whether ranking is brittle to individual filters or annotation sources. No ablation has been run.

## 13. Final submission generation

Status: not implemented

Export the organizer-required submission format from a recorded configuration and software environment. No submission file is produced by this repository yet.
