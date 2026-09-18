# Track 1 challenge specification

This document records **public** Track 1 facts from SageBio / Hugging Face sources, then our engineering interpretation, then unknowns.

It is not a diagnosis. It does not name a candidate gene or causal variant. It does not replace the organizers' rules.

Sources used: official challenge Space code and copy, the gated dataset landing page and public file tree, organizer FAQ/rules, and the Sage Bionetworks privacy clarification. Competitor repositories and leaderboard solutions were not consulted.

---

## Official facts

Facts below are stated by the organizers or their public code.

### Challenge and task

- Track 1 is a single-proband ("N-of-1") causal variant prioritization task. The public evaluator is "simplified for a single proband (`N-of-1`)" ([`evaluation.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)).
- The clinically validated answer key is **compound heterozygous**. That is stated in the evaluator module docstring and FAQ scoring text, as challenge-design information, not as a gene-level hypothesis ([`evaluation.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py); [FAQ tab source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py)).
- Participants receive the child's genomic data (FASTQ + VCF) and a clinical/symptom description, and submit ranked variant predictions plus a GitHub link and methods write-up ([Overview tab source](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py)).
- Automated scoring is not the only Track 1 judgment. The methods write-up is reviewed by a judging panel for scientific rigor ([Overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py); [FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py)).

### Dataset repository and size

- Dataset repository: [`SageBio/mva-hackathon-2026-data`](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data).
- Dataset size: approximately **85 GB** ([dataset card](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)).
- Organizers recommend **100–150 GB** of free storage for caching, indexing, and intermediates ([dataset card](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)).
- Access is gated; redistribution is prohibited; the underlying data are not an open licence ([dataset card](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data); [Official Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py); [Sage blog](https://blog.synapse.org/sagebio-blog-blog/2026/08/rare-disease-real-kid-why-we-re-running-a-hackathon-around-one-child-s-genome)).

### Public file inventory

The public Hugging Face dataset tree currently lists these objects ([dataset tree API](https://huggingface.co/api/datasets/SageBio/mva-hackathon-2026-data/tree/main)):

- one clinical phenotype Microsoft Word document (`Challenge_Clinical_Phenotype_1.docx`);
- one compressed VCF (`*.vcf.gz`);
- one VCF tabix index (`*.vcf.gz.tbi`);
- eight paired-end FASTQ files (`*_L001` through `*_L004`, each with `_R1_` and `_R2_`);
- repository metadata (`README.md`, `.gitattributes`).

No parental VCF, BAM, or CRAM files appear in that public listing.

The dataset card's download example uses `snapshot_download` / `hf download` with `ignore_patterns` / `--exclude` for `README.md` and `.gitattributes` ([dataset card](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)).

### Reference build

- Submission coordinates **must use GRCh38** ([Submit Track 1 tab](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)).
- Official field documentation uses chromosome labels with a `chr` prefix in examples (for example `chr15`) ([Submit Track 1 tab](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)). That is a format example, not a statement of the causal locus.

### Submission schema

Official CSV columns from the evaluator ([`evaluation.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)):

```csv
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr
```

The current evaluator also reads optional `finding_type` (`primary` or `secondary`; default `primary`). The public submission template and submit tab additionally document optional `notes` ([template CSV](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/static/templates/track1_submission_template.csv); [Submit Track 1 tab](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)).

Other official parser rules ([`evaluation.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)):

- one row per proposed causal variant or compound-het pair;
- second-variant fields blank for single-variant proposals;
- the first variant is required;
- up to **10** rows per proband;
- `epcr` in `(0, 1]`;
- rows are sorted by EPCR descending before scoring (ties keep original order);
- the submitter need not pre-sort; the official script sorts and is documented to warn if the file was not already ranked;
- REF/ALT are stripped and uppercased;
- the public submit handler only accepts proband id `PROBAND01` ([Submit Track 1 tab](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)).

### Scoring

Rank-point tiers from the official evaluator ([`evaluation.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py)):

```text
rank 1       -> 100
rank 2-3     -> 50
rank 4-5     -> 25
rank 6-10    -> 10
```

Additional official scoring behavior:

- a **full match** is a submitted row whose variant set equals the true set;
- if the true set has two variants (compound het) and there is no full match, **partial rank credit** is given when a row intersects the true set (one true variant alone, or one true variant paired with an incorrect second variant): half of the rank-tier points;
- **F-max** is computed at the **individual-variant** level, sweeping unique EPCR thresholds in the submission;
- secondary/incidental rows are not a separate automated metric; the FAQ says they will not hurt the automated score and are reviewed qualitatively ([FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py)).

### Quota and deliverables

- Current Track 1 submission quota: **6** per Hugging Face user ([`config.py`](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/config.py); [Overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py); [FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py)).
- Only the highest-scoring Track 1 submission appears on the leaderboard ([FAQ](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/faq.py)).
- A methods report (PDF or Markdown) and GitHub URL are required at submit time ([Submit Track 1 tab](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/submit_track1.py)).

### Data handling (organizer clarification)

Sage's Chief Privacy and Compliance Officer stated in the official Community thread ([discussion #2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)):

- material that **carries the child's genome** must be deleted according to challenge rules, including genome-scale derivatives and intermediates;
- durable research outputs may include ranked candidate variants, HPO terms, gene/pathway rankings, mechanism analyses, code, and reports, if prohibited genome-scale data are not embedded;
- third-party hosted services may be used only when they act as a processor (no training on inputs/outputs, no provider reuse rights, retention limited in time and purpose);
- provider, plan/tier, and data-handling settings should be recorded in the methods.

Hackathon Rules also require deletion of data from environments under participant control within 30 days of close, with email confirmation ([Official Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py)).

---

## Design implications

Engineering interpretation of the facts above. These are not additional organizer claims.

- Configure the project as **single-proband WGS-like inputs**: one VCF + index, eight FASTQs from four lanes, one clinical DOCX. Do not assume a sequenced trio.
- BAM/CRAM files are **not** in the current public listing. Read-level review, if done, requires local alignment from FASTQ (outputs kept outside git).
- Rank **compound-heterozygous pairs** as first-class submission rows because that is how a full match is defined. Keep broader candidate generation (including single-variant rows) so the method is scientifically defensible and not a one-gene shortcut.
- Use at most 10 rows. Spend them on high-quality primary pairs; treat `finding_type=secondary` as optional extras that must not crowd out the primary ranking.
- Sort submissions by EPCR descending before upload even though the evaluator will sort.
- Use GRCh38 coordinates and be prepared to match the evaluator's chromosome string exactly (likely `chr`-prefixed, pending local VCF header inspection).
- Store controlled data **outside** this git tree. Download tooling must refuse repo-local destinations.
- Do not send VCF body, FASTQ reads, or clinical DOCX text to external APIs or assistants unless a later, verified-safe workflow is documented.

---

## Unknowns

Not established from public organizer sources, or not inspected yet because this iteration does not open controlled files.

- Exact reference FASTA filename, source URL, and MD5/SHA-256.
- Sequencing platform / chemistry (FASTQ names look Illumina-like; that is not an official platform statement).
- Exact upstream variant-calling pipeline and VCF producer (GATK, DeepVariant, etc.).
- Whether VCF records are normalized (left-aligned, trimmed, split multi-allelics).
- Phasing quality, and whether usable phase (`|` genotypes or other phase info) exists. This cannot be known from header metadata alone without reading variant records.
- Clinical phenotype / HPO content (the DOCX must not be opened in this iteration).
- Read depth, coverage uniformity, contamination, and other QC characteristics.
- Annotation status of the provided VCF (already annotated vs sites-only).
- Contig naming in the actual VCF header (`chr1` vs `1`), sample count, INFO/FORMAT IDs.
- Whether any restricted README text on the gated dataset adds file-level documentation beyond the public card (the dataset README could not be fetched here without gated access).
- How `finding_type=secondary` is excluded from automated scoring internally, if at all (FAQ says it will not hurt the automated score; `evaluation.py` still includes every loaded row in rank and F-max calculations).
- Precise meaning of "optionally BAM/CRAM" and "HPO terms" in the Rules tab versus the current public file tree (see discrepancies below).

### Noted source discrepancies

These are differences among official pages/code, recorded so we do not silently pick one.

1. **File count wording:** the dataset card currently says "~85 GB across **11** files" in the full-download section, while some summaries say **10** files. The public tree has 13 paths including metadata, or 11 if `README.md` and `.gitattributes` are excluded.
2. **Rules vs listing:** [Official Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py) mention "raw sequencing data optionally in BAM/CRAM" and "phenotypic data as standardized HPO terms." The public tree shows FASTQ + a DOCX, not BAM/CRAM, and not a standalone HPO file.
3. **Secondary findings vs evaluator code:** FAQ text says secondary/incidental findings will not hurt automated scores; `evaluation.py` does not filter on `finding_type` when computing rank points or F-max.
4. **Partial-match comment vs code:** the evaluator comment says "exactly one of the two true variants"; the code uses set intersection (`row.variants & true_variants`).
5. **Submission window copy:** the dataset card lists 24 August–24 October 2026; the Sage blog lists submissions opening 25 August. Both close 24 October 23:59 UTC in the Overview timeline.
