# Data policy

This is **project policy** for this repository. It is not legal advice and does not replace the hackathon organizers' terms, dataset licenses, or any agreement you accepted to access the data.

## Controlled data stay local

Genomic files, alignments, variant calls, phenotype tables, and any other patient-level hackathon data must remain on private local storage (or other access-controlled locations approved for this work).

They must not be:

- committed to Git;
- pushed to GitHub or any other remote;
- added with Git LFS;
- pasted into issues, pull requests, chat logs, or documentation;
- uploaded to public object storage, public Hugging Face repos, or similar.

Prefer keeping raw inputs **outside** this repository (for example a sibling directory or a dedicated disk path). If a local `data/` directory is used, it is gitignored and still must not be force-added.

## Raw genomic data must never be committed

This includes, at minimum, FASTQ, BAM/SAM/CRAM, VCF/BCF/gVCF, and typical index companions (BAI, CRAI, TBI, CSI). Intermediate files that still contain patient-level genotypes or reads are treated the same way.

## Credentials and tokens must never be committed

Do not commit `.env` files (except the empty `.env.example` template), API tokens, SSH keys, cloud credentials, or Hugging Face / GitHub tokens.

If a secret is committed by accident, treat it as compromised: rotate it, and follow GitHub's process for removing sensitive data from history. Do not assume a later commit that "deletes" the file is enough.

## What may be published

Generated artifacts may be committed only when they are:

- aggregate or otherwise non-identifying;
- compatible with hackathon rules and the dataset terms;
- free of embedded sequences, genotypes, identifiers, and secrets.

When in doubt, keep the file local.

## Review every output before publication

`.gitignore` is **not** sufficient protection. It cannot catch force-added files, new extensions, secrets inside text, or large files that were never listed.

Before every commit:

1. Review `git status` and `git diff`.
2. Run `python scripts/check_repo_safety.py`.
3. Confirm that no path in the commit points at controlled inputs.

The safety checker is a guardrail. It does not certify that a file is legally shareable. Human review is required.

## Genome-carrying material must be deleted

This restates **organizer policy**, not additional legal advice. Sage's privacy clarification ([Community discussion #2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2)) and the [Official Rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/rules.py) require that material carrying the child's genome be deleted according to challenge rules.

Delete (examples from that organizer clarification):

- VCF, BAM, CRAM, FASTQ, and copies, subsets, slices, reformats, and index files;
- genome-scale intermediate/derivative files (filtered or annotated variant tables at genotype scale, caches, notebook state holding them);
- stored prompts or logs on systems we control that contain pasted blocks of variant data;
- model weights, embeddings, or fine-tunes trained directly on the raw genomic data.

The Official Rules additionally require deletion from environments under participant control within 30 days of hackathon close, with email confirmation.

## Permitted durable research outputs

Organizers stated that findings are meant to be published (submissions are CC BY), including things such as:

- ranked candidate variants (the Track 1 submission);
- HPO terms;
- gene/pathway rankings;
- mechanism analyses;
- code and reports;

provided prohibited genome-scale data are not embedded. A handful of named variants in a report is a finding; a genome-wide table of the child's genotypes is the dataset in another format.

## Third-party hosted services

If controlled data are ever sent to a third-party hosted service (LLM API, cloud aligner, annotation API, and similar), organizer terms require treating the service as a **processor**, not a **recipient**:

- no training on inputs or outputs, and no provider rights to reuse the content;
- retention limited in time and purpose (short-lived operational logs can be acceptable; indefinite keep-or-learn terms are not);
- record the provider, plan/tier, and data-handling settings in the methods.

Credit/grant terms can override account defaults. Opt out of training and do not rate outputs if that would re-enable learning from content.

**Default for this project:** keep controlled genomic and clinical contents out of AI assistant context unless a specific provider/tier/settings combination has been explicitly verified to meet those organizer conditions.

See [track1_challenge_spec.md](track1_challenge_spec.md) for source links.
