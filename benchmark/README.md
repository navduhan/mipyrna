# miPyRNA Benchmark Suite

This folder provides a reproducible benchmark framework to compare miPyRNA with miRDeep2 on public plant small RNA-seq studies.

## Contents

- `datasets_plant_mirna.tsv`: curated benchmark candidates with accession metadata.
- `scripts/fetch_runinfo.sh`: resolve run-level metadata from NCBI (SRA RunInfo CSV).
- `scripts/download_sra_fastq.sh`: download SRA runs (prefetch + fasterq-dump).
- `scripts/download_arabidopsis_reference.sh`: download Arabidopsis TAIR10 genome/annotation/cDNA.
- `scripts/run_benchmark_template.sh`: template benchmark runner.

## Recommended initial benchmark panel

1. `GSE13605` (Arabidopsis, AGO1-dependent small RNAs; 6 libraries).
2. `GSE12037` (Arabidopsis, AGO-related small RNAs; classic early benchmark set).
3. `PRJNA653584` (Arabidopsis high-light acclimation, modern Frontiers Plant Science dataset with multiple timepoints).

These three give:
- legacy/classic Arabidopsis datasets used widely in small RNA method papers,
- a modern stress-response dataset closer to current publication expectations.

## Quick start

```bash
cd benchmark
bash scripts/fetch_runinfo.sh GSE13605 metadata
bash scripts/fetch_runinfo.sh GSE12037 metadata
bash scripts/fetch_runinfo.sh PRJNA653584 metadata
```

Download FASTQ files from RunInfo:

```bash
bash scripts/download_sra_fastq.sh metadata/GSE13605.runinfo.csv fastq/GSE13605 8
```

Download Arabidopsis TAIR10 references:

```bash
bash scripts/download_arabidopsis_reference.sh references/arabidopsis_tair10
```

Install environment (includes `miRDeep2`):

```bash
conda env create -f ../mipyrna_environment.yaml
conda activate mipyrna-0.2
```

Then run benchmark template after downloading FASTQ files:

```bash
bash scripts/run_benchmark_template.sh \
  /path/to/project \
  benchmark/metadata/GSE13605.runinfo.csv \
  benchmark/references/arabidopsis_tair10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  ath \
  /path/to/mirdeep2_dir \
  local
```

Run with Slurm:

```bash
bash scripts/run_benchmark_template.sh \
  /path/to/project \
  benchmark/metadata/GSE13605.runinfo.csv \
  benchmark/references/arabidopsis_tair10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  ath \
  /path/to/mirdeep2_dir \
  slurm
```

## Notes

- Keep raw FASTQ files outside git; only store run manifests, logs, and summary metrics in this repo.
- For manuscript claims, report overlap, precision/recall, and rank correlation between miPyRNA and miRDeep2 outputs.
