# miPyRNA Benchmark Suite

This folder provides a reproducible benchmark framework to compare miPyRNA with miRDeep2 on public plant small RNA-seq studies.

## Contents

- `datasets_plant_mirna.tsv`: curated plant benchmark candidates with accession metadata.
- `datasets_animal_mirna.tsv`: curated animal benchmark candidates with accession metadata.
- `scripts/fetch_runinfo.sh`: resolve run-level metadata from NCBI (SRA RunInfo CSV).
- `scripts/download_sra_fastq.sh`: download SRA runs (prefetch + fasterq-dump).
- `scripts/download_arabidopsis_reference.sh`: download Arabidopsis TAIR10 genome/annotation/cDNA.
- `scripts/download_mouse_reference.sh`: download Mouse GRCm39 genome/annotation/cDNA.
- `scripts/run_benchmark_template.sh`: template benchmark runner.

## Recommended benchmark panel

1. `GSE13605` (Arabidopsis, AGO1-dependent small RNAs; 6 libraries).
2. `PRJNA1102430` (Arabidopsis, 2026 NovaSeq 6000 `Col-0` vs `amiR-mta` sRNA-seq; 6 libraries, pre-trimmed).
3. `GSE299299` / `PRJNA1273867` (Mouse MSCs, 2026 `DICER1` `EV`, `WT`, `D1320A`, `E1705K` miRNA-seq; 12 libraries, pre-trimmed).
4. `GSE12037` (Arabidopsis, AGO-related small RNAs; classic early benchmark set).
5. `PRJNA653584` (Arabidopsis high-light acclimation, Frontiers Plant Science dataset with multiple timepoints).

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
conda activate mipyrna
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

# Or set partition explicitly as the 7th argument
bash scripts/run_benchmark_template.sh \
  /path/to/project \
  benchmark/metadata/GSE13605.runinfo.csv \
  benchmark/references/arabidopsis_tair10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  ath \
  /path/to/mirdeep2_dir \
  slurm \
  compute
```

## Notes

- Keep raw FASTQ files outside git; only store run manifests, logs, and summary metrics in this repo.
- For manuscript claims, report overlap, precision/recall, and rank correlation between miPyRNA and miRDeep2 outputs.
