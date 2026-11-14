# miPyRNA

**miPyRNA** is an end-to-end Python toolkit for small RNA-seq (sRNA-seq) analysis across plants and animals — integrating read quality control, adapter trimming, non-miRNA structural RNA screening (Rfam 15.1), *de novo* and reference-guided miRNA discovery, deep-learning precursor classification, quantification, differential expression, target prediction, functional enrichment, and automated PDF reporting.

---

## Key Features

- **Cross-Kingdom Precursor Models**: Bundled SHA-256-verified deep-learning classifiers trained separately on curated plant (`plants`) and animal (`animals`) hairpin structural and sequence-thermodynamic features (129 features via ViennaRNA + forgi).
- **Biogenesis-Aware *De Novo* Discovery**: Deterministic interval clustering, multi-flank hairpin folding, Dicer 5'-end consistency, star-strand duplex geometry, and dinucleotide-shuffled `randfold` empirical *p*-values.
- **Reference-Guided Known miRNA Scan (`--known-scan`)**: Resolves known mature miRNAs (miRBase, PmiREN, MirGeneDB) onto reference precursors and merges precursor-aligned detections with *de novo* calls.
- **Bundled Rfam 15.1 Screening (`--rfam-filter`)**: Screens collapsed reads offline against all non-miRNA Rfam 15.1 families with exact Bowtie matching, annotating `rfam_reads` and `Rfam` family accessions on every call (and optionally removing structural ncRNA reads prior to discovery and quantification).
- **Downstream Quantification, DE, Targets & Report**: Automated Bowtie mature quantification, `DESeq2` differential expression via `rpy2`, optional `miRanda` target prediction + `g:Profiler` functional enrichment, and a publication-ready summary PDF (`8.Report/miPyRNA_report.pdf`) plus provenance JSON (`workflow_metadata.json`).
- **Laptop-to-HPC Execution**: Native multi-core execution (`--cpu`, `--jobs`) and Slurm job-array submission (`--slurm`).

---

## Installation

### Option 1: Conda / Mamba (Recommended)

```bash
git clone https://github.com/navduhan/mipyrna.git
cd mipyrna
conda env create -f mipyrna_environment.yaml
conda activate mipyrna
pip install --no-deps .
```

Verify external binaries, Python libraries, and bundled models:

```bash
mipyrna --check
```

### Option 2: Python `venv` / `pip`

If the required external command-line tools and R packages are already available on `PATH`:

```bash
git clone https://github.com/navduhan/mipyrna.git
cd mipyrna
pip install .
mipyrna --check
```

### Option 3: Docker

```bash
git clone https://github.com/navduhan/mipyrna.git
cd mipyrna
docker build -t mipyrna .
docker run --rm mipyrna mipyrna --check
```

### External Tool Requirements

- **Required**: `bowtie`, `bowtie-build` (Bowtie 1), `makeblastdb`, `blastn` (NCBI BLAST+), `miranda`, and `R` with `DESeq2` (via `rpy2`).
- **Optional (QC & trimming)**: `fastqc`, `trimmomatic`, `flexbar`, `trim_galore`.
- **Python**: `>= 3.10`.

---

## Sample Sheet Format

Provide a tab-delimited sample sheet (`samples.txt`). Comment lines begin with `#`; leave `File2` empty for single-end sRNA-seq libraries:

```tsv
# Arabidopsis sRNA-seq study (Col-0 vs mutant)
SampleName	Replication	Identifier	File1	File2
Col0_rep1	Col0_rep1	Col-0	SRR39833997_1.fastq.gz	
Col0_rep2	Col0_rep2	Col-0	SRR39833996_1.fastq.gz	
Col0_rep3	Col0_rep3	Col-0	SRR39833994_1.fastq.gz	
mut_rep1	mut_rep1	mutant	SRR39833993_1.fastq.gz	
mut_rep2	mut_rep2	mutant	SRR39833992_1.fastq.gz	
mut_rep3	mut_rep3	mutant	SRR39833991_1.fastq.gz	
```

---

## Quick Start

### 1. End-to-End sRNA-seq Workflow

```bash
mipyrna workflow \
  --input-file samples.txt \
  --input-path /path/to/fastq \
  --genome /path/to/genome.fa \
  --species ath \
  --species-type plants \
  --outdir mipyrna_results \
  --known-scan \
  --cpu 16
```

For pre-trimmed sRNA-seq FASTQ files, add `--skip-trim`. To also run target prediction and functional enrichment:

```bash
mipyrna workflow \
  --input-file samples.txt \
  --input-path /path/to/fastq \
  --genome /path/to/genome.fa \
  --species mmu \
  --species-type animals \
  --outdir mipyrna_mouse \
  --skip-trim \
  --known-scan \
  --mrna-file /path/to/transcripts.fa \
  --run-enrichment \
  --enrichment-organism mmusculus \
  --cpu 32
```

You can also run from an INI configuration file (see `config.ini.example`):

```bash
mipyrna -c config.ini
```

### 2. Standalone Subcommands

```bash
# Score candidate precursor hairpins directly from a FASTA file
mipyrna predict --fasta candidate_precursors.fa --species-type plants --output predictions.tsv

# Run miRanda target prediction
mipyrna targets --mirna-file mature_miRNAs.fa --mrna-file transcripts.fa --outdir results

# Run g:Profiler functional enrichment on predicted targets
mipyrna enrich --targets-file results/targets_raw/miranda_targets.tsv --organism athaliana --outdir results

# Regenerate the PDF analysis report from an existing run directory
mipyrna report --outdir mipyrna_results
```

---

## Key Workflow Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `--input-file`, `--input-path` | *(required)* | Tab-delimited sample sheet and directory containing raw/trimmed FASTQ files |
| `--genome` | *(required)* | Reference genome FASTA |
| `--species` | *(required)* | 3-letter organism code (e.g., `ath`, `osa`, `mmu`, `hsa`) |
| `--species-type` | `plants` | Deep-learning classifier and biogenesis rules (`plants` or `animals`) |
| `--outdir` | `mipyrna_results` | Output directory |
| `--paired` | `off` | Treat input libraries as paired-end |
| `--trimmer` | `trimmomatic` | Read trimmer (`trimmomatic`, `flexbar`, or `trim_galore`) |
| `--adapter` | auto-detect | 3' adapter sequence or kit name |
| `--skip-qc`, `--skip-trim` | `off` | Skip `FastQC` and/or adapter trimming (e.g., for pre-trimmed libraries) |
| `--rfam-filter` | `off` | Remove reads matching bundled non-miRNA Rfam 15.1 families prior to alignment (Rfam overlaps are always annotated) |
| `--known-scan` | `off` | Align reads against reference precursors (miRBase, PmiREN, MirGeneDB) and merge detected known arms |
| `--known-min-reads` | `10` | Minimum reads on a known precursor arm for `--known-scan` detection |
| `--known-5p-tolerance` | `2` | Maximum 5' offset (nt) from annotated mature start in `--known-scan` |
| `--known-unique` | `off` | Count only uniquely mapping precursor reads during `--known-scan` |
| `--require-star` / `--no-require-star` | `on` | Require star-strand support on the opposite arm for novel miRNAs |
| `--read-cutoff` | `20` | Minimum cluster read depth to fold candidate precursors |
| `--min-mature-reads`, `--min-star-reads` | `10`, `1` | Minimum mature and star read counts for candidate retention |
| `--max-inconsistent` | `0.25` | Maximum fraction of hairpin reads outside mature and star positions |
| `--min-5p-consistency` | `0.5` | Minimum fraction of mature-arm reads sharing the modal 5' cleavage site |
| `--min-novel-score` | `1.0` | Minimum biogenesis heuristic score for novel miRNAs |
| `--weighted` | `off` | Apportion multi-mapping reads fractionally across mature loci during quantification |
| `--mrna-file` | `None` | Transcript/cDNA FASTA to enable `miRanda` target prediction |
| `--run-enrichment`, `--enrichment-organism` | `off`, `None` | Run `g:Profiler` enrichment on predicted targets (e.g., `athaliana`, `mmusculus`) |
| `--cpu`, `--jobs`, `--mem` | `8`, auto, `20` | Total CPU cores, concurrent samples (`cpu // jobs` cores each), and Slurm memory (GB) |
| `--slurm`, `--partition`, `--account`, `--time` | `off` | Submit workflow and per-sample array tasks via Slurm |
| `--resume`, `--dryrun` | `off` | Resume completed stages/samples in `--outdir`, or preview external commands without running |

---

## Output Directory Layout

| Stage Directory | Key Outputs |
| :--- | :--- |
| `1.Quality_and_trimming/` | Raw and post-trimming `FastQC` reports and trimmed FASTQ files |
| `2.Alignment/` | Collapsed unique-read FASTAs (`collapsed_reads/`), Rfam 15.1 screening outputs (`rfam/`), and genome SAM alignments |
| `3.miRNA_Discovery/` | `known_miRNAs.tsv` / `.xlsx`, `novel_miRNAs.tsv` / `.xlsx`, `known_scan/`, and precursor/mature FASTAs |
| `4.Quantification/` | `Raw_read_counts.tsv` and `miRNA_mature.fa` |
| `5.Differential_Expression/` | `DESeq2_results.tsv` and differential expression plots |
| `6.Target_Prediction/` | `miRanda` target predictions (when `--mrna-file` is supplied) |
| `7.Functional_Annotation/` | `g:Profiler` enrichment tables (when `--run-enrichment` is enabled) |
| `8.Report/` | `miPyRNA_report.pdf` and `workflow_metadata.json` |

---

## Rfam 15.1 Alerts and Optional Filtering

Every workflow run screens collapsed reads offline against the bundled XZ-compressed Rfam 15.1 non-miRNA reference (`mipyrna/data/ncRNA_rfam.tar.xz`, verified by SHA-256 checksum at runtime with no network access required):

- **Default (`--rfam-filter` off)**: All reads remain available for discovery, `--known-scan`, and quantification, while `known_miRNAs.tsv`, `novel_miRNAs.tsv`, and `known_scan.tsv` report `rfam_reads` and matched `Rfam` family accessions (e.g., `RF00005:tRNA`).
- **Strict filtering (`--rfam-filter`)**: Reads matching non-miRNA Rfam families are removed prior to genome alignment, discovery, known scan, and quantification.

---

## Running Tests

```bash
python -m unittest discover -s tests -v
```

---

## License & Contact

Released under the **MIT License** (see `LICENSE`).

- **Author**: Naveen Duhan (`naveen.duhan@usu.edu` / `naveen.duhan@outlook.com`)
- **Corresponding / Lab**: Rakesh Kaundal (`rkaundal@usu.edu`), Kaundal Bioinformatics Lab, Utah State University
