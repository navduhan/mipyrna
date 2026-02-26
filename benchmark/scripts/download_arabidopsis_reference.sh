#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-benchmark/references/arabidopsis_tair10}"
mkdir -p "${OUTDIR}"

# Ensembl Plants "current" mirror for Arabidopsis TAIR10 references.
GENOME_URL_DEFAULT="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GTF_URL_DEFAULT="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.60.gtf.gz"
CDNA_URL_DEFAULT="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

GENOME_URL="${GENOME_URL:-${GENOME_URL_DEFAULT}}"
GTF_URL="${GTF_URL:-${GTF_URL_DEFAULT}}"
CDNA_URL="${CDNA_URL:-${CDNA_URL_DEFAULT}}"

download_file() {
  local url="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    echo "Exists: ${out}"
    return
  fi
  if command -v curl >/dev/null 2>&1; then
    curl -fL "${url}" -o "${out}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${out}" "${url}"
  else
    echo "Need curl or wget to download reference files"
    exit 2
  fi
}

echo "Downloading Arabidopsis references into ${OUTDIR}"
download_file "${GENOME_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
download_file "${GTF_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.60.gtf.gz"
download_file "${CDNA_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.60.gtf.gz"
gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

echo "Reference files ready:"
echo "  Genome: ${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
echo "  GTF:    ${OUTDIR}/Arabidopsis_thaliana.TAIR10.60.gtf"
echo "  cDNA:   ${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa"
