#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-benchmark/references/arabidopsis_tair10}"
mkdir -p "${OUTDIR}"

# Ensembl Plants "current" mirror for Arabidopsis TAIR10 references.
GENOME_URL_DEFAULT="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GTF_BASE_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gtf/arabidopsis_thaliana/"
CDNA_BASE_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/cdna/"

GENOME_URL="${GENOME_URL:-${GENOME_URL_DEFAULT}}"
GTF_URL="${GTF_URL:-}"
CDNA_URL="${CDNA_URL:-}"

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

resolve_latest_url() {
  local base_url="$1"
  local pattern="$2"
  if command -v curl >/dev/null 2>&1; then
    local listing
    listing="$(curl -fsSL "${base_url}")"
    local file
    file="$(printf "%s" "${listing}" | grep -oE "${pattern}" | sort -V | tail -n 1 || true)"
    if [[ -n "${file}" ]]; then
      printf "%s%s" "${base_url}" "${file}"
      return 0
    fi
  fi
  return 1
}

if [[ -z "${GTF_URL}" ]]; then
  GTF_URL="$(resolve_latest_url "${GTF_BASE_URL}" 'Arabidopsis_thaliana\.TAIR10\.[0-9]+\.gtf\.gz' || true)"
fi
if [[ -z "${CDNA_URL}" ]]; then
  CDNA_URL="$(resolve_latest_url "${CDNA_BASE_URL}" 'Arabidopsis_thaliana\.TAIR10\.cdna\.all\.fa\.gz' || true)"
fi
if [[ -z "${GTF_URL}" ]]; then
  echo "Could not auto-resolve GTF URL. Set GTF_URL env var and rerun."
  exit 3
fi
if [[ -z "${CDNA_URL}" ]]; then
  echo "Could not auto-resolve cDNA URL. Set CDNA_URL env var and rerun."
  exit 4
fi

echo "Downloading Arabidopsis references into ${OUTDIR}"
download_file "${GENOME_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
download_file "${GTF_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.gtf.gz"
download_file "${CDNA_URL}" "${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.gtf.gz"
gunzip -fk "${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"

echo "Reference files ready:"
echo "  Genome: ${OUTDIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
echo "  GTF:    ${OUTDIR}/Arabidopsis_thaliana.TAIR10.gtf"
echo "  cDNA:   ${OUTDIR}/Arabidopsis_thaliana.TAIR10.cdna.all.fa"
