#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-benchmark/references/mouse_grcm39}"
mkdir -p "${OUTDIR}"

# Ensembl mirror for Mus musculus GRCm39 references (primary assembly avoids alt-scaffold duplicate mappings).
GENOME_URL_DEFAULT="https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
CDNA_URL_DEFAULT="https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
ENSEMBL_PUB_URL="https://ftp.ensembl.org/pub/"

GENOME_URL="${GENOME_URL:-${GENOME_URL_DEFAULT}}"
CDNA_URL="${CDNA_URL:-${CDNA_URL_DEFAULT}}"
GTF_URL="${GTF_URL:-}"

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

resolve_latest_gtf_url() {
  if command -v curl >/dev/null 2>&1; then
    local rel
    rel="$(curl -fsSL "${ENSEMBL_PUB_URL}" | grep -oE 'release-[0-9]+/' | sort -V | tail -n 1 || true)"
    if [[ -n "${rel}" ]]; then
      local gtf_base="${ENSEMBL_PUB_URL}${rel}gtf/mus_musculus/"
      local file
      file="$(curl -fsSL "${gtf_base}" | grep -oE 'Mus_musculus\.GRCm39\.[0-9]+\.gtf\.gz' | sort -V | tail -n 1 || true)"
      if [[ -n "${file}" ]]; then
        printf "%s%s" "${gtf_base}" "${file}"
        return 0
      fi
    fi
  fi
  return 1
}

if [[ -z "${GTF_URL}" ]]; then
  GTF_URL="$(resolve_latest_gtf_url || true)"
fi
if [[ -z "${GTF_URL}" ]]; then
  echo "Could not auto-resolve GTF URL. Set GTF_URL env var and rerun."
  exit 3
fi

echo "Downloading Mus musculus GRCm39 references into ${OUTDIR}"
download_file "${GENOME_URL}" "${OUTDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
download_file "${GTF_URL}" "${OUTDIR}/Mus_musculus.GRCm39.gtf.gz"
download_file "${CDNA_URL}" "${OUTDIR}/Mus_musculus.GRCm39.cdna.all.fa.gz"

gunzip -fk "${OUTDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
gunzip -fk "${OUTDIR}/Mus_musculus.GRCm39.gtf.gz"
gunzip -fk "${OUTDIR}/Mus_musculus.GRCm39.cdna.all.fa.gz"

echo "Reference files ready:"
echo "  Genome: ${OUTDIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
echo "  GTF:    ${OUTDIR}/Mus_musculus.GRCm39.gtf"
echo "  cDNA:   ${OUTDIR}/Mus_musculus.GRCm39.cdna.all.fa"
