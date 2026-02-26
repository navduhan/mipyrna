#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <RUNINFO_CSV> <OUTDIR> [threads]"
  echo "Example: $0 benchmark/metadata/GSE13605.runinfo.csv benchmark/fastq/GSE13605 8"
  exit 1
fi

RUNINFO="$1"
OUTDIR="$2"
THREADS="${3:-8}"

if [[ ! -f "${RUNINFO}" ]]; then
  echo "RunInfo file not found: ${RUNINFO}"
  exit 2
fi

if ! command -v prefetch >/dev/null 2>&1; then
  echo "prefetch not found in PATH (install sra-tools)"
  exit 3
fi
if ! command -v fasterq-dump >/dev/null 2>&1; then
  echo "fasterq-dump not found in PATH (install sra-tools)"
  exit 4
fi

mkdir -p "${OUTDIR}" "${OUTDIR}/sra"

# Extract SRR accessions from the first CSV column "Run".
mapfile -t RUNS < <(tail -n +2 "${RUNINFO}" | cut -d',' -f1 | sed '/^$/d')

if [[ ${#RUNS[@]} -eq 0 ]]; then
  echo "No run accessions found in ${RUNINFO}"
  exit 5
fi

for SRR in "${RUNS[@]}"; do
  echo "Downloading ${SRR}"
  prefetch "${SRR}" --output-directory "${OUTDIR}/sra"

  SRA_PATH="${OUTDIR}/sra/${SRR}/${SRR}.sra"
  if [[ ! -f "${SRA_PATH}" ]]; then
    ALT_PATH="${OUTDIR}/sra/${SRR}.sra"
    if [[ -f "${ALT_PATH}" ]]; then
      SRA_PATH="${ALT_PATH}"
    else
      echo "SRA file not found for ${SRR}"
      continue
    fi
  fi

  # Use split-files to support both SINGLE and PAIRED runs.
  fasterq-dump "${SRA_PATH}" --split-files --threads "${THREADS}" --outdir "${OUTDIR}"

  if [[ -f "${OUTDIR}/${SRR}.fastq" ]]; then
    gzip -f "${OUTDIR}/${SRR}.fastq"
  fi
  if [[ -f "${OUTDIR}/${SRR}_1.fastq" ]]; then
    gzip -f "${OUTDIR}/${SRR}_1.fastq"
  fi
  if [[ -f "${OUTDIR}/${SRR}_2.fastq" ]]; then
    gzip -f "${OUTDIR}/${SRR}_2.fastq"
  fi
done

echo "FASTQ download complete in ${OUTDIR}"
