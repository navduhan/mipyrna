#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <ACCESSION> <OUTDIR>"
  echo "Example: $0 GSE13605 benchmark/metadata"
  exit 1
fi

ACCESSION="$1"
OUTDIR="$2"
mkdir -p "${OUTDIR}"

TMP_ES="${OUTDIR}/${ACCESSION}.esearch.xml"
TMP_SUMMARY="${OUTDIR}/${ACCESSION}.esummary.xml"
OUT_RUNINFO="${OUTDIR}/${ACCESSION}.runinfo.csv"
KEEP_XML="${KEEP_XML:-0}"

curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
  --data-urlencode "db=sra" \
  --data-urlencode "term=${ACCESSION}[All Fields]" \
  --data-urlencode "retmax=5000" \
  -o "${TMP_ES}"

ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_ES}" | sed -E 's#</?Id>##g' | paste -sd, -)"
if [[ -z "${ID_LIST}" ]]; then
  echo "No SRA records found for ${ACCESSION}"
  exit 2
fi

curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=${ID_LIST}&retmode=xml" -o "${TMP_SUMMARY}"

RUNS="$(grep -oE 'SRR[0-9]+' "${TMP_SUMMARY}" | sort -u | paste -sd, -)"
if [[ -z "${RUNS}" ]]; then
  echo "No SRR runs resolved for ${ACCESSION}"
  exit 3
fi

curl -fsSL "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc=${RUNS}" -o "${OUT_RUNINFO}"
if [[ "${KEEP_XML}" != "1" ]]; then
  rm -f "${TMP_ES}" "${TMP_SUMMARY}"
fi
echo "Saved: ${OUT_RUNINFO}"
