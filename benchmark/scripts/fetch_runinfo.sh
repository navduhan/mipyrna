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
TMP_GDS="${OUTDIR}/${ACCESSION}.gds.xml"
TMP_ELINK="${OUTDIR}/${ACCESSION}.elink.xml"
OUT_RUNINFO="${OUTDIR}/${ACCESSION}.runinfo.csv"
KEEP_XML="${KEEP_XML:-0}"

curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
  --data-urlencode "db=sra" \
  --data-urlencode "term=${ACCESSION}[All Fields]" \
  --data-urlencode "retmax=5000" \
  -o "${TMP_ES}"

ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_ES}" | sed -E 's#</?Id>##g' | paste -sd, -)"

# Fallback for GEO series accessions: GSE -> GDS -> linked SRA IDs.
if [[ -z "${ID_LIST}" && "${ACCESSION}" == GSE* ]]; then
  echo "Direct SRA lookup failed for ${ACCESSION}; trying GEO->SRA link fallback..."
  curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
    --data-urlencode "db=gds" \
    --data-urlencode "term=${ACCESSION}[Accession]" \
    --data-urlencode "retmax=20" \
    -o "${TMP_GDS}"

  GDS_ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_GDS}" | sed -E 's#</?Id>##g' | paste -sd, -)"
  if [[ -n "${GDS_ID_LIST}" ]]; then
    curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi" \
      --data-urlencode "dbfrom=gds" \
      --data-urlencode "db=sra" \
      --data-urlencode "id=${GDS_ID_LIST}" \
      -o "${TMP_ELINK}"
    ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_ELINK}" | sed -E 's#</?Id>##g' | paste -sd, -)"
  fi
fi

if [[ -z "${ID_LIST}" ]]; then
  echo "No SRA records found for ${ACCESSION} (direct or GEO-linked lookup)"
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
  rm -f "${TMP_ES}" "${TMP_SUMMARY}" "${TMP_GDS}" "${TMP_ELINK}"
fi
echo "Saved: ${OUT_RUNINFO}"
