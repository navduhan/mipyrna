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
TMP_GDS="${OUTDIR}/${ACCESSION}.gds.xml"
TMP_ELINK="${OUTDIR}/${ACCESSION}.elink.xml"
OUT_RUNINFO="${OUTDIR}/${ACCESSION}.runinfo.csv"
KEEP_XML="${KEEP_XML:-0}"

ncbi_curl() {
  local max_attempts=5
  local attempt=1
  local sleep_s=2
  while true; do
    if curl "$@"; then
      return 0
    fi
    if [[ ${attempt} -ge ${max_attempts} ]]; then
      return 1
    fi
    sleep "${sleep_s}"
    attempt=$((attempt + 1))
    sleep_s=$((sleep_s * 2))
  done
}

ncbi_curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
  --data-urlencode "db=sra" \
  --data-urlencode "term=${ACCESSION}[All Fields]" \
  --data-urlencode "retmax=5000" \
  -o "${TMP_ES}"

ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_ES}" 2>/dev/null | sed -E 's#</?Id>##g' | paste -sd, - || true)"

# Fallback for GEO series accessions: GSE -> GDS -> linked SRA IDs.
if [[ -z "${ID_LIST}" && "${ACCESSION}" == GSE* ]]; then
  echo "Direct SRA lookup failed for ${ACCESSION}; trying GEO->SRA link fallback..."
  ncbi_curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
    --data-urlencode "db=gds" \
    --data-urlencode "term=${ACCESSION}[Accession]" \
    --data-urlencode "retmax=20" \
    -o "${TMP_GDS}"

  GDS_ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_GDS}" 2>/dev/null | sed -E 's#</?Id>##g' | paste -sd, - || true)"
  if [[ -n "${GDS_ID_LIST}" ]]; then
    ncbi_curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi" \
      --data-urlencode "dbfrom=gds" \
      --data-urlencode "db=sra" \
      --data-urlencode "id=${GDS_ID_LIST}" \
      -o "${TMP_ELINK}"
    ID_LIST="$(grep -oE '<Id>[0-9]+</Id>' "${TMP_ELINK}" 2>/dev/null | sed -E 's#</?Id>##g' | paste -sd, - || true)"
  fi
fi

if [[ -z "${ID_LIST}" ]]; then
  echo "No SRA records found for ${ACCESSION} (direct or GEO-linked lookup)"
  exit 2
fi

ncbi_curl -fsSLG "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
  --data-urlencode "db=sra" \
  --data-urlencode "id=${ID_LIST}" \
  --data-urlencode "rettype=runinfo" \
  --data-urlencode "retmode=text" \
  -o "${OUT_RUNINFO}"

if ! grep -qE '^SRR[0-9]+' "${OUT_RUNINFO}" 2>/dev/null; then
  echo "RunInfo retrieval failed for ${ACCESSION}; output did not contain SRR rows"
  exit 3
fi
if [[ "${KEEP_XML}" != "1" ]]; then
  rm -f "${TMP_ES}" "${TMP_GDS}" "${TMP_ELINK}"
fi
echo "Saved: ${OUT_RUNINFO}"
