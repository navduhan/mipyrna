#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <PROJECT_ROOT> <RUNINFO_CSV> <GENOME_FA> <SPECIES_CODE> <MIRDEEP2_BIN_DIR>"
  exit 1
fi

PROJECT_ROOT="$1"
RUNINFO="$2"
GENOME="$3"
SPECIES="$4"
MIRDEEP2_BIN_DIR="$5"

OUTDIR="${PROJECT_ROOT}/benchmark/results/$(basename "${RUNINFO}" .runinfo.csv)"
mkdir -p "${OUTDIR}/logs" "${OUTDIR}/fastq"

echo "[1/5] Preparing sample sheet from RunInfo"
python3 - <<'PY' "${RUNINFO}" "${OUTDIR}/samples.txt"
import sys
import pandas as pd

runinfo = pd.read_csv(sys.argv[1])
out = sys.argv[2]

cols = ["SampleName", "Replication", "Identifier", "File1", "File2"]
rows = []
rep_map = {}
for _, r in runinfo.iterrows():
    srr = str(r.get("Run", "")).strip()
    if not srr:
        continue
    condition = str(r.get("LibraryName", r.get("SampleName", "group"))).replace(" ", "_")
    rep_map.setdefault(condition, 0)
    rep_map[condition] += 1
    rep = f"{condition}_rep{rep_map[condition]}"
    rows.append([rep, rep, condition, f"{srr}.fastq.gz", ""])

df = pd.DataFrame(rows, columns=cols)
with open(out, "w") as fh:
    fh.write("# benchmark sample sheet\n")
df.to_csv(out, sep="\t", index=False, mode="a")
PY

echo "[2/5] Download FASTQ runs with fasterq-dump/prefetch (user fills this step)"
echo "       RunInfo: ${RUNINFO}"
echo "       FASTQ dir: ${OUTDIR}/fastq"

echo "[3/5] Run miPyRNA workflow"
(
  cd "${PROJECT_ROOT}"
  python3 -m mipyrna workflow \
    --input-file "${OUTDIR}/samples.txt" \
    --input-path "${OUTDIR}/fastq" \
    --genome "${GENOME}" \
    --species "${SPECIES}" \
    --species-type plants \
    --outdir "${OUTDIR}/mipyrna" \
    --skip-qc \
    --run-enrichment \
    > "${OUTDIR}/logs/mipyrna.workflow.log" 2>&1
)

echo "[4/5] Run miRDeep2 baseline (template command; adjust mapper/pl options)"
echo "       Example:"
echo "       ${MIRDEEP2_BIN_DIR}/mapper.pl reads.fastq -e -h -j -m -l 18 -s reads_collapsed.fa -t reads_vs_genome.arf -p bowtie_index"
echo "       ${MIRDEEP2_BIN_DIR}/miRDeep2.pl reads_collapsed.fa genome.fa reads_vs_genome.arf mature.fa none hairpin.fa -t plant"

echo "[5/5] Save comparison metrics"
echo "       Compare:"
echo "       - known miRNA overlap"
echo "       - novel candidate overlap"
echo "       - precision/recall/F1 vs known annotations"
echo "       - score/rank correlation"

echo "Benchmark scaffold complete: ${OUTDIR}"
