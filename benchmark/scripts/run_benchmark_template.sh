#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <PROJECT_ROOT> <RUNINFO_CSV> <GENOME_FA> <SPECIES_CODE> <MIRDEEP2_BIN_DIR> [local|slurm] [slurm_partition]"
  exit 1
fi

PROJECT_ROOT="$1"
RUNINFO="$2"
GENOME="$3"
SPECIES="$4"
MIRDEEP2_BIN_DIR="$5"
RUN_MODE="${6:-local}"
SLURM_PARTITION="${7:-${SLURM_PARTITION:-}}"
PYTHON_BIN="${PYTHON_BIN:-${MIRDEEP2_BIN_DIR%/}/python}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python executable not found at ${PYTHON_BIN}"
  echo "Set PYTHON_BIN explicitly or pass a valid conda env bin path as arg5."
  exit 6
fi

OUTDIR="${PROJECT_ROOT}/benchmark/results/$(basename "${RUNINFO}" .runinfo.csv)"
mkdir -p "${OUTDIR}/logs" "${OUTDIR}/fastq"

echo "[1/5] Preparing sample sheet from RunInfo"
"${PYTHON_BIN}" - <<'PY' "${RUNINFO}" "${OUTDIR}/samples.txt"
import sys
import csv

runinfo = sys.argv[1]
out = sys.argv[2]

rep_map = {}
with open(out, "w") as fh:
    fh.write("# benchmark sample sheet\n")
    fh.write("SampleName\tReplication\tIdentifier\tFile1\tFile2\n")
    with open(runinfo, newline="") as inf:
        reader = csv.DictReader(inf)
        for r in reader:
            srr = str(r.get("Run", "")).strip()
            if not srr:
                continue
            condition = str(r.get("LibraryName", r.get("SampleName", "group"))).replace(" ", "_")
            rep_map.setdefault(condition, 0)
            rep_map[condition] += 1
            rep = f"{condition}_rep{rep_map[condition]}"
            fh.write(f"{rep}\t{rep}\t{condition}\t{srr}.fastq.gz\t\n")
PY

echo "[2/5] Download FASTQ runs with fasterq-dump/prefetch (user fills this step)"
echo "       RunInfo: ${RUNINFO}"
echo "       FASTQ dir: ${OUTDIR}/fastq"

echo "[3/5] Run miPyRNA workflow"
if [[ "${RUN_MODE}" == "slurm" ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "RUN_MODE=slurm requested but sbatch was not found"
    exit 4
  fi
  if [[ -z "${SLURM_PARTITION}" ]]; then
    echo "RUN_MODE=slurm requested but no partition provided."
    echo "Pass partition as 7th argument or set SLURM_PARTITION env var."
    exit 5
  fi
cat > "${OUTDIR}/logs/mipyrna.workflow.sbatch.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT}"
"${PYTHON_BIN}" -m mipyrna workflow \\
  --input-file "${OUTDIR}/samples.txt" \\
  --input-path "${OUTDIR}/fastq" \\
  --genome "${GENOME}" \\
  --species "${SPECIES}" \\
  --species-type plants \\
  --outdir "${OUTDIR}/mipyrna" \\
  --skip-qc \\
  --run-enrichment
EOF
  chmod +x "${OUTDIR}/logs/mipyrna.workflow.sbatch.sh"
  sbatch -p "${SLURM_PARTITION}" -J mipyrna_benchmark -o "${OUTDIR}/logs/mipyrna.workflow.out" -e "${OUTDIR}/logs/mipyrna.workflow.err" "${OUTDIR}/logs/mipyrna.workflow.sbatch.sh"
else
  (
    cd "${PROJECT_ROOT}"
    "${PYTHON_BIN}" -m mipyrna workflow \
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
fi

echo "[4/5] Run miRDeep2 baseline (template command; adjust mapper/pl options)"
echo "       Example:"
echo "       ${MIRDEEP2_BIN_DIR}/mapper.pl reads.fastq -e -h -j -m -l 18 -s reads_collapsed.fa -t reads_vs_genome.arf -p bowtie_index"
echo "       ${MIRDEEP2_BIN_DIR}/miRDeep2.pl reads_collapsed.fa genome.fa reads_vs_genome.arf mature.fa none hairpin.fa -t plant"
if [[ "${RUN_MODE}" == "slurm" ]]; then
  echo "       Tip: submit these commands as an sbatch job for direct cluster comparison."
fi

echo "[5/5] Save comparison metrics"
echo "       Compare:"
echo "       - known miRNA overlap"
echo "       - novel candidate overlap"
echo "       - precision/recall/F1 vs known annotations"
echo "       - score/rank correlation"

echo "Benchmark scaffold complete: ${OUTDIR}"
