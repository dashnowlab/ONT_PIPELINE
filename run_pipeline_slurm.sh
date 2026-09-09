#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=23:00:00
#SBATCH --partition=acpu
#SBATCH --qos cpu-normal
#SBATCH -o ./nf-runner-sv.stdout
#SBATCH -e ./nf-runner-sv.stderr
#SBATCH -J nf-long-read-sv
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user=elbay.aliyev@cuanschutz.edu

module purge
module load nextflow
module load singularity/3.7.4
module load miniforge/24.11.3-0

# Resolve repository-owned paths from this script, regardless of where sbatch
# was submitted. Absolute workflow paths and remote Nextflow identifiers are
# passed through unchanged.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -eq 0 ]]; then
  echo "ERROR: a Nextflow workflow is required" >&2
  exit 1
fi

PIPELINE="$1"
shift

if [[ "${PIPELINE}" != /* && -f "${REPO_DIR}/${PIPELINE}" ]]; then
  PIPELINE="${REPO_DIR}/${PIPELINE}"
fi

# ── Parse --work_dir from the pipeline arguments ──────────────────────────────
# All pipeline outputs go to --work_dir; the Nextflow work dir is always
# <work_dir>/work and is appended automatically — no need to pass -w manually.
WORK_DIR=""
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case "${args[$i]}" in
    --work_dir)    WORK_DIR="${args[$((i+1))]}" ;;
    --work_dir=*)  WORK_DIR="${args[$i]#--work_dir=}" ;;
  esac
done

if [[ -z "${WORK_DIR}" ]]; then
  echo "ERROR: --work_dir is required" >&2
  exit 1
fi

# Strip trailing slash for consistency
WORK_DIR="${WORK_DIR%/}"
NXF_WORK="${WORK_DIR}/work"

mkdir -p "${NXF_WORK}"

echo "Pipeline output dir : ${WORK_DIR}"
echo "Nextflow work dir   : ${NXF_WORK}"

export NXF_HOME="${WORK_DIR}/.nextflow"
export NXF_CACHE_DIR="${WORK_DIR}/.nextflow"

echo "== Starting Long-Read SV Pipeline =="
cd "${REPO_DIR}"
nextflow run "${PIPELINE}" "$@" -w "${NXF_WORK}" -profile singularity
echo "== Nextflow complete =="

# ─────────────────────────────────────────────────────────────────────────────
# Example usage (no -w needed — it is derived from --work_dir automatically):
#
#   sbatch run_pipeline_slurm.sh long_read_sv.nf \
#       --bam_list      resources/samples/test_bam.list \
#       --work_dir      /pl/active/dashnowlab/work/ealiyev/SV/SHAIKH_SV \
#       --projectname   SHAIKH_SV \
#       --ref           /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
#       --tr_bed        /gpfs/data_jrnas1/ref_data/Hsapiens/hg38/ONT/human_GRCh38_no_alt_analysis_set.trf.bed \
#       --software_path /gpfs/projects/tmedicine/KFakhroLAB/SV/software \
#       --caller        sniffles
# ─────────────────────────────────────────────────────────────────────────────
