#!/bin/bash -l
#SBATCH --job-name=ont_align
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

###############################################################################
# Usage
#
# sbatch ont_alignment.sh \
#   --sample-dir /path/to/fastq_directory \
#   --reference /path/to/reference.fasta \
#   --output-dir /path/to/output \
#   --threads 16
###############################################################################

usage() {
    cat <<EOF
Usage:
  sbatch $0 --sample-dir DIR --reference FASTA --output-dir DIR [--threads N]

Required:
  --sample-dir DIR    Directory containing ONT FASTQ files
  --reference FASTA  Reference genome FASTA
  --output-dir DIR   Root directory for results, work files, and caches

Optional:
  --threads N        Alignment threads (default: 16)
  -h, --help         Show this help message
EOF
}

###############################################################################
# Input
###############################################################################

SAMPLE_DIR=""
REFERENCE_FASTA=""
OUTPUT_DIR=""
ALIGNMENT_THREADS=16

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --sample-dir requires a value" >&2; exit 1; }
            SAMPLE_DIR="$2"
            shift 2
            ;;
        --reference)
            [[ $# -ge 2 ]] || { echo "ERROR: --reference requires a value" >&2; exit 1; }
            REFERENCE_FASTA="$2"
            shift 2
            ;;
        --output-dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --output-dir requires a value" >&2; exit 1; }
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            [[ $# -ge 2 ]] || { echo "ERROR: --threads requires a value" >&2; exit 1; }
            ALIGNMENT_THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$SAMPLE_DIR" || -z "$REFERENCE_FASTA" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --sample-dir, --reference, and --output-dir are required" >&2
    usage >&2
    exit 1
fi

if [[ ! "$ALIGNMENT_THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer" >&2
    exit 1
fi

if [[ ! -d "$SAMPLE_DIR" ]]; then
    echo "ERROR: Sample directory does not exist:"
    echo "  $SAMPLE_DIR"
    exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA does not exist:"
    echo "  $REFERENCE_FASTA"
    exit 1
fi

SAMPLE_DIR=$(realpath "$SAMPLE_DIR")
REFERENCE_FASTA=$(realpath "$REFERENCE_FASTA")
SAMPLE_ROOT=$(realpath -m "$OUTPUT_DIR")
SAMPLE_NAME=$(basename "$SAMPLE_DIR")

###############################################################################
# Main settings
###############################################################################

PIPELINE="epi2me-labs/wf-alignment"

###############################################################################
# Per-sample paths
###############################################################################

FASTQ_INPUT_DIR="${SAMPLE_ROOT}/fastq_input"
REFERENCE_INPUT_DIR="${SAMPLE_ROOT}/reference_input"

RESULTS_DIR="${SAMPLE_ROOT}/results"
WORK_DIR="${SAMPLE_ROOT}/work"
TMP_DIR="${SAMPLE_ROOT}/tmp"

NXF_HOME_DIR="${SAMPLE_ROOT}/.nextflow_home"
NXF_ASSETS_DIR="${SAMPLE_ROOT}/assets"

RUN_CONFIG="${SAMPLE_ROOT}/nextflow.config"
FASTQ_LIST="${SAMPLE_ROOT}/fastq_files.txt"

###############################################################################
# Cache paths
###############################################################################

NEXTFLOW_SINGULARITY_CACHE="${SAMPLE_ROOT}/.cache/nextflow_singularity_images"
SINGULARITY_OCI_CACHE="${SAMPLE_ROOT}/.cache/singularity_oci_cache"

###############################################################################
# Initialize modules
###############################################################################

source /etc/profile.d/modules.sh 2>/dev/null || true

module purge
module load nextflow
module load singularity/3.7.4

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: Nextflow is unavailable after module load."
    exit 1
fi

if ! command -v singularity >/dev/null 2>&1; then
    echo "ERROR: Singularity is unavailable after module load."
    exit 1
fi

echo "Nextflow version:"
nextflow -version

echo
echo "Singularity version:"
singularity --version

###############################################################################
# Create directories
###############################################################################

mkdir -p \
    "$SAMPLE_ROOT" \
    "$FASTQ_INPUT_DIR" \
    "$REFERENCE_INPUT_DIR" \
    "$RESULTS_DIR" \
    "$WORK_DIR" \
    "$TMP_DIR" \
    "$NXF_HOME_DIR" \
    "$NXF_ASSETS_DIR" \
    "$NEXTFLOW_SINGULARITY_CACHE" \
    "$SINGULARITY_OCI_CACHE"

###############################################################################
# Environment variables
###############################################################################

export TMPDIR="$TMP_DIR"

export NXF_HOME="$NXF_HOME_DIR"
export NXF_ASSETS="$NXF_ASSETS_DIR"
export NXF_WORK="$WORK_DIR"

export NXF_SINGULARITY_CACHEDIR="$NEXTFLOW_SINGULARITY_CACHE"
export SINGULARITY_CACHEDIR="$SINGULARITY_OCI_CACHE"


###############################################################################
# Remove old partial container pulls
###############################################################################

find "$NEXTFLOW_SINGULARITY_CACHE" \
    -maxdepth 1 \
    -type f \
    -name "*.pulling.*" \
    -mmin +360 \
    -print \
    -delete || true

###############################################################################
# Find FASTQ files recursively
###############################################################################

find "$SAMPLE_DIR" \
    -type f \
    \( \
        -iname "*.fastq.gz" -o \
        -iname "*.fq.gz"    -o \
        -iname "*.fastq"    -o \
        -iname "*.fq" \
    \) \
    -print0 \
    | sort -z \
    | tr '\0' '\n' > "$FASTQ_LIST"

FASTQ_COUNT=$(wc -l < "$FASTQ_LIST")
FASTQ_COUNT=${FASTQ_COUNT//[[:space:]]/}

if [[ "$FASTQ_COUNT" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found under:"
    echo "  $SAMPLE_DIR"
    exit 1
fi

###############################################################################
# Create flat FASTQ input directory
#
# This combines all FASTQ files from directories such as:
#
#   fastq_pass_36A
#   fastq_pass_36B
#
# into one logical sample.
###############################################################################

echo
echo "Preparing flat FASTQ input directory..."

find "$FASTQ_INPUT_DIR" \
    -mindepth 1 \
    -maxdepth 1 \
    \( -type l -o -type f \) \
    -delete

COUNTER=0

while IFS= read -r FASTQ_FILE; do
    [[ -n "$FASTQ_FILE" ]] || continue

    COUNTER=$((COUNTER + 1))

    ORIGINAL_NAME=$(basename "$FASTQ_FILE")
    LINK_NAME=$(printf "%06d_%s" "$COUNTER" "$ORIGINAL_NAME")

    ln -s "$FASTQ_FILE" "${FASTQ_INPUT_DIR}/${LINK_NAME}"
done < "$FASTQ_LIST"

LINK_COUNT=$(
    find "$FASTQ_INPUT_DIR" \
        -mindepth 1 \
        -maxdepth 1 \
        -type l \
        | wc -l
)

LINK_COUNT=${LINK_COUNT//[[:space:]]/}

if [[ "$LINK_COUNT" -ne "$FASTQ_COUNT" ]]; then
    echo "ERROR: FASTQ symlink count does not match source FASTQ count."
    echo "Source FASTQs: $FASTQ_COUNT"
    echo "Created links: $LINK_COUNT"
    exit 1
fi

###############################################################################
# Prepare reference directory
###############################################################################

find "$REFERENCE_INPUT_DIR" \
    -mindepth 1 \
    -maxdepth 1 \
    \( -type l -o -type f \) \
    -delete

ln -s "$REFERENCE_FASTA" \
    "${REFERENCE_INPUT_DIR}/$(basename "$REFERENCE_FASTA")"

###############################################################################
# Write Nextflow configuration
###############################################################################

cat > "$RUN_CONFIG" <<CFG
report {
    enabled = true
    file = '${RESULTS_DIR}/nextflow-report.html'
    overwrite = true
}

timeline {
    enabled = true
    file = '${RESULTS_DIR}/nextflow-timeline.html'
    overwrite = true
}

trace {
    enabled = true
    file = '${RESULTS_DIR}/nextflow-trace.txt'
    overwrite = true
}

dag {
    enabled = true
    file = '${RESULTS_DIR}/nextflow-dag.html'
    overwrite = true
}

conda.enabled = false

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '${NEXTFLOW_SINGULARITY_CACHE}'
    pullTimeout = '2h'
}

process {
    executor = 'slurm'

    /*
     * Defaults for small helper tasks.
     * Alpine requires every submitted task to include a time limit.
     */
    cpus = 1
    memory = 4.GB
    time = 2.h

    queue = 'amilan'
    clusterOptions = '--qos=normal'

    /*
     * One retry after the initial failed attempt.
     *
     * attempt 1 = original run
     * attempt 2 = one retry
     * then terminate the workflow
     */
    maxRetries = 1

    errorStrategy = {
        task.attempt <= 1 ? 'retry' : 'terminate'
    }

    /*
     * FASTQ aggregation and compression.
     */
    withName: /.*fastcat.*/ {
        cpus = 4
        memory = 32.GB
        time = 24.h

        /*
         * Do not launch a duplicate long-running fastcat task if scheduler
         * polling ever fails again.
         */
        maxRetries = 0
        errorStrategy = 'terminate'
    }

    /*
     * Minimap2 reference index.
     * CPU count is controlled by --threads.
     */
    withName: /.*makeMMIndex.*/ {
        memory = 48.GB
        time = 12.h
    }

    /*
     * Main minimap2 alignment task.
     * CPU count is controlled by --threads.
     */
    withName: /.*alignReads.*/ {
        memory = 48.GB
        time = 24.h
    }

    withName: /.*bamstats.*/ {
        memory = 32.GB
        time = 24.h
    }

    withName: /.*readDepthPerRef.*/ {
        memory = 32.GB
        time = 24.h
    }

    resourceLimits = [
        cpus: 48,
        memory: 240.GB,
        time: 72.h
    ]

    withLabel: process_long {
        resourceLimits = [
            cpus: 48,
            memory: 240.GB,
            time: 168.h
        ]
    }
}

executor {
    name = 'slurm'

    /*
     * SLURM may report the task complete before the .exitcode file becomes
     * visible to the Nextflow controller on GPFS. Wait up to 30 minutes.
     */
    exitReadTimeout = '30 min'
    queueStatInterval = '30 sec'

    queueSize = 48
    perJobMemLimit = true
    submitRateLimit = '30 / 1 min'
    killBatchSize = 50
}
CFG

###############################################################################
# Display run information
###############################################################################

echo
echo "=================================================================="
echo "ONT wf-alignment job"
echo "=================================================================="
echo "Controller SLURM job:     ${SLURM_JOB_ID:-not_available}"
echo "Sample name:              $SAMPLE_NAME"
echo "Original sample folder:   $SAMPLE_DIR"
echo "FASTQ files found:        $FASTQ_COUNT"
echo "Flat FASTQ directory:     $FASTQ_INPUT_DIR"
echo "Reference FASTA:          $REFERENCE_FASTA"
echo "Reference directory:      $REFERENCE_INPUT_DIR"
echo "Alignment threads:        $ALIGNMENT_THREADS"
echo "Output root:              $SAMPLE_ROOT"
echo "Results directory:        $RESULTS_DIR"
echo "Nextflow work directory:  $WORK_DIR"
echo "Temporary directory:      $TMP_DIR"
echo "Nextflow assets:          $NXF_ASSETS_DIR"
echo "Nextflow home:            $NXF_HOME_DIR"
echo "SIF image cache:          $NEXTFLOW_SINGULARITY_CACHE"
echo "OCI/blob cache:           $SINGULARITY_OCI_CACHE"
echo "Nextflow config:          $RUN_CONFIG"
echo "Default retries:          1"
echo "Fastcat retries:          0 (prevents duplicate orphan jobs)"
echo "Exit-read timeout:        30 minutes"
echo "=================================================================="
echo

###############################################################################
# Run wf-alignment
###############################################################################

cd "$SAMPLE_ROOT"

nextflow run "$PIPELINE" \
    --fastq "$FASTQ_INPUT_DIR" \
    --references "$REFERENCE_INPUT_DIR" \
    --sample "$SAMPLE_NAME" \
    --threads "$ALIGNMENT_THREADS" \
    --out_dir "$RESULTS_DIR" \
    -profile singularity \
    -c "$RUN_CONFIG" \
    -w "$WORK_DIR" \
    -resume

###############################################################################
# Completion
###############################################################################

echo
echo "=================================================================="
echo "Alignment completed successfully"
echo "Sample:  $SAMPLE_NAME"
echo "Results: $RESULTS_DIR"
echo "=================================================================="
