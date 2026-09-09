#!/bin/bash -l
#SBATCH --job-name=ont_variation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
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
# sbatch ont_human_variation.sh \
#   --bam /path/to/sample.bam \
#   --output-dir /path/to/output
###############################################################################

usage() {
    cat <<EOF
Usage:
  sbatch $0 --bam FILE --output-dir DIR [options]

Required:
  --bam FILE                 Input BAM file
  --output-dir DIR           Root directory for results, work files, and caches

Optional:
  --reference FASTA          Reference genome FASTA
  --sample-name NAME         Sample name (default: derived from BAM filename)
  --bam-min-coverage N       Minimum coverage (default: 5)
  --basecaller-config NAME   Basecaller configuration
  --bed FILE                 Restrict analysis to regions in a BED file
  --revision VERSION         Pin the Nextflow workflow revision
  --str                      Enable STR calling
  --no-snp                   Disable SNP calling
  --no-sv                    Disable SV calling
  --no-mod                   Disable modified-base analysis
  --no-phased                Disable phasing
  --no-resume                Start without Nextflow -resume
  -h, --help                 Show this help message
EOF
}

###############################################################################
# Input and settings
###############################################################################

BAM_INPUT=""
OUTPUT_DIR=""
REFERENCE_FASTA="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
SAMPLE_NAME=""
BAM_MIN_COVERAGE=5
BASECALLER_CONFIG="dna_r10.4.1_e8.2_400bps_hac"
REGION_BED=""
PIPELINE_REVISION=""
ENABLE_SNP=true
ENABLE_SV=true
ENABLE_STR=false
ENABLE_MOD=true
ENABLE_PHASED=true
RESUME=true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)
            [[ $# -ge 2 ]] || { echo "ERROR: --bam requires a value" >&2; exit 1; }
            BAM_INPUT="$2"
            shift 2
            ;;
        --output-dir)
            [[ $# -ge 2 ]] || { echo "ERROR: --output-dir requires a value" >&2; exit 1; }
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference)
            [[ $# -ge 2 ]] || { echo "ERROR: --reference requires a value" >&2; exit 1; }
            REFERENCE_FASTA="$2"
            shift 2
            ;;
        --sample-name)
            [[ $# -ge 2 ]] || { echo "ERROR: --sample-name requires a value" >&2; exit 1; }
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --bam-min-coverage)
            [[ $# -ge 2 ]] || { echo "ERROR: --bam-min-coverage requires a value" >&2; exit 1; }
            BAM_MIN_COVERAGE="$2"
            shift 2
            ;;
        --basecaller-config)
            [[ $# -ge 2 ]] || { echo "ERROR: --basecaller-config requires a value" >&2; exit 1; }
            BASECALLER_CONFIG="$2"
            shift 2
            ;;
        --bed)
            [[ $# -ge 2 ]] || { echo "ERROR: --bed requires a value" >&2; exit 1; }
            REGION_BED="$2"
            shift 2
            ;;
        --revision)
            [[ $# -ge 2 ]] || { echo "ERROR: --revision requires a value" >&2; exit 1; }
            PIPELINE_REVISION="$2"
            shift 2
            ;;
        --str) ENABLE_STR=true; shift ;;
        --no-snp) ENABLE_SNP=false; shift ;;
        --no-sv) ENABLE_SV=false; shift ;;
        --no-mod) ENABLE_MOD=false; shift ;;
        --no-phased) ENABLE_PHASED=false; shift ;;
        --no-resume) RESUME=false; shift ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$BAM_INPUT" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: --bam and --output-dir are required" >&2
    usage >&2
    exit 1
fi

if [[ ! "$BAM_MIN_COVERAGE" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --bam-min-coverage must be a non-negative number" >&2
    exit 1
fi

if [[ ! -f "$BAM_INPUT" ]]; then
    echo "ERROR: BAM file does not exist: $BAM_INPUT"
    exit 1
fi

BAM_PATH=$(realpath "$BAM_INPUT")

if [[ ! -r "$BAM_PATH" ]]; then
    echo "ERROR: BAM file is not readable: $BAM_PATH"
    exit 1
fi

case "${BAM_PATH,,}" in
    *.bam) ;;
    *)
        echo "ERROR: Input must be a BAM file: $BAM_PATH"
        exit 1
        ;;
esac

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA does not exist: $REFERENCE_FASTA" >&2
    exit 1
fi

REFERENCE_FASTA=$(realpath "$REFERENCE_FASTA")
SAMPLE_ROOT=$(realpath -m "$OUTPUT_DIR")

if [[ -n "$REGION_BED" ]]; then
    if [[ ! -f "$REGION_BED" ]]; then
        echo "ERROR: BED file does not exist: $REGION_BED" >&2
        exit 1
    fi
    REGION_BED=$(realpath "$REGION_BED")
fi

BAM_FILE=$(basename "$BAM_PATH")
if [[ -z "$SAMPLE_NAME" ]]; then
    SAMPLE_NAME="${BAM_FILE%.bam}"
    SAMPLE_NAME="${SAMPLE_NAME%.sorted.aligned}"
    SAMPLE_NAME="${SAMPLE_NAME%.aligned.sorted}"
    SAMPLE_NAME="${SAMPLE_NAME%.sorted}"
    SAMPLE_NAME="${SAMPLE_NAME%.aligned}"
fi

PIPELINE="elbayaliyev/wf-human-variation"

RESULTS_DIR="${SAMPLE_ROOT}/results"
WORK_DIR="${SAMPLE_ROOT}/work"
TMP_DIR="${SAMPLE_ROOT}/tmp"
NXF_HOME_DIR="${SAMPLE_ROOT}/.nextflow_home"
NXF_ASSETS_DIR="${SAMPLE_ROOT}/assets"
RUN_CONFIG="${SAMPLE_ROOT}/nextflow.config"

NEXTFLOW_SINGULARITY_CACHE="${SAMPLE_ROOT}/.cache/nextflow_singularity_images"
SINGULARITY_OCI_CACHE="${SAMPLE_ROOT}/.cache/singularity_oci_cache"

###############################################################################
# Modules and directories
###############################################################################

source /etc/profile.d/modules.sh 2>/dev/null || true
module purge
module load nextflow
module load singularity/3.7.4

command -v nextflow >/dev/null 2>&1 || {
    echo "ERROR: Nextflow is unavailable after module load."
    exit 1
}

command -v singularity >/dev/null 2>&1 || {
    echo "ERROR: Singularity is unavailable after module load."
    exit 1
}

mkdir -p \
    "$RESULTS_DIR" \
    "$WORK_DIR" \
    "$TMP_DIR" \
    "$NXF_HOME_DIR" \
    "$NXF_ASSETS_DIR" \
    "$NEXTFLOW_SINGULARITY_CACHE" \
    "$SINGULARITY_OCI_CACHE"

export TMPDIR="$TMP_DIR"
export NXF_HOME="$NXF_HOME_DIR"
export NXF_ASSETS="$NXF_ASSETS_DIR"
export NXF_WORK="$WORK_DIR"
export NXF_SINGULARITY_CACHEDIR="$NEXTFLOW_SINGULARITY_CACHE"
export SINGULARITY_CACHEDIR="$SINGULARITY_OCI_CACHE"

# Remove stale partial container downloads.
find "$NEXTFLOW_SINGULARITY_CACHE" \
    -maxdepth 1 \
    -type f \
    -name "*.pulling.*" \
    -mmin +360 \
    -delete 2>/dev/null || true

###############################################################################
# Nextflow configuration
###############################################################################

cat > "$RUN_CONFIG" <<CFG
conda.enabled = false

singularity {
    enabled = true
    autoMounts = true
    cacheDir = '${NEXTFLOW_SINGULARITY_CACHE}'
    pullTimeout = '2h'
}

process {
    executor = 'slurm'

    cpus = 1
    memory = 4.GB
    time = 2.h

    queue = 'acpu'
    clusterOptions = '--qos=cpu-normal'

    maxRetries = 1
    errorStrategy = {
        task.attempt <= 1 ? 'retry' : 'terminate'
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
    exitReadTimeout = '30 min'
    queueStatInterval = '30 sec'
    queueSize = 48
    perJobMemLimit = true
    submitRateLimit = '30 / 1 min'
    killBatchSize = 50
}
CFG

###############################################################################
# Run wf-human-variation
###############################################################################

echo "Starting wf-human-variation: ${SAMPLE_NAME}"
echo "BAM:     ${BAM_PATH}"
echo "Results: ${RESULTS_DIR}"

cd "$SAMPLE_ROOT"

NEXTFLOW_CMD=(nextflow run "$PIPELINE")
[[ -z "$PIPELINE_REVISION" ]] || NEXTFLOW_CMD+=(-r "$PIPELINE_REVISION")

NEXTFLOW_CMD+=(
    --bam "$BAM_PATH"
    --ref "$REFERENCE_FASTA"
    --sample_name "$SAMPLE_NAME"
    --bam_min_coverage "$BAM_MIN_COVERAGE"
    --override_basecaller_cfg "$BASECALLER_CONFIG"
    --out_dir "$RESULTS_DIR"
    -profile singularity
    -c "$RUN_CONFIG"
    -w "$WORK_DIR"
)

[[ "$ENABLE_SNP" == true ]] && NEXTFLOW_CMD+=(--snp)
[[ "$ENABLE_SV" == true ]] && NEXTFLOW_CMD+=(--sv)
[[ "$ENABLE_STR" == true ]] && NEXTFLOW_CMD+=(--str)
[[ "$ENABLE_MOD" == true ]] && NEXTFLOW_CMD+=(--mod)
[[ "$ENABLE_PHASED" == true ]] && NEXTFLOW_CMD+=(--phased)
[[ -z "$REGION_BED" ]] || NEXTFLOW_CMD+=(--bed "$REGION_BED")
[[ "$RESUME" == true ]] && NEXTFLOW_CMD+=(-resume)

"${NEXTFLOW_CMD[@]}"

echo "Completed wf-human-variation: ${SAMPLE_NAME}"
echo "Results: ${RESULTS_DIR}"
