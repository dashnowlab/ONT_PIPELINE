#!/bin/bash -l
#SBATCH --job-name=ont_variation
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/alpine/ealiyev@xsede.org/data/slurm_logs/%x_%j.out
#SBATCH --error=/scratch/alpine/ealiyev@xsede.org/data/slurm_logs/%x_%j.err

set -euo pipefail

###############################################################################
# Usage
#
# sbatch ont_human_variation.sh \
#   /path/to/sample.bam
###############################################################################

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 /path/to/sample.bam"
    exit 1
fi

###############################################################################
# Input and settings
###############################################################################

BAM_INPUT="$1"

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

BAM_FILE=$(basename "$BAM_PATH")
SAMPLE_NAME="${BAM_FILE%.bam}"

# Remove common alignment suffixes from the sample name.
SAMPLE_NAME="${SAMPLE_NAME%.sorted.aligned}"
SAMPLE_NAME="${SAMPLE_NAME%.aligned.sorted}"
SAMPLE_NAME="${SAMPLE_NAME%.sorted}"
SAMPLE_NAME="${SAMPLE_NAME%.aligned}"

SCRATCH_ROOT="/scratch/alpine/ealiyev@xsede.org/data"
REFERENCE_FASTA="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
PIPELINE="elbayaliyev/wf-human-variation"

SAMPLE_ROOT="${SCRATCH_ROOT}/${SAMPLE_NAME}/human_variation"
RESULTS_DIR="${SAMPLE_ROOT}/results"
WORK_DIR="${SAMPLE_ROOT}/work"
TMP_DIR="${SAMPLE_ROOT}/tmp"
NXF_HOME_DIR="${SAMPLE_ROOT}/.nextflow_home"
NXF_ASSETS_DIR="${SAMPLE_ROOT}/assets"
RUN_CONFIG="${SAMPLE_ROOT}/nextflow.config"

NEXTFLOW_SINGULARITY_CACHE="${SCRATCH_ROOT}/nextflow_singularity_images"
SINGULARITY_OCI_CACHE="${SCRATCH_ROOT}/singularity_oci_cache"

if [[ ! -f "$REFERENCE_FASTA" ]]; then
    echo "ERROR: Reference FASTA does not exist: $REFERENCE_FASTA"
    exit 1
fi

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

nextflow run "$PIPELINE" \
    --bam "$BAM_PATH" \
    --ref "$REFERENCE_FASTA" \
    --sample_name "$SAMPLE_NAME" \
    --bam_min_coverage 5 \
    --phased \
    --snp \
    --sv \
    --mod \
    --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_hac \
    --out_dir "$RESULTS_DIR" \
    -profile singularity \
    -c "$RUN_CONFIG" \
    -w "$WORK_DIR" \
    -resume

echo "Completed wf-human-variation: ${SAMPLE_NAME}"
echo "Results: ${RESULTS_DIR}"
