#!/usr/bin/env bash
module load singularity

set -u
set -o pipefail

OUTDIR="/pl/active/dashnowlab/work/ealiyev/SV/software/containers"
mkdir -p "${OUTDIR}"

IMAGES=(
  "sniffles_2.7.5_biocontainers.sif|docker://quay.io/biocontainers/sniffles:2.7.5--pyhdfd78af_0"
  "bcftools_biocontainers.sif|docker://quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0"
  "cutesv_2.1.3_biocontainers.sif|docker://quay.io/biocontainers/cutesv:2.1.3--pyhdfd78af_0"
  "severus_1.7_biocontainers.sif|docker://quay.io/biocontainers/severus:1.7--pyhdfd78af_0"
  "survivor_1.0.7_biocontainers.sif|docker://quay.io/biocontainers/survivor:1.0.7--h077b44d_7"
  "gatk4_4.6.2.0_biocontainers.sif|docker://quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
  "svtk_0.0.20190615_2_biocontainers.sif|docker://quay.io/biocontainers/svtk:0.0.20190615--py39h38f01e4_2"
  "medaka:2.2.1--py310h237e959_0.sif|docker://quay.io/biocontainers/medaka:2.2.1--py310h237e959_0"
)

echo "=============================================================="
echo "Using singularity"
echo "Output dir: ${OUTDIR}"
echo "=============================================================="
echo "Images to process:"
for entry in "${IMAGES[@]}"; do
    IFS='|' read -r sif image <<< "${entry}"
    echo "  ${image}"
    echo "    -> ${OUTDIR}/${sif}"
done
echo

cd "${OUTDIR}" || exit 1

FAILED=()
PULLED=()
SKIPPED=()

for entry in "${IMAGES[@]}"; do
    IFS='|' read -r sif image <<< "${entry}"

    echo "--------------------------------------------------------------"
    echo "Image : ${image}"
    echo "Target: ${OUTDIR}/${sif}"
    echo "--------------------------------------------------------------"

    if [[ -s "${sif}" ]]; then
        echo "SKIP: ${sif} already exists"
        SKIPPED+=("${sif}")
        echo
        continue
    fi

    tmp_sif="${sif%.sif}.tmp.sif"
    rm -f "${tmp_sif}"

    if singularity pull "${tmp_sif}" "${image}"; then
        mv -f "${tmp_sif}" "${sif}"
        echo "SUCCESS: ${image}"
        PULLED+=("${sif}")
    else
        rm -f "${tmp_sif}"
        echo "FAILED : ${image}" >&2
        FAILED+=("${image}")
    fi

    echo
done

echo "=============================================================="
echo "Summary"
echo "=============================================================="
echo "Pulled : ${#PULLED[@]}"
for x in "${PULLED[@]}"; do
    echo "  ${x}"
done

echo "Skipped: ${#SKIPPED[@]}"
for x in "${SKIPPED[@]}"; do
    echo "  ${x}"
done

echo "Failed : ${#FAILED[@]}"
for x in "${FAILED[@]}"; do
    echo "  ${x}"
done

if [[ "${#FAILED[@]}" -gt 0 ]]; then
    exit 1
fi