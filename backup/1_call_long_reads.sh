#!/usr/bin/bash
source /etc/profile.d/modules.sh
module load bcftools/1.9
module load miniforge/24.11.3-0

start_calling_sniffles() {
  set -euo pipefail

  # Ensure required dirs exist
  mkdir -p "$work_dir/sniffles/raw" "$work_dir/sniffles/filtered" "$work_dir/reports/logs"

  # $1 should be the input list file
  while IFS= read -r row || [[ -n "$row" ]]; do
    [[ -z "$row" || "$row" =~ ^[[:space:]]*# ]] && continue

    cleanline=$(readlink -f -- "$row")
    line="${cleanline%%,*}"                          # take field before comma
    echo "Text read from file: $line"

    fname=$(basename -- "$line")
    echo "$fname"

    name="${fname%%.*}"
    echo "$name"

    bam_path=$(echo "$line" | cut -d'/' -f1-11)
    echo "$bam_path"

    file_to_check="$work_dir/sniffles/raw/${name}.merged.sniffles.vcf"
    echo "$file_to_check"

    # Submit only if missing or <= 100 KB
    if [[ ! -f "$file_to_check" || $(stat -c%s "$file_to_check") -le 102400 ]]; then
      echo "File missing/small. Submitting job: $file_to_check"

      # IMPORTANT: initialize conda and activate env INSIDE the job
      job_cmd=$'set -euo pipefail
				source "$(conda info --base)/etc/profile.d/conda.sh"
				conda activate sniffles.2.6.3

				sniffles \
				  -i '"$line"$' \
				  -v '"$work_dir/sniffles/raw/${name}.merged.sniffles.vcf"$' \
				  --tandem-repeats '"$tr_bed"$' \
				  --reference '"$ref_data"$'

				bcftools view -i '"'(SVLEN>=50 | SVLEN<=-50 | SVLEN = 0 | SVLEN = 1 | SVLEN = ".")'"' \
				  '"$work_dir/sniffles/raw/${name}.merged.sniffles.vcf"$' \
				| grep -v -E '"SVTYPE=TRA|SVTYPE=BND|hs37d5|MT|chrUn|_KI2|_GL|_KB|_JH|chrM"' \
				| sed '"s/SAMPLE/'"$name"'/g"' \
				> '"$work_dir/sniffles/filtered/${name}.merged.sniffles.vcf"$''

      # Slurm submission: -p amilan, -q normal, 2 cores, 2G RAM
      sbatch \
        --job-name="sniffles_${name}" \
        --nodes=1 \
        --cpus-per-task=2 \
        --mem=2G \
        --partition=amilan \
        --qos=normal \
        --output="$work_dir/reports/logs/${name}-sniffles-out.txt" \
        --error="$work_dir/reports/logs/${name}-sniffles-err.txt" \
        --export=ALL \
        --wrap "$job_cmd"

    else
      echo "Output $file_to_check exists (>100 KB). Skipping..."
    fi
  done < "$list_file"
}

echo "Starting Long-Read SV Discovery";

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR;

uuid=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16);

echo $uuid;

bam_list=$1
work_dir=$2
projectname=$3

ref_data="/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
tr_bed="/pl/active/dashnowlab/work/ealiyev/software/blacklist_regions/human_GRCh38_no_alt_analysis_set.trf.bed"
software_path="/pl/active/dashnowlab/work/ealiyev/software";

reports_path="$work_dir/reports"

start_calling_sniffles $bam_list

