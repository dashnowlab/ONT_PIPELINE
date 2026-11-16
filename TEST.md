# Pipeline test commands

Run from the cluster. First create the test directories:

```bash
mkdir -p /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/{data,alignment,human_variation,str,sv,catalog,logs}
cd /pl/active/dashnowlab/software/ONT_PIPELINE
```

## Download human APOB FASTQ test data

This small ONT human FASTQ comes from the EPI2ME Cas9/APOB tutorial. It is
already gzip-compressed and can be passed directly to the alignment pipeline.

```bash
mkdir -p /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/data/apob_fastq
cd /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/data
wget -c -O apob_fastq/apob_reads.fastq.gz \
  'https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/human_amp_snp_calling/apob_reads.fastq.gz'
wget -c -O apob_region.bed \
  'https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/human_amp_snp_calling/apob_region.bed'
cd /pl/active/dashnowlab/software/ONT_PIPELINE
```

## Alignment (`ont_alignment.sh`)

```bash
sbatch \
  --output=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/ont_align_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/ont_align_%j.err \
  ont_alignment.sh \
  --sample-dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/data/apob_fastq \
  --reference /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --output-dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/alignment \
  --threads 16
```

## Human variation (`ont_human_variation.sh`)

Run this after alignment completes. The older FASTQ does not provide a reliable
basecaller-model identifier, so this smoke test runs structural-variant calling
only. The zero coverage threshold bypasses the whole-genome minimum-coverage
gate for this small targeted dataset.

```bash
sbatch \
  --output=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/ont_variation_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/ont_variation_%j.err \
  ont_human_variation.sh \
  --bam /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/alignment/results/apob_fastq.sorted.aligned.bam \
  --reference /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --bed /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/data/apob_region.bed \
  --sample-name TEST_APOB \
  --output-dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/human_variation \
  --bam-min-coverage 0 \
  --no-snp \
  --no-mod \
  --no-phased \
  --no-resume
```

## STR (`long_read_str_v2.nf`)

```bash
sbatch -J test_str -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str_v2.nf \
  --bam_list /pl/active/dashnowlab/software/ONT_PIPELINE/resources/samples/test_bam_str.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/str \
  --projectname TEST_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --software_dir /pl/active/dashnowlab/software/ONT_PIPELINE \
  --atarva_bed /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/STRchive-disease-loci.hg38.atarva.bed.gz \
  --atarva_json /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TR_catalog.5599658_loci.20261016_141520.json.gz"
```

## SV (`long_read_sv.nf`, Sniffles)

```bash
sbatch -J test_sv -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_sv.nf \
  --bam_list /pl/active/dashnowlab/software/ONT_PIPELINE/resources/samples/test_bam_sv \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/sv \
  --projectname TEST_SV \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --tr_bed /pl/active/dashnowlab/software/ONT_PIPELINE/blacklist_regions/human_GRCh38_no_alt_analysis_set.trf.bed \
  --software_path /pl/active/dashnowlab/software/ONT_PIPELINE \
  --caller sniffles"
```

## Catalog annotation (`annotate_bed.nf`)

```bash
sbatch -J test_catalog -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.out \
  --error=/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh annotate_bed.nf \
  --input_bed /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/STRchive-disease-loci.hg38.atarva.bed.gz \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/catalog/STRchive-disease-loci.hg38.atarva.bed.gz.tmp \
  --software_path /pl/active/dashnowlab/software/ONT_PIPELINE \
  --batch_size 10000 \
  --annotsv_mem_gb 4 \
  --java_mem_gb 16 \
  --java_heap_gb 12"
```
