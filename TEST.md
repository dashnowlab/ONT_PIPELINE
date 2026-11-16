# Pipeline test commands

Run from the cluster. First create the test directories:

```bash
mkdir -p /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/TEST/{str,sv,catalog,logs}
cd /pl/active/dashnowlab/software/ONT_PIPELINE
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
