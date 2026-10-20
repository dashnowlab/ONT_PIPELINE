# ONT_PIPELINE

ONT pipelines for alignment, human-variation analysis, structural variants, and tandem repeats.

## Alignment and human variation

Submit alignment with a directory containing ONT FASTQ files:

```bash
sbatch ont_alignment.sh \
  --sample-dir /path/to/fastq_directory \
  --reference /path/to/reference.fasta \
  --output-dir /path/to/output \
  --threads 16
```

After alignment finishes, submit human-variation analysis with its BAM output:

```bash
sbatch ont_human_variation.sh \
  --bam /path/to/aligned.bam \
  --output-dir /path/to/output
```

The reference defaults to `human_GRCh38_no_alt_analysis_set.fasta` on the
Dashnow lab filesystem. Run `bash ont_human_variation.sh --help` for optional
reference, sample, calling, region, and workflow settings.

Run the primary SV/TR workflows from the repository root, for example:

Catalog builder:

mkdir -p /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed.gz.tmp/logs

sbatch -J annotsv_catalog -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh annotate_bed.nf \
  --input_bed /pl/active/dashnowlab/projects/1k_tr_paper/data/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed.gz \
  --software_path /pl/active/dashnowlab/software/ONT_PIPELINE \
  --work_dir /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.bed.gz.tmp \
  --batch_size 10000 \
  --annotsv_mem_gb 4 \
  --java_mem_gb 16 \
  --java_heap_gb 12 \
  -resume"


```bash
bash run_pipeline_slurm.sh long_read_sv.nf \
  --bam_list resources/samples/test_bam.list \
  --work_dir /path/to/output \
  --ref /path/to/reference.fasta
```

```bash
sbatch -J long_read_str -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str_v2.nf \
  --bam_list /pl/active/dashnowlab/software/ONT_PIPELINE/resources/samples/shaikh_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_STR_PATHO_PRECOMPUTED \
  --projectname SHAIKH_STR_PATHO_PRECOMPUTED \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software \
  --annotation_mode precomputed \
  --preannotated_tsv /pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.annotated.subset.tsv \
  -resume"
```

Repository-owned paths are resolved relative to this checkout. Python utilities
are in `scripts/`, bundled datasets in `resources/`, usage examples in `docs/`,
and superseded or experimental files in `backup/`.
