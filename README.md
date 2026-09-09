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

```bash
bash run_pipeline_slurm.sh long_read_sv.nf \
  --bam_list resources/samples/test_bam.list \
  --work_dir /path/to/output \
  --ref /path/to/reference.fasta
```

Repository-owned paths are resolved relative to this checkout. Python utilities
are in `scripts/`, bundled datasets in `resources/`, usage examples in `docs/`,
and superseded or experimental files in `backup/`.
