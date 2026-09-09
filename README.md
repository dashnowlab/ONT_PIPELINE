# ONT_PIPELINE

ONT pipelines for alignment, human-variation analysis, structural variants, and tandem repeats.

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
