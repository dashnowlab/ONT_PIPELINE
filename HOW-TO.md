sbatch -J long_read_sv -p amilan --qos=normal --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_pipeline_slurm.sh long_read_sv.nf --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_SV/shaikh_bam.list --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_SV --projectname SHAIKH_SV --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta --tr_bed /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/blacklist_regions/human_GRCh38_no_alt_analysis_set.trf.bed --software_path /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software --caller sniffles -resume"

sbatch -J long_read_sv -p amilan --qos=normal --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_pipeline_slurm.sh long_read_sv.nf --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_SV/shaikh_bam.list --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_SV --projectname SHAIKH_SV --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta --tr_bed /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/blacklist_regions/human_GRCh38_no_alt_analysis_set.trf.bed --software_path /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software --caller sniffles"

Caller  Coordinates used for BED  Final normalized ID
LongTR  INFO/START + INFO/END VCF ID
Medaka  VCF ID like chr1_57367043_57367118  VCF ID
STRdust POS + INFO/END  generated chr_POS_END
ATaRVa  INFO/START + INFO/END INFO/ID

sbatch -J long_read_str -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str.nf \
  --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/deveson_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/DEVESON_STR \
  --projectname SHAIKH_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --caller longtr \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software"

  sbatch -J long_read_str -p amilan --qos=normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str.nf \
  --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/shaikh_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/SHAIKH_STR \
  --projectname SHAIKH_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --caller atarva \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software"

  sbatch -J long_read_str -p amilan --qos=normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str.nf \
  --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/deveson_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/DEVESON_STR \
  --projectname SHAIKH_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --caller atarva \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software \
  -resume"

    sbatch -J long_read_str -p amilan --qos=normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str.nf \
  --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/deveson_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/DEVESON_STR \
  --projectname SHAIKH_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --caller medaka \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software"

      sbatch -J long_read_str -p amilan --qos=normal --time=23:00:00 --mem=2G \
  --output=logs/%x_%j.out --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh long_read_str.nf \
  --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/deveson_bam.list \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/DEVESON_STR \
  --projectname SHAIKH_STR \
  --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta \
  --caller strdust \
  --software_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software"

  HPRC_genomes

  sbatch -J long_read_sv -p acpu --qos=cpu-normal --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_pipeline_slurm.sh long_read_sv.nf --bam_list /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/HPRC_SV/cram.list --work_dir /scratch/alpine/ealiyev@xsede.org/HPRC_SV --projectname HPRC_SV --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta --tr_bed /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/blacklist_regions/human_GRCh38_no_alt_analysis_set.trf.bed --software_path /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software --caller sniffles -resume"

  mkdir -p logs

sbatch \
  -J annotsv \
  -p acpu \
  --qos=cpu-normal \
  --time=23:00:00 \
  --mem=2G \
  --output=logs/%x_%j.out \
  --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh annotsv_from_tsv.nf \
  --input_tsv /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/tr_from_sv/HG002_03_04_multisample.trf.novelyFilter.tsv \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/tr_from_sv/"

  sbatch \
  -J annotsv \
  -p acpu \
  --qos=cpu-normal \
  --time=23:00:00 \
  --mem=2G \
  --output=logs/%x_%j.out \
  --error=logs/%x_%j.err \
  --wrap="bash run_pipeline_slurm.sh annotsv_from_tsv.nf \
  --input_tsv /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/tr_from_sv/hprc_multisample.trf.noveltyFiltered.tsv \
  --work_dir /pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/tr_from_sv/"

