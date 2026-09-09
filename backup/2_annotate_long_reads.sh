#!/usr/bin/bash
source /etc/profile.d/modules.sh

module load SVAFotate/v0.0.1
module load SVTK/29Apr2019
module load SURVIVOR/v1.0.7
module load ANNOTSV/3.3.8
module load bcftools/1.16
module load jdk/1.8.0_281
module load gatk/4.3.0.0



generate_dirs() {
	mkdir -p $work_dir/sv_inspector_report
	mkdir -p $reports_path/{logs,population_final,sample_final,sample_insertions,sample_translocations}
}


start_annotation_long_reads() {
	#ls $reports_path/sample_final/*.vcf > $reports_path/sample_final/sample.list;
	ls $work_dir/sniffles/filtered/*sniffles.vcf > $reports_path/sample_final/sample.list;

	./$software_path/tools/SURVIVOR merge $reports_path/sample_final/sample.list 500 1 1 0 0 0 $reports_path/population_final/$projectname.survivor.vcf;
	bcftools annotate --set-id '%CHROM\_%POS\_%INFO/END\_%INFO/SVTYPE' $reports_path/population_final/$projectname.survivor.vcf | bcftools annotate -x "INFO/SUPP_VEC" > $reports_path/population_final/$projectname.id.vcf;
	module load python/2.7.10;
	python2.7 /gpfs/software/genomics/svtyper/v0.7.0/env/bin/vcf_allele_freq.py $reports_path/population_final/$projectname.id.vcf > $reports_path/population_final/$projectname.AF.vcf;
	cat $software_path/headers/bed_header.txt > $reports_path/population_final/$projectname.bed;
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $reports_path/population_final/$projectname.AF.vcf >> $reports_path/population_final/$projectname.bed;
	grep "#CHROM" $reports_path/population_final/$projectname.AF.vcf | cut -f3,10- > $reports_path/population_final/$projectname.genotypes.vcf;
	bcftools query -f'%ID[\t%GT]\n' $reports_path/population_final/$projectname.AF.vcf >> $reports_path/population_final/$projectname.genotypes.vcf;
	sed -i 's/0\/1/0|1/g' $reports_path/population_final/$projectname.genotypes.vcf;sed -i 's/1\/1/1|1/g' $reports_path/population_final/$projectname.genotypes.vcf;
	cut -d$'\t' -f1-9 $reports_path/population_final/$projectname.AF.vcf > $reports_path/population_final/$projectname.sites.vcf
	export LC_ALL=C ;
	AnnotSV -SVinputFile $reports_path/population_final/$projectname.sites.vcf -SVminSize 0 -genomeBuild GRCh38 -overwrite 1 -bedtools /gpfs/software/genomics/bedtools/2.26/bin/bedtools -SVinputInfo 1 -outputFile $projectname -outputDir $reports_path/population_final/;
	java -jar -Xmx64G $software_path/htsSidra-1.2-jar-with-dependencies.jar processAnnotatedFileAnnotSV_3_8 $reports_path/population_final/$projectname.tsv;
	grep -v "SVTYPE=TRA" $reports_path/population_final/$projectname.sites.vcf > $reports_path/population_final/$projectname.sites.no_tra.vcf;
	svtk standardize --contigs $contigs_data --include-reference-sites $reports_path/population_final/$projectname.sites.no_tra.vcf $reports_path/population_final/$projectname.sites.no_tra.std.vcf manta;
	svafotate annotate --vcf $reports_path/population_final/$projectname.sites.no_tra.std.vcf --out $reports_path/population_final/$projectname.sites.no_tra.std.SVA_Fotate.vcf -f 0.8 -b $software_path/SVAFotate_core_SV_popAFs.GRCh38.bed.gz;
	svafotate annotate --vcf $reports_path/population_final/$projectname.sites.no_tra.std.SVA_Fotate.vcf --out $reports_path/population_final/$projectname.sites.no_tra.std.SVA_Fotate.topmed.vcf -f 0.8 -b $software_path/TOPMed.GRCh38.bed.gz;
	gatk SVAnnotate -V $reports_path/population_final/$projectname.sites.no_tra.std.SVA_Fotate.topmed.vcf --protein-coding-gtf $software_path/gencode_new_38.gtf -O $reports_path/population_final/$projectname.sites.no_tra.std.annotated.vcf;
	svtk vcf2bed --info ALL $reports_path/population_final/$projectname.sites.no_tra.std.annotated.vcf $reports_path/population_final/$projectname.sites.no_tra.std.annotated.bed;
	bcftools query -f "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%ID\n" $reports_path/population_final/$projectname.AF.vcf >> $reports_path/population_final/$projectname.bed;

}

echo "Starting SIDRA SV Annotation";


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo $DIR;

uuid=$(cat /dev/urandom | tr -cd 'a-f0-9' | head -c 16);

echo $uuid;

bam_list=$1
work_dir=$2
projectname=$3

ref_data="/gpfs/data_jrnas1/ref_data/Homo_sapiens/hs38DH/Sequences/WholeGenomeSequences/hs38DH.fa"
contigs_data="/gpfs/data_jrnas1/ref_data/Homo_sapiens/hs38DH/Sequences/WholeGenomeSequences/hs38DH.fa.fai"

#software_path=$DIR;
software_path="/gpfs/projects/tmedicine/KFakhroLAB/SV/software";

reports_path="$work_dir/reports"

#generate_dirs

start_annotation_long_reads