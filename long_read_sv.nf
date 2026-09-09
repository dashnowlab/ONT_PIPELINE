#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// Params & Help
// ─────────────────────────────────────────────────────────────────────────────
params.help              = params.help              ?: false
params.bam_list          = params.bam_list          ?: null
params.work_dir          = params.work_dir          ?: "${workDir}"
params.projectname       = params.projectname       ?: "sv_project"
params.ref               = params.ref               ?: "/gpfs/data_jrnas1/ref_data/Hsapiens/hg38/ONT/hg38_nohla.fa"
params.tr_bed            = params.tr_bed            ?: "/gpfs/data_jrnas1/ref_data/Hsapiens/hg38/ONT/human_GRCh38_no_alt_analysis_set.trf.bed"
params.software_path     = params.software_path     ?: "/pl/active/dashnowlab/work/ealiyev/SVTR_Analysis/software/"
params.caller            = params.caller            ?: "sniffles"   // sniffles | cutesv | severus

params.sample_ids_xlsx   = params.sample_ids_xlsx   ?: ""
params.acmg_gene_list    = params.acmg_gene_list    ?: "${params.software_path}/catalogs/ACMG_gene.list"
params.python_script     = params.python_script     ?: "${projectDir}/prioritize_sv_results.py"
params.python_libs       = params.python_libs       ?: "${params.software_path}/tools/python_libs"

if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║        DASHNOW_LAB Long-Read SV Discovery & Annotation Pipeline ║
    ╚══════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run long_read_sv.nf \\
            --bam_list          <bam_paths.txt>   \\
            --work_dir          <output_dir>      \\
            --projectname       <project_name>    \\
            --ref               <hg38.fa>         \\
            --tr_bed            <trf.bed>         \\
            --caller            sniffles          \\
            --software_path     <path_to_sw>      \\
            --sample_ids_xlsx   <optional.xlsx>   \\
            --acmg_gene_list    <optional.list>   \\
            --python_script     <prioritize.py>   \\
            --python_libs       <python_libs_dir>

    Callers: sniffles (default) | cutesv | severus
    """
    System.exit(0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Output sub-directories
// ─────────────────────────────────────────────────────────────────────────────
def snifflesRaw      = "${params.work_dir}/sniffles/raw"
def snifflesFiltered = "${params.work_dir}/sniffles/filtered"
def cuteSVRaw        = "${params.work_dir}/cutesv/raw"
def severusRaw       = "${params.work_dir}/severus/raw"
def popFinal         = "${params.work_dir}/reports/population_final"
def containerDir     = "${params.software_path}/containers"

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────
def sampleName = { path ->
    def leaf = (path instanceof java.nio.file.Path) ? path.fileName.toString() : path.toString()
    def bn = leaf.replaceFirst(/\.(bam|cram)$/, '')
    bn = bn.replaceFirst(/\.merged$/, '')
    return bn
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    Channel
        .fromPath(params.bam_list, checkIfExists: true)
        .splitText()
        .map  { it as String }
        .map  { it.replaceAll(/\r$/, '').trim() }
        .filter { line -> line && !line.startsWith('#') }
        .map { line ->
            def alignment = file(line, checkIfExists: true)
            def name      = sampleName(alignment)

            // Preserve the original list entry when deriving the index path.
            // For an HTTP Path, alignment.toString() can omit the URL scheme
            // and host (for example, returning only /1000g-ont/...).
            def sourceStr = line

            def indexCandidates

            if (sourceStr.endsWith('.cram')) {
                indexCandidates = [
                    file("${sourceStr}.crai"),
                    file(sourceStr.replaceFirst(/\.cram$/, '.crai'))
                ]
            }
            else if (sourceStr.endsWith('.bam')) {
                indexCandidates = [
                    file("${sourceStr}.bai"),
                    file(sourceStr.replaceFirst(/\.bam$/, '.bai'))
                ]
            }
            else {
                throw new IllegalArgumentException(
                    "Unsupported alignment format: ${alignment}. Expected BAM or CRAM."
                )
            }

            // Keep local HPC paths and remote HTTPS paths as Nextflow Path
            // objects. java.io.File can only check local files and therefore
            // incorrectly reports public remote indexes as missing.
            def idx = indexCandidates.find { candidate -> candidate.exists() }

            if (!idx) {
                throw new IllegalArgumentException(
                    "Missing index for ${alignment}. Checked: ${indexCandidates.join(', ')}"
                )
            }

            tuple(name, alignment, idx)
        }
        .set { bam_samples }

    ref     = file(params.ref,          checkIfExists: true)
    ref_fai = file("${params.ref}.fai", checkIfExists: true)

    if (params.caller == 'sniffles') {
        sniffles_raw = SNIFFLES_CALL(bam_samples, ref).raw_vcf
        called_vcfs  = SNIFFLES_FILTER(sniffles_raw).vcf
    }
    else if (params.caller == 'cutesv') {
        called_vcfs = CUTESV(bam_samples, ref).vcf
    }
    else if (params.caller == 'severus') {
        called_vcfs = SEVERUS(bam_samples, ref).vcf
    }
    else {
        error "Unknown caller '${params.caller}'. Use: sniffles | cutesv | severus"
    }

    all_vcfs       = called_vcfs.collect()

    survivor_vcf   = SURVIVOR_MERGE(all_vcfs)
    id_vcf         = BCFTOOLS_SET_ID(survivor_vcf)
    af_vcf         = ALLELE_FREQ(id_vcf)
    tables         = MAKE_TABLES(af_vcf).tables
    annotsv_out    = ANNOTATE_SV(tables).annotsv_out
    std_vcf        = SVTK_STANDARDIZE(annotsv_out)
    svafotate_out  = SVAFOTATE(std_vcf).svafotate_vcf
    annotated_vcf  = GATK_ANNOTATE(svafotate_out, ref_fai).annotated_vcf
    svtk_bed       = SVTK_VCF2BED(annotated_vcf).bed

    PRIORITIZE_SV_RESULTS(annotsv_out, tables, svtk_bed)
}

// ═════════════════════════════════════════════════════════════════════════════
// STAGE 1 — SV CALLERS
// ═════════════════════════════════════════════════════════════════════════════

process SNIFFLES_CALL {
    tag "${sample}"

    container "${containerDir}/sniffles_2.7.5_biocontainers.sif"

    cpus   8
    memory { 16.GB * task.attempt }
    time   { 6.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${snifflesRaw}", mode: 'copy', pattern: "*.raw.sniffles.vcf"

    input:
    tuple val(sample), path(bam), path(bai)
    path ref

    output:
    tuple val(sample), path("${sample}.raw.sniffles.vcf"), emit: raw_vcf

    script:
    def tr_bed = params.tr_bed
    """
    sniffles \\
        -i  ${bam} \\
        -v  ${sample}.raw.sniffles.vcf \\
        --tandem-repeats ${tr_bed} \\
        --reference      ${ref} \\
        --threads        ${task.cpus}
    """
}

process SNIFFLES_FILTER {
    tag "${sample}"

    container "${containerDir}/bcftools_biocontainers.sif"

    cpus   2
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${snifflesFiltered}", mode: 'copy', pattern: "*.merged.sniffles.vcf"

    input:
    tuple val(sample), path(raw_vcf)

    output:
    path "${sample}.merged.sniffles.vcf", emit: vcf

    script:
    """
    bcftools view \\
        -i '(SVLEN>=50 | SVLEN<=-50 | SVLEN=0 | SVLEN=1 | SVLEN=".")' \\
        ${raw_vcf} \\
    | grep -v -E 'SVTYPE=TRA|SVTYPE=BND|hs37d5|MT|chrUn|_KI2|_GL|_KB|_JH|chrM' \\
    | sed "s/SAMPLE/${sample}/" \\
    > ${sample}.merged.sniffles.vcf
    """
}

process CUTESV {
    tag "${sample}"

    container "${containerDir}/cutesv_2.1.3_biocontainers.sif"

    cpus   8
    memory { 16.GB * task.attempt }
    time   { 8.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${cuteSVRaw}", mode: 'copy', pattern: "*.cutesv.vcf"

    input:
    tuple val(sample), path(bam), path(bai)
    path ref

    output:
    path "${sample}.merged.cutesv.vcf", emit: vcf

    script:
    """
    mkdir -p tmp_${sample}

    cuteSV \\
        ${bam}  ${ref} \\
        ${sample}.merged.cutesv.vcf \\
        tmp_${sample} \\
        --threads              ${task.cpus} \\
        --min_support          10 \\
        --max_cluster_bias_INS 100 \\
        --diff_ratio_merging_INS 0.3 \\
        --max_cluster_bias_DEL 100 \\
        --diff_ratio_merging_DEL 0.3
    """
}

process SEVERUS {
    tag "${sample}"

    container "${containerDir}/severus_1.7_biocontainers.sif"

    cpus   16
    memory { 32.GB * task.attempt }
    time   { 12.h * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${severusRaw}", mode: 'copy', pattern: "*.severus.vcf"

    input:
    tuple val(sample), path(bam), path(bai)
    path ref

    output:
    path "${sample}.severus.vcf", emit: vcf

    script:
    def tr_bed = params.tr_bed
    """
    severus \\
        --target-bam ${bam} \\
        --out-dir    severus_${sample} \\
        -t           ${task.cpus} \\
        --vntr-bed   ${tr_bed}

    cp severus_${sample}/severus.vcf ${sample}.severus.vcf
    """
}

// ═════════════════════════════════════════════════════════════════════════════
// STAGE 2 — ANNOTATION CHAIN
// ═════════════════════════════════════════════════════════════════════════════

process SURVIVOR_MERGE {
    tag "survivor_merge"

    container "${containerDir}/survivor_1.0.7_biocontainers.sif"

    cpus   4
    memory { 16.GB * task.attempt }
    time   { 4.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path vcf_files

    output:
    path "${params.projectname}.survivor.vcf"

    script:
    def proj = params.projectname
    """
    ls *.vcf > sample.list
    SURVIVOR merge sample.list 500 1 1 0 0 0 ${proj}.survivor.vcf
    """
}

process BCFTOOLS_SET_ID {
    tag "bcftools_set_id"

    container "${containerDir}/bcftools_biocontainers.sif"

    cpus   2
    memory { 8.GB * task.attempt }
    time   { 1.h  * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path survivor_vcf

    output:
    path "${params.projectname}.id.vcf"

    script:
    def proj = params.projectname
    """
    bcftools annotate \\
        --set-id '%CHROM\\_%POS\\_%INFO/END\\_%INFO/SVTYPE' \\
        ${survivor_vcf} \\
    | bcftools annotate -x "INFO/SUPP_VEC" \\
    > ${proj}.id.vcf
    """
}

process ALLELE_FREQ {
    tag "allele_freq"

    cpus   2
    memory { 8.GB * task.attempt }
    time   { 1.h  * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path id_vcf

    output:
    path "${params.projectname}.AF.vcf"

    script:
    def proj = params.projectname
    def sw   = params.software_path
    """
    python3 ${sw}/vcf_allele_freq.py ${id_vcf} > ${proj}.AF.vcf
    """
}

process MAKE_TABLES {
    tag "make_tables"

    container "${containerDir}/bcftools_biocontainers.sif"

    cpus   2
    memory { 8.GB * task.attempt }
    time   { 1.h  * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path af_vcf

    output:
    tuple path("${params.projectname}.AF.vcf"),
          path("${params.projectname}.sites.vcf"),
          path("${params.projectname}.bed"),
          path("${params.projectname}.genotypes.vcf"), emit: tables

    script:
    def proj = params.projectname
    """
    echo -e "CHROM\\tSTART\\tEND\\tSVTYPE\\tSVLEN\\tID" > ${proj}.bed
    bcftools query \\
        -f "%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%ID\\n" \\
        ${af_vcf} >> ${proj}.bed

    grep "#CHROM" ${af_vcf} | cut -f3,10- > ${proj}.genotypes.vcf
    bcftools query -f'%ID[\\t%GT]\\n' ${af_vcf} >> ${proj}.genotypes.vcf
    sed -i 's/0\\/1/0|1/g; s/1\\/1/1|1/g' ${proj}.genotypes.vcf

    cut -f1-9 ${af_vcf} > ${proj}.sites.vcf
    cp ${af_vcf} ${proj}.AF.vcf
    """
}

process ANNOTATE_SV {
    tag "annotsv_htsSidra"

    container "${containerDir}/annotsv_3.5.8_biocontainers.sif"

    cpus   8
    memory { 32.GB * task.attempt }
    time   { 8.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    tuple path(af_vcf), path(sites_vcf), path(bed), path(genotypes_vcf)

    output:
    tuple path("${params.projectname}.tsv.processed.tsv"),
          path("${params.projectname}.sites.no_tra.vcf"), emit: annotsv_out

    script:
    def proj = params.projectname
    def sw   = params.software_path
    """
    AnnotSV \\
        -SVinputFile  ${sites_vcf} \\
        -SVminSize    0 \\
        -genomeBuild  GRCh38 \\
        -overwrite    1 \\
        -bedtools     \$(which bedtools) \\
        -SVinputInfo  1 \\
        -annotationsDir ${sw}/tools/AnnotSV/share/AnnotSV/ \\
        -outputFile   ${proj} \\
        -outputDir    .

    java -jar -Xmx64G ${sw}/tools/htsSidra-1.2-jar-with-dependencies.jar \\
        processAnnotatedFileAnnotSV_3_8 ${proj}.tsv

    grep -v "SVTYPE=TRA" ${sites_vcf} > ${proj}.sites.no_tra.vcf
    """
}

process SVTK_STANDARDIZE {
    tag "svtk_standardize"

    container "${containerDir}/svtk_0.0.20190615_2_biocontainers.sif"

    cpus   4
    memory { 16.GB * task.attempt }
    time   { 2.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    tuple path(annotsv_tsv), path(no_tra_vcf)

    output:
    path "${params.projectname}.sites.no_tra.std.vcf"

    script:
    def proj    = params.projectname
    def contigs = "${params.ref}.fai"
    """
    svtk standardize \\
        --contigs                 ${contigs} \\
        --include-reference-sites \\
        ${no_tra_vcf} \\
        ${proj}.sites.no_tra.std.vcf \\
        manta
    """
}

process SVAFOTATE {
    tag "svafotate"

    conda '/projects/ealiyev@xsede.org/software/anaconda/envs/svafotate-env'

    cpus   4
    memory { 16.GB * task.attempt }
    time   { 4.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path std_vcf

    output:
    path "${params.projectname}.sites.no_tra.std.SVA_Fotate.topmed.vcf", emit: svafotate_vcf

    script:
    def proj = params.projectname
    def sw   = params.software_path
    """
    svafotate annotate \\
        --vcf  ${std_vcf} \\
        --out  ${proj}.sites.no_tra.std.SVA_Fotate.vcf \\
        -f 0.8 \\
        -b ${sw}/SVAFotate_core_SV_popAFs.GRCh38.bed.gz

    svafotate annotate \\
        --vcf  ${proj}.sites.no_tra.std.SVA_Fotate.vcf \\
        --out  ${proj}.sites.no_tra.std.SVA_Fotate.topmed.vcf \\
        -f 0.8 \\
        -b ${sw}/TOPMed.GRCh38.bed.gz
    """
}

process GATK_ANNOTATE {
    tag "gatk_svAnnotate"

    container "${containerDir}/gatk4_4.6.2.0_biocontainers.sif"

    cpus   4
    memory { 16.GB * task.attempt }
    time   { 6.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path topmed_vcf
    path ref_fai

    output:
    path "${params.projectname}.sites.no_tra.std.annotated.vcf", emit: annotated_vcf

    script:
    def proj = params.projectname
    def sw   = params.software_path
    """
    gatk SVAnnotate \\
        -V ${topmed_vcf} \\
        --protein-coding-gtf ${sw}/gencode_new_38.gtf \\
        -O ${proj}.sites.no_tra.std.annotated.vcf
    """
}

process SVTK_VCF2BED {
    tag "svtk_vcf2bed"

    container "${containerDir}/svtk_0.0.20190615_biocontainers.sif"

    cpus   4
    memory { 16.GB * task.attempt }
    time   { 2.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    path annotated_vcf

    output:
    path "${params.projectname}.sites.no_tra.std.annotated.bed", emit: bed

    script:
    def proj = params.projectname
    """
    svtk vcf2bed --info ALL \\
        ${annotated_vcf} \\
        ${proj}.sites.no_tra.std.annotated.bed
    """
}

process PRIORITIZE_SV_RESULTS {
    tag "prioritize_sv_results_filtering"

    cpus   2
    memory { 8.GB * task.attempt }
    time   { 4.h  * task.attempt  }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${popFinal}", mode: 'copy'

    input:
    tuple path(processed_tsv), path(no_tra_vcf)
    tuple path(af_vcf), path(sites_vcf), path(bed), path(genotypes_vcf)
    path svtk_bed

    output:
    path "${params.projectname}.SV.filtered.tsv"
    path "${params.projectname}.rare_acmg_events.tsv", optional: true
    path "${params.projectname}.SV_variants_and_results.xlsx", optional: true

    script:
    def proj        = params.projectname
    def sampleSheet = params.sample_ids_xlsx ?: ""
    def acmgList    = params.acmg_gene_list ?: ""
    def pyScript    = params.python_script
    def pyLibs      = params.python_libs

    """
    set -euo pipefail

    PYTHONPATH_PREV=\${PYTHONPATH:-}
    export PYTHONPATH=${pyLibs}:\${PYTHONPATH_PREV}

    python3 ${pyScript} \\
        ${processed_tsv} \\
        "${sampleSheet}" \\
        "${acmgList}" \\
        ${svtk_bed} \\
        ${genotypes_vcf} \\
        ${proj}.SV.filtered.tsv \\
        ${proj}.rare_acmg_events.tsv \\
        ${proj}.SV_variants_and_results.xlsx
    """
}
