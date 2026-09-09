#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// Params & Help
// ─────────────────────────────────────────────────────────────────────────────
params.help              = params.containsKey('help')              ? params.help              : false
params.bam_list          = params.containsKey('bam_list')          ? params.bam_list          : null
params.work_dir          = params.containsKey('work_dir')          ? params.work_dir          : "${workDir}"
params.projectname       = params.containsKey('projectname')       ? params.projectname       : "str_project"
params.ref               = params.containsKey('ref')               ? params.ref               : "/gpfs/data_jrnas1/ref_data/Hsapiens/hg38/ONT/hg38_nohla.fa"
params.base_dir          = params.containsKey('base_dir')          ? params.base_dir          : null
params.caller            = params.containsKey('caller')            ? params.caller            : "atarva"   // atarva | longtr | longtr_chrom | strdust | medaka

// Software root
params.software_dir      = params.containsKey('software_dir')      ? params.software_dir      : projectDir

// Catalogs / resources
//params.atarva_bed        = params.atarva_bed        ?: "${params.software_dir}/catalogs/STRchive-disease-loci.hg38.atarva.bed.gz"
params.atarva_bed        = params.containsKey('atarva_bed')        ? params.atarva_bed        : "/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/tr_explorer_catalog/TR_catalog_for_atarva.bed.gz"
params.longtr_bed        = params.containsKey('longtr_bed')        ? params.longtr_bed        : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_192919.LongTR.bed"
params.longtr_json       = params.containsKey('longtr_json')       ? params.longtr_json       : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_193020.json.gz"
params.strdust_bed       = params.containsKey('strdust_bed')       ? params.strdust_bed       : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260428_202657.STRdust.bed"
params.medaka_bed        = params.containsKey('medaka_bed')        ? params.medaka_bed        : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260428_202124.Medaka.bed"

// Caller-specific JSON catalogs for stranno annotation
// Update these paths if your JSON catalog files have different names.
params.atarva_json       = params.containsKey('atarva_json')       ? params.atarva_json       : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_193020.json.gz"
params.strdust_json      = params.containsKey('strdust_json')      ? params.strdust_json      : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_193020.json.gz"
params.medaka_json       = params.containsKey('medaka_json')       ? params.medaka_json       : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_193020.json.gz"

// Software / envs
params.atarva_bin        = params.containsKey('atarva_bin')        ? params.atarva_bin        : "/projects/ealiyev@xsede.org/software/anaconda/envs/atarva_0.7.1/bin/atarva"
params.strdust_bin       = params.containsKey('strdust_bin')       ? params.strdust_bin       : "/pl/active/dashnowlab/software/STRdust-0.16.0/STRdust-linux-musl"
params.medaka_conda      = params.containsKey('medaka_conda')      ? params.medaka_conda      : "envs/medaka-2.1.1.yaml"
params.atarva_conda      = params.containsKey('atarva_conda')      ? params.atarva_conda      : "${params.software_dir}/tools/ATaRVa/environment.yaml"

// Annotation toolchain
params.stranno_bin       = params.containsKey('stranno_bin')       ? params.stranno_bin       : "${params.software_dir}/tools/stranno-linux-x86_64-static_0.3.0"
params.annotsv_bin       = params.containsKey('annotsv_bin')       ? params.annotsv_bin       : "${params.software_dir}/tools/AnnotSV/bin/AnnotSV"
params.htssidra_jar      = params.containsKey('htssidra_jar')      ? params.htssidra_jar      : "${params.software_dir}/tools/htsSidra-1.2-jar-with-dependencies.jar"
params.z_score_script    = params.containsKey('z_score_script')    ? params.z_score_script    : "${params.software_dir}/scripts/z_score_regularized_final.py"
params.vcf_to_bed_script = params.containsKey('vcf_to_bed_script') ? params.vcf_to_bed_script : "${params.software_dir}/scripts/vcf_to_bed.py"
params.vcf_to_tsv_script = params.containsKey('vcf_to_tsv_script') ? params.vcf_to_tsv_script : "${params.software_dir}/scripts/vcf_to_tsv.py"
params.left_join_script  = params.containsKey('left_join_script')  ? params.left_join_script  : "${params.software_dir}/scripts/left_join.py"

if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║                 DASHNOW_LAB STR Calling Pipeline                ║
    ╚══════════════════════════════════════════════════════════════════╝

    Usage:
        nextflow run long_read_str.nf \\
            --bam_list      <bam_or_cram_and_karyotype.tsv> \\
            --work_dir      <output_dir>                    \\
            --projectname   <project_name>                  \\
            --ref           <reference.fa>                  \\
            --caller        longtr                          \\
            --software_dir  <software_root>

    Input format for --bam_list:
        TAB-separated with 2 columns:
            <bam|cram_path>    <karyotype>

        Example:
            /path/sample1.cram    XY
            /path/sample2.bam     XX

    Callers:
        atarva | longtr | longtr_chrom | strdust | medaka

    Note:
        For atarva, longtr, longtr_chrom, strdust and medaka, the pipeline also runs:
        stranno -> z_score_regularized_final.py -> vcf_to_bed.py -> AnnotSV -> htsSidra -> vcf_to_tsv.py -> left_join.py
    """
    System.exit(0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Output sub-directories
// ─────────────────────────────────────────────────────────────────────────────
def atarvaOut      = "${params.work_dir}/atarva"
def longtrChromOut = "${params.work_dir}/longtr_chrom"
def strdustOut     = "${params.work_dir}/strdust"
def longtrOut      = "${params.work_dir}/longtr"
def medakaOut      = "${params.work_dir}/medaka"
def containerDir     = "${params.software_dir}/containers"

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────
def resolvePath = { String line ->
    def p = line.trim()
    if (!p) return null
    if (p.startsWith('#')) return null
    if (params.base_dir && !p.startsWith('/')) {
        return file("${params.base_dir}/${p}")
    }
    return file(p)
}

def sampleName = { x ->
    def leaf = (x instanceof java.nio.file.Path) ? x.fileName.toString()
             : (x instanceof java.io.File)       ? x.name
             :                                     x.toString()

    def bn = leaf.replaceFirst(/\.(cram|bam)$/, '')
    bn = bn.replaceFirst(/\.merged$/, '')
    return bn
}

def alignmentIndex = { alignment ->
    def pathStr = alignment.toString()
    def candidates

    if (pathStr.endsWith('.cram')) {
        candidates = [
            file("${pathStr}.crai"),
            file(pathStr.replaceFirst(/\.cram$/, '.crai'))
        ]
    }
    else if (pathStr.endsWith('.bam')) {
        candidates = [
            file("${pathStr}.bai"),
            file(pathStr.replaceFirst(/\.bam$/, '.bai'))
        ]
    }
    else {
        throw new IllegalArgumentException(
            "Unsupported alignment format: ${alignment}. Expected BAM or CRAM."
        )
    }

    def idx = candidates.find { candidate ->
        new File(candidate.toString()).exists()
    }

    if (!idx) {
        throw new IllegalArgumentException(
            "Missing index for ${alignment}. Checked: ${candidates.join(', ')}"
        )
    }

    return idx
}

// ─────────────────────────────────────────────────────────────────────────────
// Workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    def called_vcfs
    def final_reports

    Channel
        .fromPath(params.bam_list, checkIfExists: true)
        .splitText()
        .map { it as String }
        .map { it.replaceAll(/\r$/, '') }
        .filter { line -> line.trim() && !line.trim().startsWith('#') }
        .map { line ->
            def parts = line.split('\t', -1)
            if (parts.size() < 2) {
                throw new IllegalArgumentException(
                    "Expected two TAB-separated columns: <bam|cram>\\t<karyotype>; got: '${line}'"
                )
            }

            def alnPath = parts[0].trim()
            def karyo   = parts[1].trim()

            def aln = resolvePath(alnPath)
            if (aln == null) {
                throw new IllegalArgumentException("Could not resolve path from line: '${line}'")
            }

            def sample = sampleName(aln)

            if (!(karyo in ['XX', 'XY'])) {
                throw new IllegalArgumentException(
                    "Unsupported karyotype '${karyo}' for ${aln}. Expected XX or XY."
                )
            }

            def idx = alignmentIndex(aln)

            tuple(sample, aln, idx, karyo)
        }
        .set { aligned_samples }   // [sample, aln, idx, karyotype]

    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)

    if (params.caller == 'atarva') {
        called_vcfs = ATARVA(aligned_samples, ref, fai)
        final_reports = ANNOTATE(called_vcfs)
        final_reports.view { "Finished ATaRVa + annotation: ${it}" }
    }
    else if (params.caller == 'strdust') {
        called_vcfs = STRDUST(aligned_samples, ref, fai)
        final_reports = ANNOTATE(called_vcfs)
        final_reports.view { "Finished STRdust + annotation: ${it}" }
    }
    else if (params.caller == 'medaka') {
        called_vcfs = MEDAKA(aligned_samples, ref, fai)
        final_reports = ANNOTATE(called_vcfs)
        final_reports.view { "Finished Medaka Tandem + annotation: ${it}" }
    }
    else if (params.caller == 'longtr') {
        called_vcfs = LONGTR(aligned_samples, ref, fai)
        final_reports = ANNOTATE(called_vcfs)
        final_reports.view { "Finished LongTR + annotation: ${it}" }
    }
    else if (params.caller == 'longtr_chrom') {

        ch_chroms = Channel.of(
            'chr1','chr2','chr3','chr4','chr5',
            'chr6','chr7','chr8','chr9','chr10',
            'chr11','chr12','chr13','chr14','chr15',
            'chr16','chr17','chr18','chr19','chr20',
            'chr21','chr22','chrX','chrY'
        )

        ch_longtr_in = aligned_samples
            .combine(ch_chroms)
            .map { sample, aln, idx, karyotype, chrom ->
                tuple(sample, chrom, aln, idx, karyotype)
            }

        ch_longtr_out = LONGTR_PER_CHROM(ch_longtr_in, ref, fai)

        ch_for_merge = ch_longtr_out
            .map { sample, chrom, vcf -> tuple(sample, vcf) }
            .groupTuple()

        called_vcfs = MERGE_LONGTR(ch_for_merge)
        final_reports = ANNOTATE(called_vcfs)
        final_reports.view { "Finished LongTR per chromosome + merge + annotation: ${it}" }
    }
    else {
        error "Unknown caller '${params.caller}'. Use: atarva | longtr | longtr_chrom | strdust | medaka"
    }

    ch_merge_reports = final_reports
        .map { sample, tsv -> tsv }
        .collect()

    merged_report = MERGE_STR_REPORTS(ch_merge_reports)
    merged_report.view { "Finished merged STR report: ${it}" }
}

// ═════════════════════════════════════════════════════════════════════════════
// STR CALLERS
// ═════════════════════════════════════════════════════════════════════════════

process ATARVA {
    tag "${sample}"

    conda "${params.atarva_conda}"

    cpus   { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${atarvaOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.atarva.vcf"), emit: vcf

    script:
    def formatArg = aln.name.endsWith('.cram') ? 'cram' : 'bam'
    """
    ${params.atarva_bin} \\
        genotype \\
        -t ${task.cpus} \\
        -f ${ref} \\
        -b ${aln} \\
        -r ${params.atarva_bed} \\
        --format ${formatArg} \\
        --karyotype ${karyotype} \\
        -o ${sample}.atarva.vcf
    """
}

process LONGTR {
    tag "${sample}"

    container 'community.wave.seqera.io/library/longtr:1.2--3a7af9434e146eab'

    cpus 1
    memory { 8.GB * task.attempt }
    time { 8.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${longtrOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.longTR.vcf.gz")

    script:
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''

    """
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${params.longtr_bed} | sort -n | tail -n 1)";

    LongTR --help

    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --min-reads 1 \\
        --max-tr-len \$MAX_TR_LEN \\
        --min-mean-qual 1 \\
        --regions ${params.longtr_bed} \\
        ${haploid_args} \\
        --bams ${aln} \\
        --bam-samps ${sample} \\
        --bam-libs ${sample} \\
        --tr-vcf ${sample}.longTR.vcf.gz
    """
}

process LONGTR_PER_CHROM {
    tag "${sample}:${chrom}"

    container 'community.wave.seqera.io/library/longtr:1.2--3a7af9434e146eab'

    cpus   1
    memory { 16.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${longtrChromOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), val(chrom), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    tuple val(sample), val(chrom), path("${sample}.${chrom}.longTR.vcf.gz")

    script:
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''

    """
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${params.longtr_bed} | sort -n | tail -n 1)"

    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --max-tr-len \$MAX_TR_LEN \\
        --min-mean-qual 1 \\
        --regions ${params.longtr_bed} \\
        ${haploid_args} \\
        --chrom ${chrom} \\
        --bams ${aln} \\
        --bam-samps ${sample} \\
        --bam-libs ${sample} \\
        --tr-vcf ${sample}.${chrom}.longTR.vcf.gz
    """
}

process MERGE_LONGTR {
    tag "${sample}"

    container 'community.wave.seqera.io/library/bcftools_samtools:3f506bc690e52c7d'

    cpus   1
    memory { 2.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${longtrChromOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(vcfs)

    output:
    tuple val(sample), path("${sample}.longTR.vcf.gz")

    script:
    """
    for v in ${vcfs}; do
        if [ ! -f "\${v}.tbi" ] && [ ! -f "\${v}.csi" ]; then
            tabix -p vcf "\${v}"
        fi
    done

    bcftools concat -a -Oz -o ${sample}.longTR.vcf.gz ${vcfs}
    tabix -p vcf ${sample}.longTR.vcf.gz
    """
}

process STRDUST {
    tag "${sample}"

    cpus   { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${strdustOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.strdust.sorted.vcf"), emit: vcf

    script:
    def haploid = (karyotype == 'XY') ? "--haploid chrX,chrY" : ""

    """
    ${params.strdust_bin} \\
        -R ${params.strdust_bed} \\
        --unphased \\
        -t ${task.cpus} \\
        --sample ${sample} \\
        ${haploid} \\
        ${ref} ${aln} > ${sample}.strdust.vcf

    {
        grep '^#' ${sample}.strdust.vcf
        grep -v '^#' ${sample}.strdust.vcf | sort -t "\$(printf '\\t')" -k1,1V -k2,2n
    } > ${sample}.strdust.sorted.vcf
    """
}

process MEDAKA {
    tag "${sample}"

    container "${containerDir}/medaka:2.2.1--py310h237e959_0.sif"

    cpus   { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${medakaOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.medaka.vcf"), emit: vcf

    script:
    def sex = (karyotype == 'XX') ? 'female' : (karyotype == 'XY') ? 'male' : 'unknown'

    """
    set -euo pipefail

    medaka tandem \\
        ${aln} \\
        ${ref} \\
        ${params.medaka_bed} \\
        ${sex} \\
        ${sample}.medaka.output \\
        --workers ${task.cpus} \\
        --sample_name ${sample} \\
        --debug

    if [ ! -f "${sample}.medaka.output/medaka_to_ref.TR.vcf" ]; then
        echo "ERROR: Expected Medaka output VCF not found:"
        echo "${sample}.medaka.output/medaka_to_ref.TR.vcf"
        echo "Available files:"
        find "${sample}.medaka.output" -maxdepth 2 -type f -print
        exit 1
    fi

    cp "${sample}.medaka.output/medaka_to_ref.TR.vcf" "${sample}.medaka.vcf"
    """
}

// ═════════════════════════════════════════════════════════════════════════════
// UNIVERSAL ANNOTATION
// ═════════════════════════════════════════════════════════════════════════════

process ANNOTATE {
    tag "${sample}"

    container "${containerDir}/annotsv_3.5.8_biocontainers.sif"

    cpus 2
    memory { 8.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 0

    publishDir (
        params.caller == 'longtr_chrom' ? longtrChromOut :
        params.caller == 'longtr'       ? longtrOut :
        params.caller == 'atarva'       ? atarvaOut :
        params.caller == 'strdust'      ? strdustOut :
        params.caller == 'medaka'       ? medakaOut :
        "${params.work_dir}/annotate"
    ), mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(vcf)

    output:
    tuple val(sample), path("${sample}.${params.caller}.final.tsv")

    script:
    def caller_label =
        params.caller == 'longtr_chrom' ? 'longTR' :
        params.caller == 'longtr'       ? 'longTR' :
        params.caller == 'atarva'       ? 'atarva' :
        params.caller == 'strdust'      ? 'strdust' :
        params.caller == 'medaka'       ? 'medaka' :
                                          params.caller

    def json_catalog =
        params.caller == 'longtr_chrom' ? params.longtr_json :
        params.caller == 'longtr'       ? params.longtr_json :
        params.caller == 'atarva'       ? params.atarva_json :
        params.caller == 'strdust'      ? params.strdust_json :
        params.caller == 'medaka'       ? params.medaka_json :
                                          params.longtr_json

    """
    set -euo pipefail

    input_vcf="${vcf}"
    filtered_vcf="${sample}.${caller_label}.filtered.vcf"

    if [[ "\${input_vcf}" == *.gz ]]; then
        if [ ! -f "\${input_vcf}.tbi" ] && [ ! -f "\${input_vcf}.csi" ]; then
            tabix -p vcf "\${input_vcf}"
        fi
    fi

    # Keep VCF header and remove records where the sample GT is missing/reference-only:
    # ., ./., .|., 0, 0/0, or 0|0
    {
        if [[ "\${input_vcf}" == *.gz ]]; then
            gzip -cd "\${input_vcf}"
        else
            cat "\${input_vcf}"
        fi
    } | awk -F'\t' -v OFS='\t' '
    BEGIN {
        kept = 0
        removed = 0
    }
    (/^#/) {
        print
        next
    }
    {
        if (NF < 10) {
            print
            kept++
            next
        }

        nfmt = split(\$9, fmt, ":")
        gt_idx = 0

        for (i = 1; i <= nfmt; i++) {
            if (fmt[i] == "GT") {
                gt_idx = i
                break
            }
        }

        if (gt_idx == 0) {
            print
            kept++
            next
        }

        nsample = split(\$10, sample_fields, ":")
        gt = sample_fields[gt_idx]

        if (gt == "." || gt == "./." || gt == ".|." || gt == "0" || gt == "0/0" || gt == "0|0") {
            removed++
            next
        }

        print
        kept++
    }
    END {
        print "VCF GT filter kept " kept " records and removed " removed " records with missing/reference-only GT" > "/dev/stderr"
    }
    ' > "\${filtered_vcf}"

    ${params.stranno_bin} anno \\
        "\${filtered_vcf}" \\
        ${json_catalog} \\
        --log-level debug \\
        --max-dist 10 \\
        -o ${sample}.${caller_label}.annotated.vcf

    python3 ${params.z_score_script} \\
        ${sample}.${caller_label}.annotated.vcf \\
        ${sample} \\
        ${sample}.${caller_label}.annotated.z_score.vcf

    python3 ${params.vcf_to_bed_script} \\
        ${sample}.${caller_label}.annotated.z_score.vcf \\
        ${sample}.${caller_label}.bed

    AnnotSV \\
        -SVinputFile ${sample}.${caller_label}.bed \\
        -svtBEDcol 4 \\
        -SVminSize 1 \\
        -SVinputInfo 1 \\
        -includeCI 0 \\
        -overlap 80 \\
        -outputDir . \\
        -annotationsDir ${params.software_dir}/tools/AnnotSV/share/AnnotSV/ \\
        -outputFile ${sample}.${caller_label}.annotated.tsv

    java -jar -Xmx16G ${params.htssidra_jar} \\
        processAnnotatedFileAnnotSV_3_8 \\
        ${sample}.${caller_label}.annotated.tsv

    python3 ${params.vcf_to_tsv_script} \\
        -i ${sample}.${caller_label}.annotated.z_score.vcf \\
        -o ${sample}.${caller_label}.annotated.z_score.tsv

    python3 ${params.left_join_script} \\
        ${sample}.${caller_label}.annotated.z_score.tsv \\
        ${sample}.${caller_label}.annotated.tsv.processed.tsv \\
        --left_key ID \\
        --right_key ID \\
        --out ${sample}.${params.caller}.final.tsv
    """
}

// ═════════════════════════════════════════════════════════════════════════════
// MERGE ALL SAMPLE FINAL TSV REPORTS
// ═════════════════════════════════════════════════════════════════════════════

process MERGE_STR_REPORTS {
    tag "${params.projectname}:${params.caller}"

    cpus 1
    memory '1 GB'
    time '1h'
    errorStrategy 'terminate'

    publishDir (
        params.caller == 'longtr_chrom' ? longtrChromOut :
        params.caller == 'longtr'       ? longtrOut :
        params.caller == 'atarva'       ? atarvaOut :
        params.caller == 'strdust'      ? strdustOut :
        params.caller == 'medaka'       ? medakaOut :
        "${params.work_dir}/merged"
    ), mode: 'copy'

    input:
    path final_tsvs

    output:
    path "${params.projectname}.merged.tsv", emit: merged

    script:
    """
    set -euo pipefail

    OUTFILE="${params.projectname}.merged.tsv"
    rm -f "\${OUTFILE}"

    first_file=true
    n_files=0

    for tsv in ${final_tsvs}; do
        [ -n "\${tsv}" ] || continue

        if [ ! -s "\${tsv}" ]; then
            echo "WARNING: Missing or empty final TSV: \${tsv}" >&2
            continue
        fi

        sample="\$(basename "\${tsv}")"
        sample="\${sample%.${params.caller}.final.tsv}"

        if [ "\${first_file}" = true ]; then
            awk -v sample="\${sample}" 'BEGIN{OFS="\t"} NR==1 {print "sample", \$0; next} {print sample, \$0}' "\${tsv}" > "\${OUTFILE}"
            first_file=false
        else
            awk -v sample="\${sample}" 'BEGIN{OFS="\t"} NR>1 {print sample, \$0}' "\${tsv}" >> "\${OUTFILE}"
        fi

        n_files=\$((n_files + 1))
    done

    if [ "\${n_files}" -eq 0 ]; then
        echo "ERROR: No final TSV files were available to merge" >&2
        echo "PWD: \$(pwd)" >&2
        echo "Directory contents:" >&2
        ls -lah >&2 || true
        exit 1
    fi

    echo "Merged \${n_files} final TSV files into: \${OUTFILE}"
    """
}
