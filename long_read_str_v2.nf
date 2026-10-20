#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// Params & Help
// ─────────────────────────────────────────────────────────────────────────────
params.help              = params.containsKey('help')              ? params.help              : false
params.bam_list          = params.containsKey('bam_list')          ? params.bam_list          : null
params.work_dir          = params.containsKey('work_dir')          ? params.work_dir          : launchDir
params.projectname       = params.containsKey('projectname')       ? params.projectname       : "str_project"
params.ref               = params.containsKey('ref')               ? params.ref               : "/gpfs/data_jrnas1/ref_data/Hsapiens/hg38/ONT/hg38_nohla.fa"
params.base_dir          = params.containsKey('base_dir')          ? params.base_dir          : null

// Supplying --preannotated_tsv bypasses AnnotSV and htsSidra.
params.preannotated_tsv  = params.containsKey('preannotated_tsv') ? params.preannotated_tsv : null

// Software root
params.software_dir      = params.containsKey('software_dir')      ? params.software_dir      : projectDir

// Catalogs / resources
params.atarva_bed        = params.containsKey('atarva_bed')        ? params.atarva_bed        : "/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/tr_explorer_catalog/TR_catalog_for_atarva.bed.gz"

// Caller-specific JSON catalogs for stranno annotation
// Update these paths if your JSON catalog files have different names.
params.atarva_json       = params.containsKey('atarva_json')       ? params.atarva_json       : "${params.software_dir}/catalogs/TR_catalog.63_loci.20260327_193020.json.gz"

// Software / envs
params.atarva_bin        = params.containsKey('atarva_bin')        ? params.atarva_bin        : "/projects/ealiyev@xsede.org/software/anaconda/envs/atarva-0.9.3/bin/atarva"
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
    DASHNOW LAB ATaRVa STR Calling Pipeline

    Usage:
        nextflow run atarva_str_annotation_modes_fixed.nf \
            --bam_list <bam_or_cram_and_karyotype.tsv> \
            --work_dir <output_dir> \
            --projectname <project_name> \
            --ref <reference.fa> \
            --software_dir <software_root> \
            --annotation_mode <run|precomputed> \
            --preannotated_tsv <catalog.annotated.tsv>  # optional

    Input TSV: <bam|cram_path>\t<XX|XY> (no header; # comments allowed)
    Workflow: ATaRVa -> stranno -> z-score -> [AnnotSV + htsSidra OR catalog TSV]
              -> VCF-to-TSV -> left join -> merge reports
    """
    System.exit(0)
}

// ─────────────────────────────────────────────────────────────────────────────
// Output sub-directories
// ─────────────────────────────────────────────────────────────────────────────
def atarvaOut      = "${params.work_dir}/atarva"
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
    if (!params.bam_list) error 'Missing --bam_list'
    Channel.fromPath(params.bam_list, checkIfExists: true)
        .splitText()
        .map { it.replaceAll(/\r$/, '') }
        .filter { it.trim() && !it.trim().startsWith('#') }
        .map { line ->
            def parts = line.split('\t', -1)
            if (parts.size() < 2) error "Expected BAM/CRAM and karyotype: ${line}"
            def aln = resolvePath(parts[0])
            def karyo = parts[1].trim()
            if (!(karyo in ['XX','XY'])) error "Invalid karyotype ${karyo}"
            tuple(sampleName(aln), aln, alignmentIndex(aln), karyo)
        }.set { aligned_samples }

    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)
    calls = ATARVA(aligned_samples, ref, fai)
    prepared = PREPARE_VARIANTS(calls)

    if (params.preannotated_tsv) {
        annotation = Channel.value(file(params.preannotated_tsv, checkIfExists: true))
    } else {
        beds = VCF_TO_BED(prepared)
        annotated = ANNOTSV(beds)
        annotation = HTSSIDRA(annotated)
    }

    // Each sample receives the same catalog file in precomputed mode.
    // In run mode, match per-sample annotation by sample ID.
    if (params.preannotated_tsv) {
        join_inputs = prepared.combine(annotation)
            .map { sample, vcf, tsv, catalog -> tuple(sample, tsv, catalog) }
    } else {
        join_inputs = prepared.map { sample, vcf, tsv -> tuple(sample, tsv) }
            .join(annotation, by: 0)
            .map { sample, tsv, catalog -> tuple(sample, tsv, catalog) }
    }
    reports = LEFT_JOIN(join_inputs)
    MERGE_STR_REPORTS(reports.map { sample, tsv -> tsv }.collect())
}

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
        --decompose \\
        --strict \\
        --aln-format ${formatArg} \\
        --karyotype ${karyotype} \\
        -o ${sample}.atarva.vcf
    """
}

process PREPARE_VARIANTS {
    tag "${sample}"
    cpus 2
    memory '8 GB'
    time '4h'
    publishDir "${atarvaOut}", mode: 'copy', saveAs: { filename -> "${sample}/${filename}" }
    input:
    tuple val(sample), path(vcf)
    output:
    tuple val(sample), path("${sample}.atarva.annotated.z_score.vcf"), path("${sample}.atarva.annotated.z_score.tsv")
    script:
    """
    set -euo pipefail
    input_vcf="${vcf}"
    filtered_vcf="${sample}.atarva.filtered.vcf"

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
        ${params.atarva_json} \\
        --log-level debug \\
        --max-dist 10 \\
        -o ${sample}.atarva.annotated.vcf

    python3 ${params.z_score_script} \\
        ${sample}.atarva.annotated.vcf \\
        ${sample} \\
        ${sample}.atarva.annotated.z_score.vcf

    python3 ${params.vcf_to_tsv_script} \
        -i ${sample}.atarva.annotated.z_score.vcf \
        -o ${sample}.atarva.annotated.z_score.tsv
    """
}

process VCF_TO_BED {
    tag "${sample}"
    cpus 1
    memory '2 GB'
    time '2h'
    input:
    tuple val(sample), path(vcf), path(tsv)
    output:
    tuple val(sample), path("${sample}.atarva.bed")
    script:
    """
    set -euo pipefail
    python3 ${params.vcf_to_bed_script} ${vcf} ${sample}.atarva.bed
    """
}

process ANNOTSV {
    tag "${sample}"
    container "${containerDir}/annotsv_3.5.8_biocontainers.sif"
    cpus 2
    memory '8 GB'
    time '8h'
    input:
    tuple val(sample), path(bed)
    output:
    tuple val(sample), path("${sample}.atarva.annotated.tsv")
    script:
    """
    set -euo pipefail
    AnnotSV \
        -SVinputFile ${bed} \
        -svtBEDcol 4 \
        -SVminSize 1 \
        -SVinputInfo 1 \
        -includeCI 0 \
        -overlap 80 \
        -genomeBuild GRCh38 \
        -outputDir . \
        -annotationsDir ${params.software_dir}/tools/AnnotSV/share/AnnotSV/ \
        -outputFile ${sample}.atarva.annotated.tsv
    """
}

process HTSSIDRA {
    tag "${sample}"
    cpus 2
    memory '20 GB'
    time '8h'
    input:
    tuple val(sample), path(annotated_tsv)
    output:
    tuple val(sample), path("${sample}.atarva.annotated.tsv.processed.tsv")
    script:
    """
    set -euo pipefail
    java -Xmx16G -jar ${params.htssidra_jar} \
        processAnnotatedFileAnnotSV_3_8 ${annotated_tsv}
    test -s ${sample}.atarva.annotated.tsv.processed.tsv
    """
}

process LEFT_JOIN {
    tag "${sample}"

    cpus 2
    memory '8 GB'
    time '4h'

    publishDir "${atarvaOut}", mode: 'copy',
        saveAs: { filename -> "${sample}/${filename}" }

    input:
    tuple val(sample), path(sample_tsv), path(annotation_tsv)

    output:
    tuple val(sample), path("${sample}.atarva.final.tsv")

    script:
    """
    set -euo pipefail

    python3 ${params.left_join_script} \\
        ${sample_tsv} \\
        ${annotation_tsv} \\
        --left_key ID \\
        --right_key ID \\
        --out ${sample}.atarva.final.tsv
    """
}

process MERGE_STR_REPORTS {
    tag "${params.projectname}:atarva"

    cpus 1
    memory '1 GB'
    time '1h'
    errorStrategy 'terminate'

    publishDir "${atarvaOut}", mode: 'copy'

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
    expected_header=""

    for tsv in ${final_tsvs}; do
        [ -n "\${tsv}" ] || continue

        if [ ! -s "\${tsv}" ]; then
            echo "WARNING: Missing or empty final TSV: \${tsv}" >&2
            continue
        fi

        sample="\$(basename "\${tsv}")"
        sample="\${sample%.atarva.final.tsv}"

        current_header="\$(head -n 1 "\${tsv}")"
        if [ "\${first_file}" = true ]; then
            expected_header="\${current_header}"
            awk -v sample="\${sample}" 'BEGIN{OFS="\t"} NR==1 {print "sample", \$0; next} {print sample, \$0}' "\${tsv}" > "\${OUTFILE}"
            first_file=false
        else
            if [ "\${current_header}" != "\${expected_header}" ]; then
                echo "ERROR: Header mismatch in \${tsv}" >&2
                exit 1
            fi
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
