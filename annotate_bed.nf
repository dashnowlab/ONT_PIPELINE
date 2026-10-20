#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TRExplorer BED6(.gz) -> BED5 -> AnnotSV per batch -> htsSidra per batch -> merge
params.input_bed       = params.containsKey('input_bed')       ? params.input_bed       : null
params.work_dir        = params.containsKey('work_dir')        ? params.work_dir        : launchDir
params.software_path   = params.containsKey('software_path')   ? params.software_path   : projectDir
params.batch_size      = (params.containsKey('batch_size')      ? params.batch_size      : 10000) as Integer
params.annotsv_mem_gb  = (params.containsKey('annotsv_mem_gb')  ? params.annotsv_mem_gb  : 4) as Integer
params.java_mem_gb     = (params.containsKey('java_mem_gb')     ? params.java_mem_gb     : 16) as Integer
params.java_heap_gb    = (params.containsKey('java_heap_gb')    ? params.java_heap_gb    : 12) as Integer
params.java_max_forks  = (params.containsKey('java_max_forks')  ? params.java_max_forks  : 8) as Integer

if (!params.input_bed) error 'Missing --input_bed <catalog.bed.gz>'
if (params.batch_size < 1) error '--batch_size must be positive'
if (params.java_max_forks < 1) error '--java_max_forks must be positive'
if (params.annotsv_mem_gb < 1 || params.java_mem_gb < 1 || params.java_heap_gb < 1) error 'Memory values must be positive'
if (params.java_heap_gb >= params.java_mem_gb) error '--java_mem_gb must exceed --java_heap_gb'

def containerDir = "${params.software_path}/containers"
// --work_dir is already the catalog-specific .tmp directory from the SLURM runner.
def tmpDir = params.work_dir.toString()
def finalDir = file(params.work_dir).parent.toString()

workflow {
    def catalog = file(params.input_bed, checkIfExists: true)
    def project = catalog.name.replaceFirst(/\.gz$/, '').replaceFirst(/\.bed$/, '')

    prepared = PREPARE_BED(Channel.of(tuple(project, catalog)))
    chunks = SPLIT_BED(prepared)
    batch_inputs = chunks.flatMap { project_id, beds ->
        beds.collect { bed -> tuple(project_id, bed.baseName, bed) }
    }
    annotated = ANNOTATE_SV(batch_inputs)
    processed = PROCESS_ANNOTSV_BATCH(annotated)
    grouped = processed.groupTuple(by: 0)
    MERGE_PROCESSED(grouped)
}

process PREPARE_BED {
    tag "${project}"
    cpus 1
    memory '4 GB'
    time '4h'
    publishDir "${tmpDir}", mode: 'copy', pattern: '*.bed'

    input:
    tuple val(project), path(input_bed)

    output:
    tuple val(project), path("${project}.bed")

    script:
    """
    set -euo pipefail
    case '${input_bed}' in
        *.gz) gzip -cd '${input_bed}' ;;
        *) cat '${input_bed}' ;;
    esac | awk -F '\t' 'BEGIN { OFS="\t"; print "#CHROM", "START", "END", "SVTYPE", "ID" }
        /^#/ || /^track/ || /^browser/ || NF == 0 { next }
        {
            if (NF < 6 || \$1 == "" || \$6 == "" || \$6 == "." ||
                \$2 == "" || \$3 == "" || \$2 ~ /[^0-9]/ || \$3 ~ /[^0-9]/ || \$3 <= \$2) {
                print "ERROR: invalid BED6 record at input line " NR ": " \$0 > "/dev/stderr"
                failed=1; exit 1
            }
            print \$1, \$2, \$3, "DEL", \$6
            count++
        }
        END { if (failed) exit 1; if (count == 0) { print "ERROR: no BED records" > "/dev/stderr"; exit 1 } }
    ' > '${project}.bed'
    awk 'END { print "Prepared " NR-1 " loci" > "/dev/stderr" }' '${project}.bed'
    """
}

process SPLIT_BED {
    tag "${project}"
    cpus 1
    memory '2 GB'
    time '4h'
    publishDir "${tmpDir}/batches", mode: 'copy', pattern: 'batch_*.bed'

    input:
    tuple val(project), path(bed)

    output:
    tuple val(project), path('batch_*.bed', arity: '1..*')

    script:
    """
    set -euo pipefail
    awk -v n=${params.batch_size} '
        NR == 1 { header=\$0; next }
        {
            batch=int((NR-2)/n)+1
            f=sprintf("batch_%06d.bed", batch)
            if (f != previous) {
                if (previous != "") close(previous)
                print header > f
                previous = f
            }
            print >> f
        }
    ' '${bed}'
    ls batch_*.bed >/dev/null
    """
}

process ANNOTATE_SV {
    tag "${project}:${batch_id}"
    container "${containerDir}/annotsv_3.5.8_biocontainers.sif"
    cpus 2
    memory { params.annotsv_mem_gb.GB }
    time '8h'
    errorStrategy 'terminate'
    maxForks 32
    publishDir "${tmpDir}/annotations", mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(project), val(batch_id), path(batch_bed)

    output:
    tuple val(project), val(batch_id), path("${batch_id}.tsv")

    script:
    def sw = params.software_path
    """
    set -euo pipefail
    AnnotSV \\
        -SVinputFile '${batch_bed}' \\
        -svtBEDcol 4 \\
        -SVminSize 1 \\
        -SVinputInfo 1 \\
        -includeCI 0 \\
        -overlap 80 \\
        -genomeBuild GRCh38 \\
        -overwrite 1 \\
        -outputDir . \\
        -annotationsDir '${sw}/tools/AnnotSV/share/AnnotSV/' \\
        -outputFile '${batch_id}.tsv'
    test -s '${batch_id}.tsv'
    """
}

process PROCESS_ANNOTSV_BATCH {
    tag "${project}:${batch_id}"
    cpus 2
    memory { params.java_mem_gb.GB }
    time '8h'
    errorStrategy 'terminate'
    maxForks params.java_max_forks
    publishDir "${tmpDir}/processed", mode: 'copy', pattern: '*.processed.tsv'

    input:
    tuple val(project), val(batch_id), path(annotated_tsv)

    output:
    tuple val(project), val(batch_id), path("${batch_id}.processed.tsv")

    script:
    def sw = params.software_path
    """
    set -euo pipefail
    java -Xmx${params.java_heap_gb}G -jar \\
        '${sw}/tools/htsSidra-1.2-jar-with-dependencies.jar' \\
        processAnnotatedFileAnnotSV_3_8 \\
        '${annotated_tsv}'
    test -s '${batch_id}.tsv.processed.tsv'
    mv '${batch_id}.tsv.processed.tsv' '${batch_id}.processed.tsv'
    """
}

process MERGE_PROCESSED {
    tag "${project}"
    cpus 1
    memory '4 GB'
    time '8h'
    publishDir "${finalDir}", mode: 'copy', pattern: '*.annotated.tsv'

    input:
    tuple val(project), val(batch_ids), path(processed_tsvs)

    output:
    path("${project}.annotated.tsv")

    script:
    // Batch identifiers are zero-padded, so lexical sorting preserves input order.
    def ordered = [batch_ids, processed_tsvs].transpose().sort { a, b -> a[0] <=> b[0] }
    def files = ordered.collect { "'${it[1].name}'" }.join(' ')
    """
    set -euo pipefail
    files=( ${files} )
    test \${#files[@]} -gt 0
    header=\$(head -n 1 "\${files[0]}")
    printf '%s\\n' "\$header" > '${project}.annotated.tsv'
    for f in "\${files[@]}"; do
        test -s "\$f" || { echo "ERROR: empty processed batch \$f" >&2; exit 1; }
        test "\$(head -n 1 "\$f")" = "\$header" || { echo "ERROR: header mismatch \$f" >&2; exit 1; }
        tail -n +2 "\$f" >> '${project}.annotated.tsv'
    done
    echo "Merged \${#files[@]} processed batches" >&2
    """
}
