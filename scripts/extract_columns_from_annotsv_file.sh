INPUT="/pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.annotated.tsv"
OUTPUT="/pl/active/dashnowlab/software/ONT_PIPELINE/catalogs/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.atarva.annotated.subset.tsv"

awk -F'\t' -v OFS='\t' '
BEGIN {
    cols = "AnnotSV_ID SV_chrom SV_start SV_end SV_length SV_type ID CytoBand Gene_name Tx_merged Tx_version_merged Tx_start_merged Tx_end_merged Overlapped_tx_length_merged Overlapped_CDS_length_merged Overlapped_CDS_percent_merged Frameshift_merged Exon_count_merged Location_merged Location2_merged ACMG_merged HI_merged TS_merged DDD_HI_percent_merged ExAC_delZ_merged ExAC_dupZ_merged ExAC_cnvZ_merged ExAC_synZ_merged ExAC_misZ_merged GenCC_disease_merged GenCC_moi_merged GenCC_classification_merged GenCC_pmid_merged NCBI_gene_ID_merged HGNC_gene_ID_merged OMIM_ID_merged OMIM_phenotype_merged OMIM_inheritance_merged OMIM_morbid_merged OMIM_morbid_candidate_merged LOEUF_bin_merged GnomAD_pLI_merged ExAC_pLI_merged highest_PLI Genic Exonic Intronic Novel pLOF DISEASE_GENES"

    n = split(cols, wanted, " ")
}
NR == 1 {
    sub(/\r$/, "")
    for (i = 1; i <= NF; i++) {
        index_of[$i] = i
    }

    for (i = 1; i <= n; i++) {
        if (!(wanted[i] in index_of)) {
            print "ERROR: Missing column: " wanted[i] > "/dev/stderr"
            missing = 1
        }
    }

    if (missing) exit 1

    for (i = 1; i <= n; i++) {
        printf "%s%s", wanted[i], (i < n ? OFS : ORS)
    }
    next
}
{
    for (i = 1; i <= n; i++) {
        printf "%s%s", $(index_of[wanted[i]]), (i < n ? OFS : ORS)
    }
}
' "$INPUT" > "$OUTPUT"