
#' Import CCLE hybrid capture sequencing data
#'
#' This function imports the information in CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf into the ccle_hybap table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_hybcap <- function ( fn , con ) {


  #fn <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf'

  data <- readr::read_tsv(fn,
                    col_names=c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
                                "Start_position", "End_position", "Strand", "Variant_Classification", "Variant_Type",
                                "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",  "dbSNP_RS", "dbSNP_Val_Status",
                                "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Transcript_Strand", "cDNA_Change", "Codon_Change",
                                "Protein_Change"),
                    col_types=paste(c('cicicii', rep('c', 10), rep('_', 17), rep('c', 4), rep('_', 13)), collapse=''),
                    skip=1)


  DBI::dbWriteTable(con, "ccle_hybcap", as.data.frame(data), overwrite=TRUE)

}
