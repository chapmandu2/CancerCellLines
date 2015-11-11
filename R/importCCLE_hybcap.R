
#' Import CCLE hybrid capture sequencing data
#'
#' This function imports the information in CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx into the ccle_hybap table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_hybcap <- function ( fn , con ) {

  require(readxl)
  #fn <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx'

  data <- read_excel(fn, col_types=c('text', 'numeric', 'text', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                                     'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text'))

  dbWriteTable(con, "ccle_hybcap", as.data.frame(data), overwrite=TRUE)

}
