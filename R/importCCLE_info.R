#' Import CCLE cell line information
#'
#' This function imports the information in CCLE_sample_info_file_2012-10-18.txt into the database.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_info <- function ( fn , con )  {

  data <- read.table(fn,header=T, sep='\t')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Cell_line_aliases', 'Gender', 'Site_primary', 'Histology', 'Hist_subtype1', 'Notes', 'Source', 'Expression_arrays', 'SNP_arrays', 'Oncomap', 'Hybrid_capture_sequencing')

  dbWriteTable(con, "ccle_sampleinfo", data, overwrite=TRUE)

}
