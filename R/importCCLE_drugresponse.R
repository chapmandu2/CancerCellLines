#' Import CCLE drug response data
#'
#' This function imports the information in CCLE_NP24.2009_Drug_data_2012.02.20.csv into the ccle_drugresponse table in the database.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_drugresponse <- function (fn , con) {

  data <- read.table(fn, header=T, sep=',')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Compound', 'Target', 'Doses_uM', 'Activity_median', 'Activity_sd', 'Fit_type', 'Num_data', 'EC50_uM', 'IC50_uM', 'Amax', 'Act_area')

  dbWriteTable(con, "ccle_drug_data", data, overwrite=TRUE)
}
