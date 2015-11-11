#' Import cell line interconversion data
#'
#' This function imports the information in the CellLineIDNormalisationOct15.txt file into the cell_line_ids table in the database.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCellLineIDs <- function (fn , con) {

  data <- read.table(fn, header=T, sep='\t')

  dbWriteTable(con, "cell_line_ids", data, overwrite=TRUE)
}
