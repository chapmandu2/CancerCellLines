
#' Merge different data types into a single wide data frame
#'
#' This function creates a \code{data.frame} containing data for different data types in a form suitable for further statistical modelling in R.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vector of cell line identifiers
#' @param drugs A vector of compound identifiers
#' @param data_types A vector with default \code{c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp')} to specify which data types should be returned.
#' @param drug_df A \code{data.frame} containing input for the \code{getDrugData_custom} function
#' @return A wide \code{data.frame} containing the requested data
#' @export
makeWideDataFrame <- function(con, genes, cell_lines, drugs, data_types=c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp'), drug_df=NULL) {

  #run makeTallDataFrame function
  all_data <- makeTallDataFrame(con, genes, cell_lines, drugs, data_types, drug_df)

  #create matrix
  output.df <- all_data %>% makeWideFromTallDataFrame()

  return(output.df)

}
