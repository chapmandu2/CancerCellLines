#' Generic function to get response data from multiple sources
#'
#' This function creates a \code{data.frame} containing the response data from a variety of sourses.
#'
#' @param src A \code{SQLiteConnection} object to the database (ccle) or \code{data.frame} object (custom)
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @param resp_type The type of response data: \code{ccle, custom}
#' @return A \code{data.frame} containing the response data for the requested compounds and cell lines
#' @export
getResponseData <- function(src, drugs, cell_lines, resp_type='ccle') {

  if (class(src) == 'SQLiteConnection' & resp_type == 'ccle') {
    getDrugData_CCLE(src, drugs, cell_lines)
  } else if (class(src) == 'data.frame' & resp_type == 'custom') {
    getDrugData_custom(src, drugs, cell_lines)
  } else {
    stop('Check ?getResponseData to ensure you are passing the correct options into src')
  }

}
