#' Extract CCLE drug response data
#'
#' This function creates a \code{data.frame} containing the CCLE drug response data from the database.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
getDrugData_CCLE <- function(con, drugs, cell_lines) {

  drugs.sql <- paste(drugs, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select CCLE_name as unified_id, Compound as assayed_id, 'resp' as data_type, EC50_uM as original
               from ccle_drug_data
               where CCLE_name IN ('%s') and Compound IN ('%s')", cell_lines.sql, drugs.sql)
  out <- DBI::dbGetQuery(con, sql)
  out <- out %>% dplyr::mutate(value=6-log10(original)) %>%
    dplyr::mutate_each(dplyr::funs(as.character), -value) %>%
    dplyr::mutate_each(dplyr::funs(as.numeric), value)
  return(out)

}
