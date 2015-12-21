#' Extract copy number data
#'
#' This function creates a \code{data.frame} containing the by gene copy number data from the database for the requested cell lines and genes.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vectore of cell line identifiers
#' @return A \code{data.frame} containing the by gene copy number data for the requested compounds and cell lines
#' @export
getCopyNumberData <- function(con, genes, cell_lines) {

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select CCLE_name, Symbol as ID, 'cn' as Type, log2cn as original, log2cn as value
                 from ccle_cn
                 where CCLE_name IN ('%s') and Symbol IN ('%s')", cell_lines.sql, genes.sql)

  data <- DBI::dbGetQuery(con, sql)

  #duplicate remove
  data <- data %>% dplyr::group_by(ID) %>%
    dplyr::mutate(N=n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(N==median(N)) %>%
    dplyr::select(-N) %>%
    as.data.frame

  data <- data %>% dplyr::mutate_each(dplyr::funs(as.character), -value) %>%
    dplyr::mutate_each(dplyr::funs(as.numeric), value)

  return(data)

}
