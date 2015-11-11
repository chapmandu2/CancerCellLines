#' Extract affymetrix expression data
#'
#' This function creates a \code{data.frame} containing the affymetrix  gene expression data from the database for the requested cell lines and genes.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vectore of cell line identifiers
#' @return A \code{data.frame} containing the affymetrix gene expression data for the requested compounds and cell lines
#' @export
getAffyData <- function(con, genes, cell_lines) {

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select CCLE_name, Symbol as ID, 'affy' as Type, Signal as original, Signal as value
               from ccle_affy
               where CCLE_name IN ('%s') and Symbol IN ('%s')", cell_lines.sql, genes.sql)

  data <- dbGetQuery(con, sql)

  #duplicate remove
  data <- data %>% group_by(ID) %>% mutate(N=n()) %>% ungroup %>% filter(N==median(N)) %>% select(-N) %>% as.data.frame

  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)

}
