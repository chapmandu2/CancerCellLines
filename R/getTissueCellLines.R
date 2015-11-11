#' Convert tissue to cell lines
#'
#' This function creates a \code{vector} of cell lines from a tissue
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{vector} of cell lines for the provided tissue
#' @export
getTissueCellLines <- function (con, tissue_id) {

  cls.df <- src_sqlite(con@dbname) %>%
    tbl('ccle_sampleinfo') %>%
    transmute(CCLE_name, tissue=Site_primary) %>%
    collect %>%
    filter(tissue %in% tissue_id)


  if(nrow(cls.df) == 0)  {
    warning('No cell lines found for that tissue')
    return(NULL)
  } else {
    return(cls.df$CCLE_name)
  }


}
