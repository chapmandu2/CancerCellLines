#' Convert tissue to cell lines
#'
#' This function creates a \code{vector} of cell lines from a tissue for a given type of tissue ifo
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param tissue_id A vector of tissue identifiers
#' @param tissue_info Source of tissue information: \code{ccle} uses the ccle sample information, \code{crude} ensures that all cell lines are returned but mixes data from different sources, \code{custom} can be one of any id_type present in the cell_line_ids table for extra control
#' @return A \code{vector} of cell lines for the provided tissue(s)
#' @export
getTissueCellLines <- function (con, tissue_id, tissue_info='ccle') {

  cls.df <- getTissueInfo(con, tissue_info) %>%
    dplyr::select(CCLE_name, tissue) %>%
    dplyr::filter(tissue %in% tissue_id)

  if(nrow(cls.df) == 0)  {
    warning('No cell lines found for that tissue')
    return(NULL)
  } else {
    return(cls.df$CCLE_name)
  }


}
