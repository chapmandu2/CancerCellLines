#' Get tissue info for cell lines
#'
#' This function creates a \code{data.frame} of cell lines with associated tissue information which is derived in different ways depending on the tissue_info parameter
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param tissue_info Source of tissue information: \code{ccle} uses the ccle sample information, \code{crude} ensures that all cell lines are returned but mixes data from different sources, \code{custom} can be one of any id_type present in the cell_line_ids table for extra control
#' @return A \code{data.frame} containing tissue information for available cell lines
#' @export
getTissueInfo <- function (con, tissue_info='ccle') {

  if (tissue_info == 'ccle') {
    cls_df <- dplyr::src_sqlite(con@dbname) %>%
      dplyr::tbl('ccle_sampleinfo') %>%
      dplyr::transmute(CCLE_name, native_id=Primary_cell_name, tissue=Site_primary, subtype1=Histology,  subtype2=Hist_subtype1) %>%
      dplyr::collect()
  } else  if (tissue_info == 'crude') {
    id_type_order_df <- data.frame(id_type=c('CCLE', 'cosmic_clp'), id_type_order=c(1,2), stringsAsFactors = FALSE)

    cls_df <- dplyr::src_sqlite(con@dbname) %>%
      dplyr::tbl('cell_line_ids') %>%
      dplyr::collect() %>%
      dplyr::transmute(CCLE_name=unified_id, native_id, tissue, subtype1=hist_primary, subtype2=hist_secondary, id_type) %>%
      dplyr::left_join(id_type_order_df, by='id_type') %>%
      dplyr::group_by(CCLE_name) %>%
      dplyr::mutate(rank=dplyr::min_rank(id_type_order)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(rank == 1) %>%
      dplyr::select(CCLE_name, native_id, tissue, subtype1, subtype2) %>%
      dplyr::arrange(CCLE_name)

  } else {
    cls_df <- dplyr::src_sqlite(con@dbname) %>%
      dplyr::tbl('cell_line_ids') %>%
      dplyr::filter(id_type == tissue_info) %>%
      dplyr::transmute(CCLE_name=unified_id, native_id, tissue, subtype1=hist_primary, subtype2=hist_secondary) %>%
      dplyr::collect()
  }

  if(nrow(cls_df)==0) {
    stop('tissue_info not valid - see ?getTissueInfo for details')
  }

  return(as.data.frame(cls_df))

}
