#' Make a RespVsResp data frame
#'
#' This function creates a \code{data.frame} containing the response data suitable for the RespVsResp plotting functions
#'
#' @param src A \code{SQLiteConnection} object to the database (ccle) or \code{data.frame} object (custom)
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @param resp_type The type of response data: \code{ccle, custom}
#' @param tissue_info Source of tissue information: custom can be one of any id_type present in the cell_line_ids table: \code{ccle, custom}
#' @param drug_df A \code{data.frame} containing input for the \code{getDrugData_custom} function
#' @return A \code{data.frame} containing the response data for the requested compounds and cell lines
#' @export
makeRespVsRespDataFrame <- function(con, drugs, cell_lines, resp_type='ccle', tissue_info='ccle', drug_df=NULL) {

  #use the getResponseData function to retrieve the data initially
  if (resp_type == 'ccle') {
    df <- getResponseData(con, drugs, cell_lines, resp_type)
  } else if (resp_type == 'custom') {
    df <- getResponseData(drug_df, drugs, cell_lines, resp_type)
  } else {
    stop('resp_type not valid - see ?makeRespVsRespDataFrame')
  }

  #now get the tissue info depending on selected option
  if (tissue_info == 'ccle') {
    cls_df <- src_sqlite(con@dbname) %>%
      tbl('ccle_sampleinfo') %>%
      transmute(CCLE_name, native_id=Primary_cell_name, tissue=Site_primary, subtype1=Histology,  subtype2=Hist_subtype1) %>%
      collect
  } else  {
    cls_df <- src_sqlite(con@dbname) %>%
      tbl('cell_line_ids') %>%
      filter(id_type == tissue_info) %>%
      transmute(CCLE_name=unified_id, native_id, tissue, subtype1=hist_primary, subtype2=hist_secondary) %>%
      collect
  }

  if(nrow(cls_df)==0) {
    stop('tissue_info not valid - see ?makeRespVsRespDataFrame')
  }

  #now combine
  out_df <- df %>% left_join(cls_df, by='CCLE_name')
  return(out_df)

}
