#' Get RespVsGenetic data
#'
#' This function creates a \code{data.frame} suitable for the RespVsGenetic plots and shiny app
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param gene The gene symbol for the gene of interest
#' @param cell_lines A vector of cell line identifiers
#' @param drug The drug of interest
#' @param data_types A vector with default \code{c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp')} to specify which data types should be returned.
#' @param drug_df A \code{data.frame} containing input for the \code{getDrugData_custom} function
#' @param tissue_info Source of tissue information - see ?getTissueInfo
#' @return A \code{data.frame} containing the drug response and genetic data for the requested compounds and cell lines
#' @export
makeRespVsGeneticDataFrame <- function(con, gene, cell_lines, drug, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'), drug_df=NULL, tissue_info='ccle') {

  gene_data <- makeTallDataFrame(con,
                                 gene,
                                 cell_lines,
                                 drug,
                                 data_types) %>%
    dplyr::transmute(unified_id, assayed_id, feature_type=data_type, feature_name=paste(assayed_id, data_type, sep="_"), feature_value=value, feature_original=original)

  if (is.null(drug_df)) {
    #if no drug data frame provided just use CCLE
    ccle_drugs.df <- dplyr::src_sqlite(con@dbname) %>%
      dplyr::tbl("ccle_drug_data") %>%
      dplyr::select(Compound) %>%
      dplyr::distinct() %>%
      dplyr::collect()
    resp_data <- getDrugData_CCLE(con, drug, cell_lines) %>%
      dplyr::transmute(unified_id, resp_id=assayed_id, resp_value=value)
  } else {
    #if drug data provided then use that
    resp_data <- getDrugData_custom(drug_df, drug, cell_lines ) %>%
      dplyr::transmute(unified_id, resp_id=assayed_id, resp_value=value)
  }

  #now get the tissue info depending on selected option
  cls_df <- getTissueInfo(con, tissue_info)

  plot_data <- gene_data %>%
    dplyr::inner_join(resp_data, by='unified_id') %>%
    dplyr::inner_join(cls_df, by='unified_id')

  return(plot_data)
}
