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
#' @param tissue_info Source of tissue information: custom can be one of any id_type present in the cell_line_ids table: \code{ccle, custom}
#' @return A \code{data.frame} containing the drug response and genetic data for the requested compounds and cell lines
#' @export
makeRespVsGeneticDataFrame <- function(con, gene, cell_lines, drug, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'), drug_df=NULL, tissue_info='ccle') {

  gene_data <- makeTallDataFrame(con,
                                 gene,
                                 cell_lines,
                                 drug,
                                 data_types) %>%
    transmute(CCLE_name, ID, feature_type=Type, feature_name=paste(ID, Type, sep="_"), feature_value=value, feature_original=original)

  if (is.null(drug_df)) {
    #if no drug data frame provided just use CCLE
    ccle_drugs.df <- src_sqlite(con@dbname) %>% tbl("ccle_drug_data") %>% select(Compound) %>% distinct %>% collect
    resp_data <- getDrugData_CCLE(con, drug, cell_lines) %>%
      transmute(CCLE_name, resp_id=ID, resp_value=value)
  } else {
    #if drug data provided then use that
    resp_data <- getDrugData_custom(drug_df, drug, cell_lines ) %>%
      transmute(CCLE_name, resp_id=ID, resp_value=value)
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

  plot_data <- gene_data %>%
    inner_join(resp_data, by='CCLE_name') %>%
    inner_join(cls_df, by='CCLE_name')

  return(plot_data)
}
