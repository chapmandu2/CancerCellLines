#' Get data for shiny
#'
#' This function creates a \code{data.frame} suitable for the shiny app
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
makeRespVsGeneticDataFrame <- function(con, gene, cell_lines, drug, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'), drug_df=NULL) {

  #gene <- 'TP53'
  #drug <- 'GSKMI00000714_pIC50'
  #cell_lines <- cls.df$CCLE_name
  #data_types <- c('affy', 'hybcap')
  #drug_df <- dnmt1_data

  gene_data <- makeTallDataFrame(con,
                                 gene,
                                 cell_lines,
                                 drug,
                                 data_types) %>%
    transmute(CCLE_name, ID, feature_type=Type, feature_name=paste(ID, Type, sep="_"), feature_value=value)

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

  cls.df <- src_sqlite(con@dbname) %>%
    tbl('ccle_sampleinfo') %>%
    transmute(CCLE_name, tissue=Site_primary) %>%
    collect

  plot_data <- gene_data %>%
    inner_join(resp_data, by='CCLE_name') %>%
    inner_join(cls.df, by='CCLE_name')

  return(plot_data)
}
