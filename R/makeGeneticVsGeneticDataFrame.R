#' Get GeneticVsGenetic data
#'
#' This function creates a \code{data.frame} suitable for the GeneticVsGenetic plots and shiny app
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param cell_lines A vector of cell line identifiers
#' @param gene1 The gene symbol of the primary gene
#' @param gene2 The gene symbol of the secondary gene
#' @param data_type1 The data type required for gene1: should be one of \code{'affy', 'cn', 'hybcap', 'cosmicclp'}
#' @param data_type2 The data type required for gene1: should be one of \code{'affy', 'cn', 'hybcap', 'cosmicclp'}
#' @return A \code{data.frame} containing the combined genetic data for the requested cell lines
#' @export
makeGeneticVsGeneticDataFrame <- function(con, cell_lines, gene1, gene2, data_type1, data_type2) {

  gene1_data <- makeTallDataFrame(con,
                                 gene1,
                                 cell_lines,
                                 drugs=NULL,
                                 data_type1) %>%
    transmute(CCLE_name, gene1=ID, feature_type1=Type, feature_name1=paste(ID, Type, sep="_"), feature_value1=value, feature_original1=original)

  gene2_data <- makeTallDataFrame(con,
                                  gene2,
                                  cell_lines,
                                  drugs=NULL,
                                  data_type2) %>%
    transmute(CCLE_name, gene2=ID, feature_type2=Type, feature_name2=paste(ID, Type, sep="_"), feature_value2=value, feature_original2=original)


  cls.df <- src_sqlite(con@dbname) %>%
    tbl('ccle_sampleinfo') %>%
    transmute(CCLE_name, tissue=Site_primary) %>%
    collect

  plot_data <- gene1_data %>%
    inner_join(gene2_data, by='CCLE_name') %>%
    left_join(cls.df, by='CCLE_name')

  return(plot_data)
}
