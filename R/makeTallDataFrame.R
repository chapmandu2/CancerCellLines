#' Merge different data types into a single tall data frame
#'
#' This function creates a \code{data.frame} containing data for different data types in a form suitable for further visualisation in R.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vector of cell line identifiers
#' @param drugs A vector of compound identifiers
#' @param data_types A vector with default \code{c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp')} to specify which data types should be returned.
#' @param drug_df A \code{data.frame} containing input for the \code{getDrugData_custom} function
#' @return A tall \code{data.frame} containing the requested data
#' @export
makeTallDataFrame <- function(con, genes, cell_lines, drugs, data_types=c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp'), drug_df=NULL) {

  #make a data container list
  all_data.list <- list()

  #get the data
  if ('affy' %in% data_types) {
    all_data.list[['affy_data']] <- getAffyData(con, genes, cell_lines)
  }

  if ('hybcap' %in% data_types) {
    all_data.list[['hybcap_data']] <- getHybcapData(con, genes, cell_lines)
  }

  if ('cosmicclp' %in% data_types) {
    all_data.list[['cosmicclp_data']] <- getCosmicCLPData(con, genes, cell_lines)
  }

  if ('cn' %in% data_types) {
    all_data.list[['cn_data']] <- getCopyNumberData(con, genes, cell_lines)
  }

  if ('resp' %in% data_types) {
    if(is.null(drug_df)) {
      all_data.list[['resp_data']] <- getResponseData(con, drugs, cell_lines, resp_type='ccle')
    } else {
      all_data.list[['resp_data']] <- getResponseData(drug_df, drugs, cell_lines, resp_type='custom')
    }
    #check that all compounds in drugs have been returned, give warning if not
    missing_resps <- setdiff(drugs,unique(all_data.list[['resp_data']]$ID))
    if (length(missing_resps) != 0) {
      warning(sprintf('No data returned for: %s', paste(missing_resps, collapse=', ')))
    }

    #combine
    all_data <- dplyr::bind_rows(all_data.list)

    #get rid of cell lines with no response data
    present_cls <- unique(all_data.list[['resp_data']]$CCLE_name)
    missing_cls <- setdiff(cell_lines, present_cls)
    all_data <- all_data %>% filter(CCLE_name %in% present_cls)
    if (length(missing_cls) != 0) {
      warning(sprintf('No response data for following cell lines: %s', paste(missing_cls, collapse=', ')))
    }
  } else {
    #combine
    all_data <- dplyr::bind_rows(all_data.list)
  }

  return(all_data)

}
