#functions to export data from sqlite db

#' Extract CCLE drug response data
#'
#' This function creates a \code{data.frame} containing the CCLE drug response data from the database.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
getDrugData_CCLE <- function(con, drugs, cell_lines) {

drugs.sql <- paste(drugs, collapse="','")
cell_lines.sql <- paste(cell_lines, collapse="','")
sql <- sprintf("select CCLE_name, Compound as ID, 'resp' as Type, EC50_uM as original, EC50_uM as value
               from ccle_drug_data
               where CCLE_name IN ('%s') and Compound IN ('%s')", cell_lines.sql, drugs.sql)
out <- dbGetQuery(con, sql)
out <- out %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)
return(out)

}

#' Extract custom drug response data
#'
#' This function creates a \code{data.frame} containing the drug response data from a user provided data frame
#'
#' @param df A \code{data.frame} object containing the data which should include columns names compound_id, unified_id, endpoint, original and value
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
getDrugData_custom <- function(df, drugs, cell_lines) {

  required_colnames <- c('compound_id', 'unified_id', 'endpoint', 'original', 'value')
  missing_columns <- setdiff(required_colnames,colnames(df))

  if(length(missing_columns) > 0) {
    stop(sprintf('Following required column names not present in supplied data frame: %s', sprintf(paste(missing_columns, collapse=', '))))
    }

  outdata <- df %>% dplyr::filter(compound_id %in% drugs & unified_id %in% cell_lines) %>%
                    dplyr::transmute(CCLE_name=unified_id,
                                     ID=paste(compound_id, endpoint, sep='_'),
                                     Type='resp',
                                     original=original,
                                     value) %>% as.data.frame

  outdata <- outdata %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  if(nrow(outdata)==0) {
    stop('No data found matching supplied drugs/cell lines')
  } else {
  return(outdata)
  }
}

#' Generic function to get response data from multiple sources
#'
#' This function creates a \code{data.frame} containing the response data from a variety of sourses.
#'
#' @param src A \code{SQLiteConnection} object to the database (ccle) or \code{data.frame} object (custom)
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @param resp_type The type of response data: \code{ccle, custom}
#' @return A \code{data.frame} containing the response data for the requested compounds and cell lines
#' @export
getResponseData <- function(src, drugs, cell_lines, resp_type='ccle') {

  if (class(src) == 'SQLiteConnection' & resp_type == 'ccle') {
    getDrugData_CCLE(src, drugs, cell_lines)
  } else if (class(src) == 'data.frame' & resp_type == 'custom') {
    getDrugData_custom(src, drugs, cell_lines)
  } else {
    stop('Check ?getResponseData to ensure you are passing the correct options into src')
  }

}

#' Extract hybrid capture sequencing data
#'
#' This function creates a \code{data.frame} containing the hybrid capture sequencing from the database for the requested cell lines and genes.  Gene/sample pairs where there exists at least one missense, nonsense or framshift mutation are categorised as 1, those that don't are categorised as 0.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the hybrid capture sequencing data for the requested compounds and cell lines
#' @export
getHybcapData <- function(con, genes, cell_lines) {
  require(dplyr)

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select Tumor_Sample_Barcode as CCLE_name, Hugo_Symbol as ID,  Protein_Change, Variant_Classification
                 from ccle_hybcap
                 where CCLE_name IN ('%s') and Hugo_Symbol IN ('%s')", cell_lines.sql, genes.sql)
  data <- dbGetQuery(con, sql)

  data <- data %>% filter(grepl('Missense|Nonsense|Frame_Shift', Variant_Classification)) %>%
      select(-Variant_Classification) %>% group_by(CCLE_name, ID) %>%
    summarise(original=paste(Protein_Change, collapse='|'), value=1) %>% ungroup()

  blank_data <- merge(data.frame(CCLE_name = cell_lines, stringsAsFactors = FALSE),
                      data.frame(ID = genes, stringsAsFactors = FALSE))

  data <- blank_data %>% left_join(data, by=c('CCLE_name', 'ID')) %>%
    mutate(original=ifelse(is.na(original), 'wt', original),
           value = ifelse(is.na(value), 0, value),
           Type = 'hybcap') %>%
    select (CCLE_name, ID, Type, original, value)

  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)

}

#' Extract Cosmic CLP Exome sequencing data
#'
#' This function creates a \code{data.frame} containing the Cosmic CLP exome sequencing datafrom the database for the requested cell lines and genes.  Gene/sample pairs where there exists at least one missense, nonsense or framshift mutation are categorised as 1, those that don't are categorised as 0.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the hybrid capture sequencing data for the requested compounds and cell lines
#' @export
getCosmicCLPData <- function(con, genes, cell_lines) {
  require(dplyr)

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select t2.unified_id as CCLE_name, t1.gene_name as ID,  t1.mutation_aa, t1.mutation_description
                 from cosmicclp_exome t1
                  inner join cell_line_ids t2 on t1.sample_name = t2.native_id
                 where t2.id_type = 'cosmic_clp' and
                 t2.unified_id IN ('%s') and gene_name IN ('%s')", cell_lines.sql, genes.sql)
  data <- dbGetQuery(con, sql)

  data <- data %>% filter(grepl('Missense|Nonsense|Frameshift', mutation_description)) %>%
    select(-mutation_description) %>% group_by(CCLE_name, ID) %>%
    summarise(original=paste(mutation_aa, collapse='|'), value=1) %>% ungroup()

  blank_data <- merge(data.frame(CCLE_name = cell_lines, stringsAsFactors = FALSE),
                      data.frame(ID = genes, stringsAsFactors = FALSE))

  data <- blank_data %>% left_join(data, by=c('CCLE_name', 'ID')) %>%
    mutate(original=ifelse(is.na(original), 'wt', original),
           value = ifelse(is.na(value), 0, value),
           Type = 'cosmicclp') %>%
    select (CCLE_name, ID, Type, original, value)

  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)


}

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
  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)

}

#' Extract copy number data
#'
#' This function creates a \code{data.frame} containing the by gene copy number data from the database for the requested cell lines and genes.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vectore of cell line identifiers
#' @return A \code{data.frame} containing the by gene copy number data for the requested compounds and cell lines
#' @export
getCopyNumberData <- function(con, genes, cell_lines) {

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select CCLE_name, Symbol as ID, 'cn' as Type, log2cn as original, log2cn as value
               from ccle_cn
               where CCLE_name IN ('%s') and Symbol IN ('%s')", cell_lines.sql, genes.sql)

  data <- dbGetQuery(con, sql)
  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)

}

#' Merge different data types into a single data frame
#'
#' This function creates a \code{data.frame} containing data for different data types in a form suitable for further statistical modelling in R.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param genes A vector of gene symbols
#' @param cell_lines A vectore of cell line identifiers
#' @param drugs A vector of compound identifiers
#' @param data_types A vector with default \code{c('affy', 'hybcap', 'resp')} to specify which data types should be returned.
#' @return A \code{data.frame} containing the affymetrix gene expression data for the requested compounds and cell lines
#' @export
make_df <- function(con, genes, cell_lines, drugs, data_types=c('affy', 'cn', 'hybcap', 'resp', 'cosmicclp'), drug_df=NULL) {
  require(tidyr)
  require(dplyr)

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
    #check that all compounds in params have been returned, give warning if not
    missing_resps <- setdiff(drugs,unique(all_data.list[['resp_data']]$ID))
    if (length(missing_resps) != 0) {
      warning(sprintf('No data returned for: %s', paste(missing_resps, collapse=', ')))
    }

    #combine
    all_data <- bind_rows(all_data.list)

    #get rid of cell lines with no response data
    present_cls <- unique(all_data.list[['resp_data']]$CCLE_name)
    missing_cls <- setdiff(cell_lines, present_cls)
    all_data <- all_data %>% filter(CCLE_name %in% present_cls)
    if (length(missing_cls) != 0) {
      warning(sprintf('No response data for following cell lines: %s', paste(missing_cls, collapse=', ')))
    }
  } else {
    #combine
    all_data <- bind_rows(all_data.list)
  }

  #create matrix
  output.df <- all_data %>% transmute(CCLE_name, fn=paste(ID,Type,sep='_'), value) %>% spread(fn, value)

  #reorder columns
  output.df <- output.df %>% dplyr::select(CCLE_name, ends_with('_resp'), everything())

  #make mutation fields into factors
  output.df <- output.df %>% mutate_each(funs(as.factor), ends_with('_hybcap|_cosmicclp'))

  return(output.df)


}
