#functions to export data from sqlite db

#' Extract drug response data
#'
#' This function creates a \code{data.frame} containing the drug response data from the database.
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drugs A vector of compound identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
getDrugData <- function(con, drugs, cell_lines) {

drugs.sql <- paste(drugs, collapse="','")
cell_lines.sql <- paste(cell_lines, collapse="','")
sql <- sprintf("select CCLE_name, Compound as ID, 'resp' as Type, EC50_uM as original, EC50_uM as value
               from ccle_drug_data
               where CCLE_name IN ('%s') and Compound IN ('%s')", cell_lines.sql, drugs.sql)
return(dbGetQuery(con, sql))

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
  sql <- sprintf("select sample_name as CCLE_name, gene_name as ID,  mutation_aa, mutation_description
                 from cosmicclp_exome
                 where sample_name IN ('%s') and gene_name IN ('%s')", cell_lines.sql, genes.sql)
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
  return(dbGetQuery(con, sql))

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
  return(dbGetQuery(con, sql))

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
make_df <- function(con, genes, cell_lines, drugs, data_types=c('affy', 'cn', 'hybcap', 'resp')) {
  require(reshape2)
  require(dplyr)
  affy.df <- getAffyData(con, genes, cell_lines)
  cn.df <- getCopyNumberData(con, genes, cell_lines)
  hybcap.df <- getHybcapData(con, genes, cell_lines)
  drug.df <- getDrugData(con, drugs, cell_lines)
  all.df <- rbind(affy.df, cn.df, hybcap.df, drug.df) %>% filter(Type %in% data_types)
  out <- dcast(all.df , CCLE_name ~ ID + Type, value.var='value' )
  return(out)

}
