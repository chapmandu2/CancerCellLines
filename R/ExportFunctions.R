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
sql <- sprintf("select CCLE_name, Compound as ID, 'resp' as Type, EC50_uM as original
               from ccle_drug_data
               where CCLE_name IN ('%s') and Compound IN ('%s')", cell_lines.sql, drugs.sql)
out <- dbGetQuery(con, sql)
out <- out %>% mutate(value=6-log10(original)) %>%
  mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)
return(out)

}

#' Extract custom drug response data
#'
#' This function creates a \code{data.frame} containing the drug response data from a user provided data frame
#'
#' @param df A \code{data.frame} object containing the data which should include columns names compound_id, unified_id, endpoint, original and value
#' @param drugs A vector of compound identifiers or compound/endpoint identifiers
#' @param cell_lines A vector of cell line identifiers
#' @return A \code{data.frame} containing the drug response data for the requested compounds and cell lines
#' @export
getDrugData_custom <- function(df, drugs, cell_lines) {

  required_colnames <- c('compound_id', 'unified_id', 'endpoint', 'original', 'value')
  missing_columns <- setdiff(required_colnames,colnames(df))

  if(length(missing_columns) > 0) {
    stop(sprintf('Following required column names not present in supplied data frame: %s', sprintf(paste(missing_columns, collapse=', '))))
    }

  outdata <- df %>% dplyr::transmute(CCLE_name=unified_id,
                                     compound_id,
                                     ID=paste(compound_id, endpoint, sep='_'),
                                     Type='resp',
                                     original,
                                     value) %>%
                    dplyr::filter((compound_id %in% drugs | ID %in% drugs ) & CCLE_name %in% cell_lines) %>%
                    dplyr::select(-compound_id) %>% as.data.frame

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

  #make sql
  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select Tumor_Sample_Barcode as CCLE_name, Hugo_Symbol as ID,  Protein_Change, Variant_Classification
                 from ccle_hybcap
                 where CCLE_name IN ('%s') and Hugo_Symbol IN ('%s')", cell_lines.sql, genes.sql)

  #get data
  data <- dbGetQuery(con, sql)

  #process variant classifications - only certain types counted as variants
  data <- data %>% filter(grepl('Missense|Nonsense|Frame_Shift', Variant_Classification)) %>%
      select(-Variant_Classification) %>% group_by(CCLE_name, ID) %>%
    summarise(original=paste(Protein_Change, collapse='|'), value=1) %>% ungroup()

  #which samples were actually sequenced?
  tested <- dbGetQuery(con, 'select distinct Tumor_Sample_Barcode as CCLE_name from ccle_hybcap')
  tested <- tested %>% filter(CCLE_name %in% cell_lines)
  sequenced_ids <- tested$CCLE_name
  notsequenced_ids <- setdiff(cell_lines, sequenced_ids)

  #generate dataframes for sequenced and not sequenced cell lines
  sequenced.df <- data.frame ( CCLE_name=rep(sequenced_ids, length(genes)) ,
                                   ID=rep(genes, each= length(sequenced_ids) ),
                                   original='-',
                                   value=0, stringsAsFactors=FALSE)

  if (length(notsequenced_ids) > 0 ) {

    notsequenced.df <- data.frame ( CCLE_name=rep(notsequenced_ids, length(genes)) ,
                                        ID=rep(genes, each= length(notsequenced_ids) ),
                                        original=NA,
                                        value=NA, stringsAsFactors=FALSE )

  } else {
    notsequenced.df <- data.frame ()
  }

  #get rid of rows in sequenced.df which are duplicated in data
  sequenced.df <- sequenced.df %>% filter(!( paste(CCLE_name, ID) %in% paste(data$CCLE_name, data$ID) ))

  #combine hybcap dataframes and add additional standard columns
  outdata <- bind_rows(data, sequenced.df, notsequenced.df) %>%
    transmute(CCLE_name, ID, Type='hybcap', original, value)

  #make sure data types correct
  outdata <- outdata %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(outdata)

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

  #make sql - process cell line ids into unified ids
  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select t2.unified_id as CCLE_name, t1.gene_name as ID,  t1.mutation_aa, t1.mutation_description
                 from cosmicclp_exome t1
                  inner join cell_line_ids t2 on t1.sample_name = t2.native_id
                 where t2.id_type = 'cosmic_clp' and
                 t2.unified_id IN ('%s') and gene_name IN ('%s')", cell_lines.sql, genes.sql)

  #get data
  data <- dbGetQuery(con, sql)

  #process variant classifications - only certain types counted as variants
  data <- data %>% filter(grepl('Missense|Nonsense|Frameshift', mutation_description)) %>%
    select(-mutation_description) %>% group_by(CCLE_name, ID) %>%
    summarise(original=paste(mutation_aa, collapse='|'), value=1) %>% ungroup()

  #which samples were actually sequenced?
  tested <- dbGetQuery(con, "select distinct t2.unified_id as CCLE_name
                       from cosmicclp_exome t1
                       inner join cell_line_ids t2 on t1.sample_name = t2.native_id
                       where t2.id_type = 'cosmic_clp'")
  tested <- tested %>% filter(CCLE_name %in% cell_lines)
  sequenced_ids <- tested$CCLE_name
  notsequenced_ids <- setdiff(cell_lines, sequenced_ids)

  #generate dataframes for sequenced and not sequenced cell lines
  sequenced.df <- data.frame ( CCLE_name=rep(sequenced_ids, length(genes)) ,
                               ID=rep(genes, each= length(sequenced_ids) ),
                               original='-',
                               value=0, stringsAsFactors=FALSE)

  if (length(notsequenced_ids) > 0 ) {

    notsequenced.df <- data.frame ( CCLE_name=rep(notsequenced_ids, length(genes)) ,
                                    ID=rep(genes, each= length(notsequenced_ids) ),
                                    original=NA,
                                    value=NA, stringsAsFactors=FALSE )

  } else {
    notsequenced.df <- data.frame ()
  }

  #get rid of rows in sequenced.df which are duplicated in data
  sequenced.df <- sequenced.df %>% filter(!( paste(CCLE_name, ID) %in% paste(data$CCLE_name, data$ID) ))

  #combine hybcap dataframes and add additional standard columns
  outdata <- bind_rows(data, sequenced.df, notsequenced.df) %>%
    transmute(CCLE_name, ID, Type='cosmicclp', original, value)

  #make sure data types correct
  outdata <- outdata %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(outdata)

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

  #duplicate remove
  data <- data %>% group_by(ID) %>% mutate(N=n()) %>% ungroup %>% filter(N==median(N)) %>% select(-N) %>% as.data.frame

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

  #duplicate remove
  data <- data %>% group_by(ID) %>% mutate(N=n()) %>% ungroup %>% filter(N==median(N)) %>% select(-N) %>% as.data.frame

  data <- data %>% mutate_each(funs(as.character), -value) %>% mutate_each(funs(as.numeric), value)

  return(data)

}
