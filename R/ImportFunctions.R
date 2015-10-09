#functions to import raw text files to a SQLite database

#' Connect to or create a RSQlite connection
#'
#' This function connects to or creates if necessary a SQLite database and returns a SQLiteConnection object.
#'
#' @param dbpath The path of the database
#' @return A \code{SQLiteConnection} object
#' @export
setupSQLite <- function ( dbpath=system.file('extdata/toy.db', package="CancerCellLines") ) {
  require(RSQLite)
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = dbpath)
  return(con)
}

#' Import CCLE cell line information
#'
#' This function imports the information in CCLE_sample_info_file_2012-10-18.txt into the database.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_info <- function ( fn , con )  {

  data <- read.table(fn,header=T, sep='\t')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Cell_line_aliases', 'Gender', 'Site_primary', 'Histology', 'Hist_subtype1', 'Notes', 'Source', 'Expression_arrays', 'SNP_arrays', 'Oncomap', 'Hybrid_capture_sequencing')

  dbWriteTable(con, "ccle_sampleinfo", data, overwrite=TRUE)

}

#' Import CCLE affy data
#'
#' This function imports the information in CCLE_Expression_Entrez_2012-09-29.gct into the ccle_affy table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_affy <- function ( fn , con ) {
  require(reshape2)

  message('Parse the gene expression data file')
  data <- read.table(fn,header=T, sep='\t', skip=2)

  ###unpivot the expression data
  data_tall <- melt(data, id.vars =  colnames(data)[1:2])
  colnames(data_tall) <- c('ProbeID', 'Symbol', 'CCLE_name', 'Signal')

  #write to db
  message('Writing to database')
  dbWriteTable(con, "ccle_affy", data_tall, overwrite=TRUE)

  ##index
  message('Indexing the table')
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_affy_symbol` ON `%s` (`Symbol` ASC); ', 'ccle_affy' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_affy_CCLE_name` ON `%s` (`CCLE_name` ASC); ', 'ccle_affy' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_affy_symbol_AND_CCLE_name` ON `%s` (`Symbol`,`CCLE_name` ASC); ', 'ccle_affy' )   )

  message('Finished importing affy data')

}

#' Import CCLE copy number data
#'
#' This function imports the information in CCLE_copynumber_byGene_2012-09-29.txt into the ccle_cn table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_cn <- function ( fn , con ) {
  require(reshape2)

  message('Parse the copy number data file')
  data <- read.table(fn,header=T, sep='\t')

  ###unpivot the copy number data
  data_tall <- melt(data, id.vars =  colnames(data)[1:4])
  colnames(data_tall) <- c('Symbol', 'Chr', 'txStart', 'txEnd', 'CCLE_name', 'log2cn')

  #write to db
  message('Writing to database')
  dbWriteTable(con, "ccle_cn", data_tall, overwrite=TRUE)

  ##index
  message('Indexing the table')
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_cn_symbol` ON `%s` (`Symbol` ASC); ', 'ccle_cn' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_cn_CCLE_name` ON `%s` (`CCLE_name` ASC); ', 'ccle_cn' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `ccle_cn_symbol_AND_CCLE_name` ON `%s` (`Symbol`,`CCLE_name` ASC); ', 'ccle_cn' )   )

  message('Finished importing cn data')

}

#' Import CCLE hybrid capture sequencing data
#'
#' This function imports the information in CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx into the ccle_hybap table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_hybcap <- function ( fn , con ) {

  require(readxl)
  #fn <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx'

  data <- read_excel(fn, col_types=c('text', 'numeric', 'text', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                                     'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text'))

  dbWriteTable(con, "ccle_hybcap", as.data.frame(data), overwrite=TRUE)

}

#' Import CCLE drug response data
#'
#' This function imports the information in CCLE_NP24.2009_Drug_data_2012.02.20.csv into the ccle_drugresponse table in the database.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCCLE_drugresponse <- function (fn , con) {

  data <- read.table(fn, header=T, sep=',')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Compound', 'Target', 'Doses_uM', 'Activity_median', 'Activity_sd', 'Fit_type', 'Num_data', 'EC50_uM', 'IC50_uM', 'Amax', 'Act_area')

  dbWriteTable(con, "ccle_drug_data", data, overwrite=TRUE)
}

#' Set up a toy database
#'
#' This function imports the example toy dataset included with the package into a temporary database.
#'
#' @param none
#' @return A \code{SQLiteConnection} object
#' @export
makeToyDB <- function() {

  toy_db <- tempfile()
  toy_con <- setupSQLite(toy_db)
  importCCLE_info(system.file("extdata", "CCLE_sample_info_file_2012-10-18_toy.txt", package = "CancerCellLines") , toy_con)
  importCCLE_affy(system.file("extdata", "CCLE_Expression_Entrez_2012-09-29_toy.gct", package = "CancerCellLines") , toy_con)
  importCCLE_cn(system.file("extdata", "CCLE_copynumber_byGene_2012-09-29_toy.txt", package = "CancerCellLines") , toy_con)
  importCCLE_hybcap(system.file("extdata", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_toy.xlsx", package = "CancerCellLines") , toy_con)
  importCCLE_drugresponse(system.file("extdata", "CCLE_NP24.2009_Drug_data_2012.02.20_toy.txt", package = "CancerCellLines") , toy_con)

  print('A database of toy data has successfully been created and a connector returned')
  return(toy_con)

}
