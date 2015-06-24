#functions to import raw text files to a SQLite database

#setupSQLite <- function ( dbpath=system.file('extdata/toy.db', package="CancerCellLines") ) {
setupSQLite <- function ( dbpath ) {
  require(RSQLite)
  #dbpath <- '~/BigData/CellLineData/CancerCellLines.db'
  dbpath <- system.file('extdata/toy.db', package="CancerCellLines")
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = dbpath)
  return(con)
}

importCCLE_info <- function ( fn , con )  {

  #fn <- '~/BigData/CellLineData/RawData/CCLE_sample_info_file_2012-10-18.txt'

  data <- read.table(fn,header=T, sep='\t')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Cell_line_aliases', 'Gender', 'Site_primary', 'Histology', 'Hist_subtype1', 'Notes', 'Source', 'Expression_arrays', 'SNP_arrays', 'Oncomap', 'Hybrid_capture_sequencing')

  dbWriteTable(con, "ccle_sampleinfo", data, overwrite=TRUE)


}

importCCLE_affy <- function ( fn , con ) {
  require(reshape2)
  #fn <- '~/BigData/CellLineData/RawData/CCLE_Expression_Entrez_2012-09-29.gct'

  data <- read.table(fn,header=T, sep='\t', skip=2)

  ###unpivot the expression data
  data_tall <- melt(data, id.vars =  colnames(data)[1:2])
  colnames(data_tall) <- c('ProbeID', 'Symbol', 'CCLE_name', 'Signal')

  #write to db
  dbWriteTable(con, "ccle_affy", data_tall, overwrite=TRUE)

  ##index
  dbGetQuery ( con , sprintf(' CREATE INDEX `ccle_affy_symbol` ON `%s` (`Symbol` ASC); ', 'ccle_affy' )   )
  dbGetQuery ( con , sprintf(' CREATE INDEX `ccle_affy_CCLE_name` ON `%s` (`CCLE_name` ASC); ', 'ccle_affy' )   )


}

importCCLE_hybcap <- function ( fn , con ) {

  require(readxl)
  #fn <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx'

  data <- read_excel(fn, col_types=c('text', 'numeric', 'text', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                                     'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text'))

  dbWriteTable(con, "ccle_hybcap", as.data.frame(data), overwrite=TRUE)

}

importCCLE_drugresponse <- function (fn , con) {
  #fn <- '~/BigData/CellLineData/RawData/CCLE_NP24.2009_Drug_data_2012.02.20.csv'

  data <- read.table(fn, header=T, sep=',')
  colnames(data) <- c('CCLE_name', 'Primary_cell_name', 'Compound', 'Target', 'Doses_uM', 'Activity_median', 'Activity_sd', 'Fit_type', 'Num_data', 'EC50_uM', 'IC50_uM', 'Amax', 'Act_area')

  dbWriteTable(con, "ccle_drug_data", data, overwrite=TRUE)
}

makeToyDB <- function() {

  toy_db <- tempfile()
  toy_con <- setupSQLite(toy_db)
  importCCLE_info(system.file("extdata", "CCLE_sample_info_file_2012-10-18_toy.txt", package = "CancerCellLines") , toy_con)
  importCCLE_affy(system.file("extdata", "CCLE_Expression_Entrez_2012-09-29_toy.gct", package = "CancerCellLines") , toy_con)
  importCCLE_hybcap(system.file("extdata", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_toy.xlsx", package = "CancerCellLines") , toy_con)
  importCCLE_drugresponse(system.file("extdata", "CCLE_NP24.2009_Drug_data_2012.02.20_toy.txt", package = "CancerCellLines") , toy_con)

  print('A database of toy data has successfully been created and a connector returned')
  return(toy_con)

}
