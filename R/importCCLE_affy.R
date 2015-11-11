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
