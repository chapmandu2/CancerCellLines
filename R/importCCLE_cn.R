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
