#' Import Cosmic Cell Lines Project exome sequencing data
#'
#' This function imports the information in CosmicCLP_CompleteExport_v74.tsv into the cosmicclp_exome table in the database.
#' Also indexes the table for fast retrieval.
#'
#' @param fn The path of the data file
#' @param con A \code{SQLiteConnection} object to the database
#' @return TRUE or FALSE depending on whether the data has been written successfully
#' @export
importCosmicCLP_exome <- function ( fn , con ) {

  require(readr)
  #fn <- '~/BigData/CellLineData/RawData/CosmicCLP_CompleteExport_v74.tsv'

  message('Parse the Cosmic CLP exome data file')
  data <- read_tsv(fn,
                   col_names=c('gene_name', 'accession_number', 'hgnc_id', 'sample_name', 'id_sample',
                               'id_tumour', 'mutation_id', 'mutation_cds', 'mutation_aa', 'mutation_description',
                               'mutation_zygosity', 'loh', 'grch', 'mutation_genome_position', 'strand',
                               'snp', 'fathmm_prediction', 'fathmm_score', 'mutation_somatic_status'),
                   col_types='cc_iccc_________cccccccccccdc________',
                   skip=1)

  message('Write the data to the database')
  dbWriteTable(con, "cosmicclp_exome", as.data.frame(data), overwrite=TRUE)

  ##index
  message('Indexing the table')
  dbSendQuery ( con , sprintf(' CREATE INDEX `cosmicclp_exome_gene_name` ON `%s` (`gene_name` ASC); ', 'cosmicclp_exome' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `cosmicclp_exome_sample_name` ON `%s` (`sample_name` ASC); ', 'cosmicclp_exome' )   )
  dbSendQuery ( con , sprintf(' CREATE INDEX `cosmicclp_exome_gene_name_AND_sample_name` ON `%s` (`gene_name`,`sample_name` ASC); ', 'cosmicclp_exome' )   )

  message('Finished importing Cosmic CLP exome sequencing data')

}
