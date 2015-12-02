#' Set up a toy database
#'
#' This function imports the example toy dataset included with the package into a temporary database.
#'
#' @return A \code{SQLiteConnection} object
#' @export
makeToyDB <- function() {

  toy_db <- tempfile()
  toy_con <- setupSQLite(toy_db)
  importCCLE_info(system.file("extdata", "CCLE_sample_info_file_2012-10-18_toy.txt", package = "CancerCellLines") , toy_con)
  importCCLE_affy(system.file("extdata", "CCLE_Expression_Entrez_2012-09-29_toy.gct", package = "CancerCellLines") , toy_con)
  importCCLE_cn(system.file("extdata", "CCLE_copynumber_byGene_2012-09-29_toy.txt", package = "CancerCellLines") , toy_con)
  importCCLE_hybcap(system.file("extdata", "CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_toy.maf", package = "CancerCellLines") , toy_con)
  importCosmicCLP_exome(system.file("extdata", "CosmicCLP_CompleteExport_v74_toy.tsv", package = "CancerCellLines") , toy_con)
  importCCLE_drugresponse(system.file("extdata", "CCLE_NP24.2009_Drug_data_2012.02.20_toy.txt", package = "CancerCellLines") , toy_con)
  importCellLineIDs(system.file("extdata", "CellLineIDNormalisationNov15_toy.txt", package = "CancerCellLines") , toy_con)

  print('A database of toy data has successfully been created and a connector returned')
  return(toy_con)

}
