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
  data <- DBI::dbGetQuery(con, sql)

  #process variant classifications - only certain types counted as variants
  data <- data %>% dplyr::filter(grepl('Missense|Nonsense|Frame_Shift', Variant_Classification)) %>%
    dplyr::select(-Variant_Classification) %>%
    dplyr::group_by(CCLE_name, ID) %>%
    dplyr::summarise(original=paste(Protein_Change, collapse='|'), value=1) %>%
    dplyr::ungroup()

  #which samples were actually sequenced?
  tested <- DBI::dbGetQuery(con, 'select distinct Tumor_Sample_Barcode as CCLE_name from ccle_hybcap')
  tested <- tested %>% dplyr::filter(CCLE_name %in% cell_lines)
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
  sequenced.df <- sequenced.df %>% dplyr::filter(!( paste(CCLE_name, ID) %in% paste(data$CCLE_name, data$ID) ))

  #combine hybcap dataframes and add additional standard columns
  outdata <- dplyr::bind_rows(data, sequenced.df, notsequenced.df) %>%
    dplyr::transmute(CCLE_name, ID, Type='hybcap', original, value)

  #make sure data types correct
  outdata <- outdata %>% dplyr::mutate_each(dplyr::funs(as.character), -value) %>%
    dplyr::mutate_each(dplyr::funs(as.numeric), value)

  return(outdata)

}
