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
