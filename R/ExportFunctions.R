#functions to export data from sqlite db

getDrugData <- function(con, drugs, cell_lines) {

drugs.sql <- paste(drugs, collapse="','")
cell_lines.sql <- paste(cell_lines, collapse="','")
sql <- sprintf("select CCLE_name, Compound as ID, 'resp' as Type, EC50_uM as original, EC50_uM as value
               from ccle_drug_data
               where CCLE_name IN ('%s') and Compound IN ('%s')", cell_lines.sql, drugs.sql)
return(dbGetQuery(con, sql))

}

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

getAffyData <- function(con, genes, cell_lines) {

  genes.sql <- paste(genes, collapse="','")
  cell_lines.sql <- paste(cell_lines, collapse="','")
  sql <- sprintf("select CCLE_name, Symbol as ID, 'affy' as Type, Signal as original, Signal as value
               from ccle_affy
               where CCLE_name IN ('%s') and Symbol IN ('%s')", cell_lines.sql, genes.sql)
  return(dbGetQuery(con, sql))

}

make_df <- function(con, genes, cell_lines, drugs, data_types=c('affy', 'hybcap', 'resp')) {
  require(reshape2)
  require(dplyr)
  affy.df <- getAffyData(con, genes, cell_lines)
  hybcap.df <- getHybcapData(con, genes, cell_lines)
  drug.df <- getDrugData(con, drugs, cell_lines)
  all.df <- rbind(affy.df, hybcap.df, drug.df) %>% filter(Type %in% data_types)
  out <- dcast(all.df , CCLE_name ~ ID + Type, value.var='value' )
  return(out)

}
