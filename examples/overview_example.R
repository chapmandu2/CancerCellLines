library(CancerCellLines)
library(dplyr)

#make a toy database
test_db <- makeToyDB()
test_db
test_db@dbname
dbListTables(test_db)

#use normal DBI functions
dbGetQuery(test_db, "select * from ccle_affy limit 10")
dbGetQuery(test_db, "select * from ccle_sampleinfo limit 10")[,1:5]
dbGetQuery(test_db, "select Symbol, t1.CCLE_name, Signal, Site_primary, Hist_subtype1 from ccle_affy as t1
           inner join ccle_sampleinfo t2 on t1.CCLE_name = t2.CCLE_name
           where t2.Hist_subtype1 == 'ductal_carcinoma'
           order by Symbol desc
           limit 10")

#database is indexed
#retrieval syntax becomes inconvenient
dbGetQuery(test_db, "select * from ccle_affy
                      where symbol IN ('PTEN', 'TP53', 'BRAF' ) and
                            CCLE_name IN ('BT474_BREAST', 'MDAMB468_BREAST')
                      limit 10")

symbols <- c('PTEN', 'TP53', 'BRAF')
cell_lines <- c('BT474_BREAST', 'MDAMB468_BREAST')
symbols.sql <- paste(symbols, collapse="','")
cell_lines.sql <- paste(cell_lines, collapse="','")

dbGetQuery(test_db, sprintf("select * from ccle_affy
                      where symbol IN ('%s' ) and
                            CCLE_name IN ('%s')
                      limit 10", symbols.sql, cell_lines.sql))

#dplyr makes things nicer
con <- src_sqlite(test_db@dbname)
ccle_affy <- con %>% tbl('ccle_affy')
ccle_affy
ccle_sampleinfo <- con %>% tbl('ccle_sampleinfo')
ccle_sampleinfo

ccle_sampleinfo %>% dplyr::select(CCLE_name, Site_primary, Hist_subtype1) %>%
  dplyr::filter(Hist_subtype1 == 'ductal_carcinoma') %>%
  dplyr::inner_join(ccle_affy, by='CCLE_name') %>%
  dplyr::arrange(desc(Symbol))

ccle_affy %>% dplyr::filter(symbol %in% symbols & CCLE_name %in% cell_lines)

#convenience functions for different data types are available to make things easier
getAffyData(test_db, symbols, cell_lines)
getCopyNumberData(test_db, symbols, cell_lines)
getHybcapData(test_db, symbols, cell_lines)
getCosmicCLPData(test_db, symbols, cell_lines)

#cell line id's are pain - built in functionality to handle this
con %>% tbl('cell_line_ids') %>% filter(unified_id %in% cell_lines)

#handle drug response data elegantly, from CCLE
drugs <- c('Lapatinib', 'AZD6244', 'Nilotinib' )
getDrugData_CCLE(test_db, drugs, cell_lines)

#and custom data
data(dietlein_data)
head(dietlein_data)
getDrugData_custom(dietlein_data, drugs = 'KU60648_pGI50', cell_lines = c('DMS114_LUNG', 'A549_LUNG'))

#all conveniece functions have a comme table format, so can be combined.
makeTallDataFrame(test_db, symbols, cell_lines, drugs)

#convert to a wide format
my_df <- makeTallDataFrame(test_db, symbols, cell_lines, drugs)
makeWideFromTallDataFrame(my_df)

#all in one function
makeWideDataFrame(test_db, symbols, cell_lines, drugs)

#specify which data types you want to be included
makeWideDataFrame(test_db, symbols, cell_lines, drugs, data_types=c('hybcap', 'affy', 'resp'))

#####
# FULL DATASET
#####

#now let's work with the full dataset - more in the vignette about how to set this up
dbpath <- '~/BigData/CellLineData/CancerCellLines.db'
full_con <- setupSQLite(dbpath)

#show off the speed of the indexed database
dplyr_con <- src_sqlite(full_con@dbname)

#get 2000 random genes
random_genes <- dplyr_con %>% tbl('ccle_affy') %>% group_by(Symbol) %>% summarise(N=n()) %>%
  ungroup() %>% collect %>%
  dplyr::filter(N < mean(N)) %>% sample_n(2000) %>% as.data.frame
random_genes <- random_genes$Symbol

#get 200 random cell lines
random_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% dplyr::select(CCLE_name) %>%
  distinct %>% collect %>% sample_n(200) %>% as.data.frame
random_cell_lines <- random_cell_lines$CCLE_name

#get 10 random compounds
random_drugs <- dplyr_con %>% tbl('ccle_drug_data') %>% dplyr::select(Compound) %>%
  distinct %>% collect %>% sample_n(10) %>% as.data.frame
random_drugs <- random_drugs$Compound

#retrieve the data
test_affy <- getAffyData(full_con, random_genes, random_cell_lines)
test_cn <- getCopyNumberData(full_con, random_genes, random_cell_lines)
test_hybcap <- getHybcapData(full_con, random_genes, random_cell_lines)
test_cosmicclp <- getCosmicCLPData(full_con, random_genes, random_cell_lines)

#make a big data frame
big_df <- makeWideDataFrame(full_con, random_genes, random_cell_lines, random_drugs)

#without resp data
big_df <- makeWideDataFrame(full_con, random_genes, random_cell_lines, drugs=NULL, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'))

#with custom resp data
big_df <- makeWideDataFrame(full_con, random_genes, cell_lines = c('DMS114_LUNG', 'A549_LUNG'), drugs = 'KU60648_pGI50', drug_df = dietlein_data)

