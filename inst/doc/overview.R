## ---- echo=FALSE, message=FALSE------------------------------------------
library(CancerCellLines)

## ------------------------------------------------------------------------
list.files(system.file("extdata", package = "CancerCellLines"))

## ----eval=TRUE-----------------------------------------------------------
test_db <- makeToyDB()
test_db
test_db@dbname
dbListTables(test_db)

## ----eval=TRUE-----------------------------------------------------------
test_db <- setupSQLite(system.file('extdata/toy.db', package="CancerCellLines"))
test_db
test_db@dbname
dbListTables(test_db)

## ------------------------------------------------------------------------
dbGetQuery(test_db, "select * from ccle_affy limit 10")
dbGetQuery(test_db, "select * from ccle_sampleinfo limit 10")[,1:5]
dbGetQuery(test_db, "select Symbol, t1.CCLE_name, Signal, Site_primary, Hist_subtype1 from ccle_affy as t1 
                      inner join ccle_sampleinfo t2 on t1.CCLE_name = t2.CCLE_name
                      where t2.Hist_subtype1 == 'ductal_carcinoma'
                      order by Symbol desc
                      limit 10")

## ------------------------------------------------------------------------
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


## ------------------------------------------------------------------------
con <- src_sqlite(test_db@dbname) 
ccle_affy <- con %>% tbl('ccle_affy')
ccle_affy
ccle_sampleinfo <- con %>% tbl('ccle_sampleinfo')
ccle_sampleinfo

ccle_sampleinfo %>% dplyr::select(CCLE_name, Site_primary, Hist_subtype1) %>% 
  dplyr::filter(Hist_subtype1 == 'ductal_carcinoma') %>%
  dplyr::inner_join(ccle_affy, by='CCLE_name') %>%
  dplyr::arrange(desc(Symbol))

ccle_affy %>% filter(symbol %in% symbols & CCLE_name %in% cell_lines)


## ------------------------------------------------------------------------
getAffyData(test_db, symbols, cell_lines)
getCopyNumberData(test_db, symbols, cell_lines)

## ------------------------------------------------------------------------
getHybcapData(test_db, symbols, cell_lines)
getCosmicCLPData(test_db, symbols, cell_lines)

## ------------------------------------------------------------------------
con %>% tbl('cell_line_ids') %>% filter(unified_id %in% cell_lines)

## ------------------------------------------------------------------------
drugs <- c('Lapatinib', 'AZD6244', 'Nilotinib' )
getDrugData_CCLE(test_db, drugs, cell_lines)

## ------------------------------------------------------------------------
data(dietlein_data)
head(dietlein_data)
getDrugData_custom(dietlein_data, drugs = 'KU60648_pGI50', cell_lines = c('DMS114_LUNG', 'A549_LUNG'))


## ------------------------------------------------------------------------
makeTallDataFrame(test_db, symbols, cell_lines, drugs)

## ------------------------------------------------------------------------
my_df <- makeTallDataFrame(test_db, symbols, cell_lines, drugs)
makeWideFromTallDataFrame(my_df)


## ------------------------------------------------------------------------
makeWideDataFrame(test_db, symbols, cell_lines, drugs)

## ------------------------------------------------------------------------
makeWideDataFrame(test_db, symbols, cell_lines, drugs, data_types=c('hybcap', 'affy', 'resp'))


## ----eval=FALSE----------------------------------------------------------
#  dbpath <- '~/BigData/CellLineData/CancerCellLines.db'
#  infopath <- '~/BigData/CellLineData/RawData/CCLE_sample_info_file_2012-10-18.txt'
#  affypath <- '~/BigData/CellLineData/RawData/CCLE_Expression_Entrez_2012-09-29.gct'
#  cnpath <- '~/BigData/CellLineData/RawData/CCLE_copynumber_byGene_2012-09-29.txt'
#  hybcappath <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf'
#  cosmicclppath <- '~/BigData/CellLineData/RawData/CosmicCLP_CompleteExport_v74.tsv'
#  drugpath <- '~/BigData/CellLineData/RawData/CCLE_NP24.2009_Drug_data_2012.02.20.csv'
#  idspath <- system.file("extdata", "CellLineIDNormalisationNov15.txt", package = "CancerCellLines")
#  
#  

## ----eval=FALSE----------------------------------------------------------
#    full_con <- setupSQLite(dbpath)
#    importCCLE_info(infopath , full_con)
#    importCCLE_hybcap(hybcappath , full_con)
#    importCosmicCLP_exome(cosmicclppath, full_con)
#    importCCLE_drugresponse(drugpath , full_con)
#    importCCLE_affy(affypath , full_con)
#    importCCLE_cn(cnpath, full_con)
#    importCellLineIDs(idspath, full_con)
#  

## ----eval=FALSE----------------------------------------------------------
#      dplyr_con <- src_sqlite(full_con@dbname)
#  
#      #get 2000 random genes
#      random_genes <- dplyr_con %>% tbl('ccle_affy') %>% group_by(Symbol) %>% summarise(N=n()) %>%
#        ungroup() %>% collect %>%
#        dplyr::filter(N < mean(N)) %>% sample_n(2000) %>% as.data.frame
#      random_genes <- random_genes$Symbol
#  
#      #get 200 random cell lines
#      random_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% dplyr::select(CCLE_name) %>%
#        distinct %>% collect %>% sample_n(200) %>% as.data.frame
#      random_cell_lines <- random_cell_lines$CCLE_name
#  
#      #get 10 random compounds
#      random_drugs <- dplyr_con %>% tbl('ccle_drug_data') %>% dplyr::select(Compound) %>%
#        distinct %>% collect %>% sample_n(10) %>% as.data.frame
#      random_drugs <- random_drugs$Compound
#  
#      #retrieve the data
#      test_affy <- getAffyData(full_con, random_genes, random_cell_lines)
#      test_cn <- getCopyNumberData(full_con, random_genes, random_cell_lines)
#      test_hybcap <- getHybcapData(full_con, random_genes, random_cell_lines)
#      test_cosmicclp <- getCosmicCLPData(full_con, random_genes, random_cell_lines)
#  
#      #make a big data frame
#      big_df <- makeWideDataFrame(full_con, random_genes, random_cell_lines, random_drugs)
#  
#      #without resp data
#      big_df <- makeWideDataFrame(full_con, random_genes, random_cell_lines, drugs=NULL, data_types=c('affy', 'cn', 'hybcap', 'cosmicclp'))
#  
#      #with custom resp data
#      big_df <- makeWideDataFrame(full_con, random_genes, cell_lines = c('DMS114_LUNG', 'A549_LUNG'), drugs = 'KU60648_pGI50', drug_df = dietlein_data)
#  
#  

## ------------------------------------------------------------------------
   sessionInfo() 

