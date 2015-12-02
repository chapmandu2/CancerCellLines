## ---- echo=FALSE, message=FALSE------------------------------------------
library(CancerCellLines)

## ------------------------------------------------------------------------
dbpath <- '~/BigData/CellLineData/CancerCellLines.db'
#dbpath <- system.file('extdata/toy.db', package="CancerCellLines")
full_con <- setupSQLite(dbpath)
dplyr_con <- src_sqlite(full_con@dbname)

## ------------------------------------------------------------------------
    #specify the genes
    ex1_genes <- c('BRAF', 'NRAS', 'CRAF', 'TP53')
  
    #get the melanoma cell lines
    ex1_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% dplyr::filter(Site_primary=='skin') %>%
       collect %>% as.data.frame
    ex1_cell_lines <- ex1_cell_lines$CCLE_name
    ex1_cell_lines[1:10]
    
    #get BRAF and MEK inhibitors
    ex1_drugs <- c('AZD6244','PLX4720','PD-0325901')
    

## ----fig.width=6, fig.height=6-------------------------------------------
    #make a tall frame
    ex1_tall_df <- makeTallDataFrame(full_con, ex1_genes, ex1_cell_lines, ex1_drugs)
    ex1_tall_df
    
    #convert this into a wide data frame
    ex1_wide_df <- ex1_tall_df %>% makeWideFromTallDataFrame
    ex1_wide_df
    
    #compare the drug activities
    pairs(~AZD6244_resp+PLX4720_resp+`PD-0325901_resp`, ex1_wide_df)
    

## ----fig.width=6, fig.height=6-------------------------------------------
    #make a heatmap!
    plotHeatmap(ex1_tall_df)
    

## ----fig.width=6, fig.height=6-------------------------------------------
    plotHeatmap(ex1_tall_df, order_feature='PLX4720_resp')

## ----fig.width=6, fig.height=4-------------------------------------------
   
    #get all cell lines
    ex2_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% 
       collect %>% as.data.frame
    ex2_cell_lines <- ex2_cell_lines$CCLE_name
    
    #make a data frame for the affy analysis
    df <- makeRespVsGeneticDataFrame(full_con, gene='EGFR',
                               cell_lines=ex2_cell_lines,
                               drug='Erlotinib',
                               data_types = 'affy',
                               drug_df = NULL) 
    
    #scatter plot of EGFR expression vs Erlotinib response
    plotRespVsGeneticHist(df, 'affy', FALSE)
    
    #histogram of Erlotinib response coloured by EGFR expression
    plotRespVsGeneticPoint(df, 'affy', FALSE)
    

## ----fig.width=6, fig.height=4-------------------------------------------

    #make a data frame for the affy analysis
    df <- makeRespVsGeneticDataFrame(full_con, gene='BRAF',
                               cell_lines=ex2_cell_lines,
                               drug='PLX4720',
                               data_types = 'hybcap',
                               drug_df = NULL) 
    
    #scatter plot of EGFR expression vs Erlotinib response
    plotRespVsGeneticHist(df, 'hybcap', FALSE)
    
    #histogram of Erlotinib response coloured by EGFR expression
    plotRespVsGeneticPoint(df, 'hybcap', FALSE)
    

## ----fig.width=6, fig.height=6-------------------------------------------
    
    #get lung cell lines
    ex4_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% filter(Site_primary == 'lung') %>%
       collect %>% as.data.frame
    ex4_cell_lines <- ex4_cell_lines$CCLE_name

    #make the data frame
    gvg.df <- makeGeneticVsGeneticDataFrame(full_con, 
                                            cell_lines=ex4_cell_lines,
                                            gene1='SMARCA4',
                                            data_type1='hybcap',
                                            gene2='SMARCA4',
                                            data_type2='affy') 
    
    #view the data frame
    head(gvg.df)
    
    #do the plot
    plotGeneticVsGeneticPoint(gvg.df)
    
    #all in one go with axes swapped
    makeGeneticVsGeneticDataFrame(full_con, cell_lines=ex4_cell_lines, gene1='SMARCA4', data_type1='affy',
                                            gene2='SMARCA4', data_type2='hybcap') %>% plotGeneticVsGeneticPoint()
    
    #two continuous
    makeGeneticVsGeneticDataFrame(full_con, cell_lines=ex4_cell_lines, gene1='SMARCA4', data_type1='affy',
                                            gene2='SMARCA4', data_type2='cn') %>% plotGeneticVsGeneticPoint()
    
    #two discrete
    makeGeneticVsGeneticDataFrame(full_con, cell_lines=ex4_cell_lines, gene1='SMARCA4', data_type1='hybcap',
                                            gene2='KRAS', data_type2='hybcap') %>% plotGeneticVsGeneticPoint()
    
    #also plot by cell line with one feature a y axis and another as fill colour
    #continous + discrete
    makeGeneticVsGeneticDataFrame(full_con, cell_lines=ex4_cell_lines, gene1='SMARCA4', data_type1='affy',
                                            gene2='SMARCA4', data_type2='hybcap') %>% plotGeneticVsGeneticHist()
    
    #continous + continous
    makeGeneticVsGeneticDataFrame(full_con, cell_lines=ex4_cell_lines[1:25], gene1='SMARCA4', data_type1='affy',
                                            gene2='SMARCA4', data_type2='cn') %>% plotGeneticVsGeneticHist(label_option = TRUE)
    


## ----eval=FALSE----------------------------------------------------------
#  
#      dietlein_data_fn <- system.file("extdata", "Dietlein2014_supp_table_1.txt", package = "CancerCellLines")
#      dietlein_data <- read.table(dietlein_data_fn, header=T, sep='\t', stringsAsFactors=F)
#      head(dietlein_data)
#      dietlein_data <- dietlein_data %>%
#        filter(nchar(CCLE_name) > 1) %>%
#        transmute(unified_id=CCLE_name, compound_id='KU60648', endpoint='pGI50', original=GI50, value=9-log10(GI50))
#      head(dietlein_data)
#  
#      full_con <- setupSQLite('~/BigData/CellLineData/CancerCellLines.db')
#      shinyRespVsGeneticApp(con=full_con, drug_df=dietlein_data)
#  

## ----eval=FALSE----------------------------------------------------------
#      shinyRespVsGeneticApp(con=full_con)

## ----eval=FALSE----------------------------------------------------------
#      shinyGeneticVsGeneticApp(con=full_con)

## ----fig.width=8, fig.height=8-------------------------------------------
    #get all cell lines
    ex5_cell_lines <- dplyr_con %>% tbl('ccle_sampleinfo') %>% 
       collect %>% as.data.frame
    ex5_cell_lines <- ex5_cell_lines$CCLE_name
    
    #make a data frame
    df <- makeRespVsRespDataFrame(full_con, 
                               cell_lines=ex5_cell_lines,
                               drugs=c('Erlotinib', 'AZD6244'),
                               tissue_info = 'ccle')
    head(df)
    
    #makes a wide data frame
    wide.df <- df %>% makeWideFromRespVsRespDataFrame()
    head(wide.df)
     
    #now do some plots
    plotRespVsRespWaterfall(filter(df, grepl('Erlotinib', ID)))
    plotRespVsRespDensity(df)
    plotRespVsRespPairs(df)

## ----eval=FALSE----------------------------------------------------------
#      shinyRespVsRespApp(con=full_con)

## ------------------------------------------------------------------------
   sessionInfo() 

