#script used to make the toy data

#location of files
infopath <- '~/BigData/CellLineData/RawData/CCLE_sample_info_file_2012-10-18.txt'
affypath <- '~/BigData/CellLineData/RawData/CCLE_Expression_Entrez_2012-09-29.gct'
cnpath <- '~/BigData/CellLineData/RawData/CCLE_copynumber_byGene_2012-09-29.txt'
hybcappath <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf'
cosmicclppath <- '~/BigData/CellLineData/RawData/CosmicCLP_CompleteExport_v74.tsv'
drugpath <- '~/BigData/CellLineData/RawData/CCLE_NP24.2009_Drug_data_2012.02.20.csv'
idspath <- system.file("extdata", "CellLineIDNormalisationNov15.txt", package = "CancerCellLines")

#location of toy data directory
toydir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

#define cell line pattern
toypattern <- '_LUNG|_BREAST'

#run the functions to make the toy data
CancerCellLines:::maketoyCCLE_info(infopath, toypattern, toydir)
CancerCellLines:::maketoyCCLE_affy(affypath, toypattern, toydir)
CancerCellLines:::maketoyCCLE_cn(cnpath, toypattern, toydir)
CancerCellLines:::maketoyCCLE_hybcap(hybcappath, toypattern, toydir)
CancerCellLines:::maketoyCosmicCLP_exome(cosmicclppath, toypattern, toydir)
CancerCellLines:::maketoyCCLE_drugresponse(drugpath, toypattern, toydir)
CancerCellLines:::maketoyCellLineIDs(idspath, toypattern, toydir)

#make the toy database
dest <- makeToyDB()
file.copy(dest@dbname, paste0(toydir, 'toy.db'))


