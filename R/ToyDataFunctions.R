#functions to make the toy data

maketoyCCLE_info <- function (fn, pattern='_BREAST', outdir=getwd()) {

  fn <- '~/BigData/CellLineData/RawData/CCLE_sample_info_file_2012-10-18.txt'
  outdir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- readLines(fn)
  header <- data[1]
  patternMatch <- data [ grepl(pattern, data)  ]

  out <- c(header, patternMatch)
  outfn <- gsub('.txt', '_toy.txt', basename(fn))
  write( out, paste0(outdir,outfn))

}

maketoyCCLE_affy <- function (fn, pattern='_BREAST', outdir=getwd()) {

  fn <- '~/BigData/CellLineData/RawData/CCLE_Expression_Entrez_2012-09-29.gct'
  outdir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- read.table(fn,header=T, sep='\t', skip=2, stringsAsFactors = FALSE)
  reqcols <- grepl(pattern, colnames(data))
  reqcols[1:2] <- TRUE
  data <- data [  , reqcols ]

  header <- readLines(fn, n=2)

  outfn <- gsub('.gct', '_toy.gct', basename(fn))
  outpath <- paste0(outdir,outfn)
  write( header, outpath)
  write.table(data, outpath, sep='\t', append=TRUE, row.names=FALSE, quote = FALSE)

}

maketoyCCLE_cn <- function (fn, pattern='_BREAST', outdir=getwd()) {

  fn <- '~/BigData/CellLineData/RawData/CCLE_copynumber_byGene_2012-09-29.txt'
  outdir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- read.table(fn,header=T, sep='\t', stringsAsFactors = FALSE)
  reqcols <- grepl(pattern, colnames(data))
  reqcols[1:4] <- TRUE
  data <- data [  , reqcols ]

  outfn <- gsub('.txt', '_toy.txt', basename(fn))
  outpath <- paste0(outdir,outfn)

  write.table(data, outpath, sep='\t', append=FALSE, row.names=FALSE, quote = FALSE)

}

maketoyCCLE_hybcap <- function (fn, pattern='_BREAST', outdir=getwd()) {

  require(readxl)
  require(xlsx)
  fn <- '~/BigData/CellLineData/RawData/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.xlsx'
  outdir <- '/Users/pchapman/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- read_excel(fn, col_types=c('text', 'numeric', 'text', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'text',
                                     'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text', 'text'))

  out <- data [ grepl(pattern, data$Tumor_Sample_Barcode ) , ]

  outfn <- gsub('.xlsx', '_toy.xlsx', basename(fn))
  write.xlsx( out, paste0(outdir,outfn), row.names = FALSE)

}

maketoyCCLE_drugresponse <- function (fn, pattern='_BREAST', outdir=getwd()) {

  fn <- '~/BigData/CellLineData/RawData/CCLE_NP24.2009_Drug_data_2012.02.20.csv'
  outdir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- readLines(fn)
  header <- data[1]
  patternMatch <- data [ grepl(pattern, data)  ]

  out <- c(header, patternMatch)
  outfn <- gsub('.csv', '_toy.txt', basename(fn))
  write( out, paste0(outdir,outfn))

}
