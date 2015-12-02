#functions to make the toy data

maketoyIDs <- function() {
  genes <- c('PTEN', 'TP53', 'BRAF', 'NRAS',  'CRAF', 'EGFR', 'SMARCA4', 'KRAS')
  drugs <- c('AZD6244','PLX4720','PD-0325901', 'Erlotinib', 'Nilotinib', 'Lapatinib')
  return(list(genes=genes, drugs=drugs))
}

maketoyCCLE_info <- function (fn, pattern='_BREAST', outdir=getwd()) {

  data <- readLines(fn)
  header <- data[1]
  patternMatch <- data [ grepl(pattern, data)  ]

  out <- c(header, patternMatch)
  outfn <- gsub('.txt', '_toy.txt', basename(fn))
  write( out, paste0(outdir,outfn))

}

maketoyCCLE_affy <- function (fn, pattern='_BREAST', outdir=getwd()) {

  data <- read.table(fn,header=T, sep='\t', skip=2, stringsAsFactors = FALSE)
  reqcols <- grepl(pattern, colnames(data))
  reqcols[1:2] <- TRUE
  reqrows <- data$Description %in% CancerCellLines:::maketoyIDs()$genes
  data <- data [ reqrows , reqcols ]

  header <- readLines(fn, n=2)

  outfn <- gsub('.gct', '_toy.gct', basename(fn))
  outpath <- paste0(outdir,outfn)
  write( header, outpath)
  write.table(data, outpath, sep='\t', append=TRUE, row.names=FALSE, quote = FALSE)

}

maketoyCCLE_cn <- function (fn, pattern='_BREAST', outdir=getwd()) {

  data <- read.table(fn,header=T, sep='\t', stringsAsFactors = FALSE)
  reqcols <- grepl(pattern, colnames(data))
  reqcols[1:4] <- TRUE
  reqrows <- data$geneName %in% CancerCellLines:::maketoyIDs()$genes
  data <- data [ reqrows , reqcols ]

  outfn <- gsub('.txt', '_toy.txt', basename(fn))
  outpath <- paste0(outdir,outfn)

  write.table(data, outpath, sep='\t', append=FALSE, row.names=FALSE, quote = FALSE)

}

maketoyCCLE_hybcap <- function (fn, pattern='_BREAST', outdir=getwd()) {

  require(readr)

  data <- read_tsv(fn,
                   col_names=c("Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
                               "Start_position", "End_position", "Strand", "Variant_Classification",
                               "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                               "dbSNP_RS", "dbSNP_Val_Status", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                               "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Tumor_Validation_Allele1",
                               "Tumor_Validation_Allele2", "Match_Norm_Validation_Allele1",
                               "Match_Norm_Validation_Allele2", "Verification_Status", "Validation_Status",
                               "Mutation_Status", "Sequencing_Phase", "Sequence_Source", "Validation_Method",
                               "Score", "BAM_file", "Sequencer", "Genome_Change", "Annotation_Transcript",
                               "Transcript_Strand", "cDNA_Change", "Codon_Change", "Protein_Change",
                               "Other_Transcripts", "Refseq_mRNA_Id", "Refseq_prot_Id", "SwissProt_acc_Id",
                               "SwissProt_entry_Id", "Description", "UniProt_AApos", "UniProt_Region",
                               "UniProt_Site", "Alternative_allele_reads_count", "Reference_allele_reads_count",
                               "X46vertebrates_AA_alignment_column", "Method"),
                   col_types=paste(rep('c', 51), collapse=''),
                   skip=1)

  out <- data [ grepl(pattern, data$Tumor_Sample_Barcode ) & data$Hugo_Symbol %in% CancerCellLines:::maketoyIDs()$genes , ]

  outfn <- gsub('.maf', '_toy.maf', basename(fn))
  write.table( out, paste0(outdir,outfn), row.names = FALSE, sep='\t')

}

maketoyCosmicCLP_exome <- function (fn, pattern='_BREAST', outdir=getwd()) {

  require(readr)
  data <- read_tsv(fn)

  out <- data [ grepl(pattern, paste0('_', toupper(data$`Primary site`))) & data$`Gene name` %in% CancerCellLines:::maketoyIDs()$genes, ]

  outfn <- gsub('.tsv', '_toy.tsv', basename(fn))
  write.table( out, paste0(outdir,outfn), row.names = FALSE, sep='\t' )

}

maketoyCCLE_drugresponse <- function (fn, pattern='_BREAST', outdir=getwd()) {

  data <- readLines(fn)
  header <- data[1]
  patternMatch <- data [ grepl(pattern, data) & grepl(paste0(CancerCellLines:::maketoyIDs()$drugs, collapse='|'), data)  ]

  out <- c(header, patternMatch)
  outfn <- gsub('.csv', '_toy.txt', basename(fn))
  write( out, paste0(outdir,outfn))

}

maketoyCellLineIDs <- function (fn, pattern='_BREAST', outdir=getwd()) {

  fn <- '~/BigData/CellLineData/RawData/CellLineIDNormalisationNov15.txt'
  outdir <- '~/Documents/2015 Projects/20150622 CancerCellLine Package/ToyData/'

  data <- read.table(fn, header=T, sep='\t')

  out <- data [ grepl(pattern, data$unified_id) , ]

  outfn <- gsub('.txt', '_toy.txt', basename(fn))
  write.table( out, paste0(outdir,outfn), row.names = FALSE, sep='\t' )

}
