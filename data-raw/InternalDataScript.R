#script used to create internal data

#example dataset
dietlein_data_fn <- system.file("extdata", "Dietlein2014_supp_table_1.txt", package = "CancerCellLines")
dietlein_data <- read.table(dietlein_data_fn, header=T, sep='\t', stringsAsFactors=F)
head(dietlein_data)
dietlein_data <- dietlein_data %>%
  filter(nchar(CCLE_name) > 1) %>%
  transmute(unified_id=CCLE_name, compound_id='KU60648', endpoint='pGI50', original=GI50, value=-log10(GI50))

devtools::use_data(dietlein_data)
