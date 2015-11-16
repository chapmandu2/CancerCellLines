#' Compare two genetic features (histogram)
#'
#' Generate a histogram type plot that compares two genetic features in a set of cell lines
#'
#' @param df A data frame from the \code{makeGeneticVsGeneticDataFrame} function
#' @return A ggplot2 object
#' @export
plotGeneticVsGeneticHist <- function(df, facet_option=FALSE) {

  plot_data <- df %>% filter(!is.na(feature_value1)) %>% arrange(desc(feature_value1))
  dt1 <- ifelse(unique(plot_data$feature_type1) %in% c('hybcap', 'cosmicclp'), 'discrete', 'cont')
  dt2 <- ifelse(unique(plot_data$feature_type2) %in% c('hybcap', 'cosmicclp'), 'discrete', 'cont')

  if(dt1 =='discrete' ) {

    stop('Feature 1 must be a continous variable for this plot type')

  } else if (dt1 =='cont' & dt2 == 'cont' ) {

    p <- ggplot(plot_data, aes(x=CCLE_name, y=feature_value1, fill=scale(feature_value2)) ) +
      geom_point(shape=21, size=rel(5)) +
      scale_x_discrete(limits=plot_data$CCLE_name) +
      scale_fill_gradient2(low='blue', mid='white', high='red', name=unique(plot_data$feature_name2)) +
      ylab(unique(plot_data$feature_name1)) +
      theme_bw() + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))


  } else if (dt1 == 'cont' & dt2 == 'discrete') {
    plot_data <- plot_data %>% mutate(feature_value2=as.factor(feature_value2))
    p <- ggplot(plot_data, aes(x=CCLE_name, y=feature_value1, fill=feature_value2) ) +
      geom_point(shape=21, size=rel(5)) +
      scale_x_discrete(limits=plot_data$CCLE_name) +
      scale_fill_discrete(name=unique(plot_data$feature_name2)) +
      ylab(unique(plot_data$feature_name1)) +
      theme_bw() + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))

  } else {
    stop("I'm confused")
  }

  if(nrow(plot_data) > 30 | facet_option) {
    p <- p + theme(axis.text.x=element_text(size=0))  #remove cell line id's if more than 30 cell lines
  }

  if (facet_option) {
    p <- p + facet_wrap(~tissue)
  }

  return(p)

}
