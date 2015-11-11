#' Compare response against a genetic feature (histogram)
#'
#' Generate a plot that compares response in cell lines with different genetic features
#'
#' @param df A data frame from the \code{makeRespVsGeneticDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsGeneticHist <- function(df, data_type=NULL, facet_option=FALSE) {

  plot_data <- df %>% filter(!is.na(resp_value)) %>% arrange(desc(resp_value))

  if(data_type %in% c('hybcap', 'cosmicclp')) {
    p <- ggplot(plot_data, aes(x=CCLE_name, y=resp_value, fill=as.factor(feature_value)) ) +
      geom_bar(stat='identity') +
      scale_x_discrete(limits=plot_data$CCLE_name) +
      ylab(unique(plot_data$resp_id)) +
      theme_bw() + theme(axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))

    if(nrow(plot_data) > 30 | facet_option) {
      p <- p + theme(axis.text.x=element_text(size=0))  #remove cell line id's if more than 30 cell lines
    }

    if (facet_option) {
      p <- p + facet_wrap(~tissue)
    }

    return(p)

  } else {
    p <- ggplot(plot_data, aes(x=CCLE_name, y=resp_value, fill=scale(feature_value)) ) +
      geom_bar(stat='identity') +
      scale_x_discrete(limits=plot_data$CCLE_name) +
      scale_fill_gradient2(low='blue', mid='white', high='red') +
      ylab(unique(plot_data$resp_id)) +
      theme_bw() + theme(legend.position='none', axis.text.x = element_text(size=rel(1), angle=330, hjust=0, vjust=1))

    if(nrow(plot_data) > 30 | facet_option) {
      p <- p + theme(axis.text.x=element_text(size=0))  #remove cell line id's if more than 30 cell lines
    }

    if (facet_option) {
      p <- p + facet_wrap(~tissue)
    }

    return(p)
  }


}
