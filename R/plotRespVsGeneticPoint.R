#' Compare response against a genetic feature (scatter/boxplot)
#'
#' Generate a plot that compares response in cell lines with different genetic features
#'
#' @param df A data frame from the \code{makeRespVsGeneticDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsGeneticPoint <- function(df, data_type=NULL, facet_option=FALSE) {

  plot_data <- df %>% dplyr::arrange(dplyr::desc(resp_value))

  if(data_type %in% c('hybcap', 'cosmicclp')) {
    p <- ggplot(plot_data, aes(x=as.factor(feature_value), y=resp_value) ) +
      geom_boxplot(outlier.size=0, aes(colour=as.factor(feature_value))) +
      geom_point(aes(fill=as.factor(feature_value)), shape=21, size=rel(2), position = position_jitter(width=0.1)) +
      xlab(unique(plot_data$feature_name)) + ylab(unique(plot_data$resp_id)) +
      theme_bw()

    if (facet_option) {
      p <- p + facet_wrap(~tissue)
    }

    return(p)

  } else {
    p <- ggplot(plot_data, aes(x=feature_value, y=resp_value) ) +
      geom_point(aes(fill=scale(feature_value)), shape=21, size=rel(2)) +
      stat_smooth(method = 'lm') +
      scale_fill_gradient2(low='blue', mid='white', high='red') +
      xlab(unique(plot_data$feature_name)) + ylab(unique(plot_data$resp_id)) +
      theme_bw() + theme(legend.position='none')

    if (facet_option) {
      p <- p + facet_wrap(~tissue)
    }
    return(p)

  }

}

