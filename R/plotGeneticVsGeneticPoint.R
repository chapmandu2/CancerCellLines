#' Compare two genetic features (scatter/boxplot)
#'
#' Generate a plot that compares two genetic features in a set of cell lines
#'
#' @param df A data frame from the \code{makeGeneticVsGeneticDataFrame} function
#' @return A ggplot2 object
#' @export
plotGeneticVsGeneticPoint <- function(df, facet_option=FALSE) {

  plot_data <- df
  dt1 <- ifelse(unique(plot_data$feature_type1) %in% c('hybcap', 'cosmicclp'), 'discrete', 'cont')
  dt2 <- ifelse(unique(plot_data$feature_type2) %in% c('hybcap', 'cosmicclp'), 'discrete', 'cont')

  if(dt1 =='discrete' & dt2 == 'cont' ) {
    plot_data <- plot_data %>% mutate(feature_value1=as.factor(feature_value1))
    p <- ggplot(plot_data, aes(x=feature_value1, y=feature_value2) ) +
      geom_boxplot(outlier.size=0, aes(colour=feature_value1)) +
      geom_point(aes(fill=feature_value1), shape=21, size=rel(3), position = position_jitter(width=0.3)) +
      xlab(unique(plot_data$feature_name1)) + ylab(unique(plot_data$feature_name2)) +
      theme_bw()

  } else if (dt1 =='discrete' & dt2 == 'discrete' ) {

    plot_data <- plot_data %>% mutate(feature_value1=as.factor(feature_value1),
                                      feature_value2=as.factor(feature_value2))
    p <- ggplot(plot_data, aes(x=feature_value1, y=feature_value2) ) +
      geom_point(aes(fill=feature_value1), shape=21, size=rel(3), position = position_jitter(width=0.3, height=0.3)) +
      xlab(unique(plot_data$feature_name1)) + ylab(unique(plot_data$feature_name2)) +
      theme_bw()

  } else if (dt1 =='cont' & dt2 == 'cont' ) {

    p <- ggplot(plot_data, aes(x=feature_value1, y=feature_value2) ) +
      geom_point(aes(fill=scale(feature_value1)), shape=21, size=rel(3)) +
      stat_smooth(method = 'lm') +
      scale_fill_gradient2(low='blue', mid='white', high='red') +
      xlab(unique(plot_data$feature_name1)) + ylab(unique(plot_data$feature_name2)) +
      theme_bw() + theme(legend.position='none')


  } else if (dt1 == 'cont' & dt2 == 'discrete') {
    plot_data <- plot_data %>% mutate(feature_value2=as.factor(feature_value2))
    p <- ggplot(plot_data, aes(x=feature_value2, y=feature_value1) ) +
      geom_boxplot(outlier.size=0, aes(colour=feature_value2)) +
      geom_point(aes(fill=feature_value2), shape=21, size=rel(3), position = position_jitter(width=0.3)) +
      xlab(unique(plot_data$feature_name2)) + ylab(unique(plot_data$feature_name1)) + coord_flip() +
      theme_bw()

  } else {
    stop("I'm confused")
  }

  if (facet_option) {
    p <- p + facet_wrap(~tissue)
  }

  return(p)

}

