#' Correlation plots for selected response variables
#'
#' Generate a set of scatter plots for the selected response variables
#'
#' @param df A data frame from the \code{makeRespVsRespDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsRespPairs <- function(df) {

  #check that there are least 2 ids!
  if(length(unique(df$assayed_id)) < 2) {
    stop('Need at least one response variable in the supplied data frame')
  }

  plot_data <- df %>%
    dplyr::inner_join(df, by=c('unified_id', 'tissue')) %>%
    dplyr::select(-starts_with('data_type'), -starts_with('original'), -starts_with('subtype')) %>%
    dplyr::filter(assayed_id.x != assayed_id.y)

  p <- ggplot(plot_data, aes(value.x, value.y)) +
    geom_point(shape=21, size=rel(3), aes(fill=tissue)) +
    stat_smooth(method='lm') +
    facet_grid(assayed_id.x ~ assayed_id.y) +
    theme_bw()
  p

}
