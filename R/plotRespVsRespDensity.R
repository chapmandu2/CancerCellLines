#' Density plots for selected response variables
#'
#' Generate a set of density plots for the selected response variables
#'
#' @param df A data frame from the \code{makeRespVsRespDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsRespDensity <- function(df) {

  p <- ggplot(df, aes(x=value, fill=assayed_id, colour=assayed_id)) + geom_density(alpha=0.7) + facet_grid(assayed_id~.) + theme_bw() + theme(legend.position='none')
  return(p)


}
