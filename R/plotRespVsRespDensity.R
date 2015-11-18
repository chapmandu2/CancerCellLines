#' Density plots for selected response variables
#'
#' Generate a set of density plots for the selected response variables
#'
#' @param df A data frame from the \code{makeRespVsRespDataFrame} function
#' @return A ggplot2 object
#' @export
plotRespVsRespDensity <- function(df) {

  p <- ggplot(df, aes(x=value, fill=ID, colour=ID)) + geom_density(alpha=0.7) + facet_grid(ID~.) + theme_bw() + theme(legend.position='none')
  return(p)


}
