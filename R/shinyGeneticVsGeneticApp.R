
#' Genetic vs Genetic Shiny visualisation
#'
#' Interactive plots from user defined genes and data types in Shiny
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @return Launches an interactive Shiny application
#' @export
#' @import shiny
shinyGeneticVsGeneticApp <- function(con) {

  shiny::shinyApp(
    ui = shinyGeneticVsGeneticUI(),
    server = function(input, output) {
      shinyGeneticVsGeneticServer(input, output, con)
    }
  )
}
