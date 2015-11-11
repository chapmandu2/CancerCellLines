
#' Runs the Shiny visualisation
#'
#' Takes a user supplied data frame of drug response data and combines it with user defined gene feature data
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @param drug_df A data frame containing the drug data
#' @return Launches an interactive Shiny application
#' @export
shinyRespVsGeneticApp <- function(con, drug_df=NULL) {

  shinyApp(
    ui = shinyRespVsGeneticUI(),
    server = function(input, output) {
      shinyRespVsGeneticServer(input, output, con, drug_df)
    }
  )
}
