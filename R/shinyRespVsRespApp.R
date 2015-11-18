
#' Response vs Response Shiny visualisation
#'
#' Interactive plots comparing and visualising different response variables in Shiny
#'
#' @param con A \code{SQLiteConnection} object to the database
#' @return Launches an interactive Shiny application
#' @export
shinyRespVsRespApp <- function(con, drug_df=NULL) {

  shinyApp(
    ui = shinyRespVsRespUI(),
    server = function(input, output) {
      shinyRespVsRespServer(input, output, con, drug_df)
    }
  )
}
