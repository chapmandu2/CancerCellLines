#functions to import raw text files to a SQLite database

#' Connect to or create a RSQlite connection
#'
#' This function connects to or creates if necessary a SQLite database and returns a SQLiteConnection object.
#'
#' @param dbpath The path of the database
#' @return A \code{SQLiteConnection} object
#' @export
setupSQLite <- function ( dbpath=system.file('extdata/toy.db', package="CancerCellLines") ) {
  require(RSQLite)
  drv <- dbDriver("SQLite")
  con <- dbConnect(drv, dbname = dbpath)
  return(con)
}
