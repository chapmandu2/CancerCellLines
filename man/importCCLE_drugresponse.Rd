% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importCCLE_drugresponse.R
\name{importCCLE_drugresponse}
\alias{importCCLE_drugresponse}
\title{Import CCLE drug response data}
\usage{
importCCLE_drugresponse(fn, con)
}
\arguments{
\item{fn}{The path of the data file}

\item{con}{A \code{SQLiteConnection} object to the database}
}
\value{
TRUE or FALSE depending on whether the data has been written successfully
}
\description{
This function imports the information in CCLE_NP24.2009_Drug_data_2012.02.20.csv into the ccle_drugresponse table in the database.
}

