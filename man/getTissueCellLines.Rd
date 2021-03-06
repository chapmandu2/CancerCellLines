% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTissueCellLines.R
\name{getTissueCellLines}
\alias{getTissueCellLines}
\title{Convert tissue to cell lines}
\usage{
getTissueCellLines(con, tissue_id, tissue_info = "ccle")
}
\arguments{
\item{con}{A \code{SQLiteConnection} object to the database}

\item{tissue_id}{A vector of tissue identifiers}

\item{tissue_info}{Source of tissue information: \code{ccle} uses the ccle sample information, \code{crude} ensures that all cell lines are returned but mixes data from different sources, \code{custom} can be one of any id_type present in the cell_line_ids table for extra control}
}
\value{
A \code{vector} of cell lines for the provided tissue(s)
}
\description{
This function creates a \code{vector} of cell lines from a tissue for a given type of tissue ifo
}

