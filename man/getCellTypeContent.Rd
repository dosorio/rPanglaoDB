% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCellTypeContent.R
\name{getCellTypeContent}
\alias{getCellTypeContent}
\title{Get the cell-type content for each sample from the panglaoDB database.}
\usage{
getCellTypeContent(
  sra = "All",
  srs = "All",
  tissue = "All",
  protocol = "All",
  specie = "All",
  verbose = TRUE
)
}
\arguments{
\item{sra}{A}

\item{srs}{A}

\item{tissue}{A}

\item{protocol}{A}

\item{specie}{A}

\item{verbose}{A}
}
\description{
Get the cell-type content for each sample from the panglaoDB database.
}
\examples{
# From PanglaoDB
# https://panglaodb.se/list_clusters_and_cell_types.html?sra=SRA689041&srs=SRS3166675

\dontrun{
SRS3166675 <- getCellTypeContent(srs = 'SRS3166675')
SRS3166675
}
}
