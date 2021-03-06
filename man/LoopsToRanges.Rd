% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoopsToRanges.R
\name{LoopsToRanges}
\alias{LoopsToRanges}
\title{Create a list of ranges objects from looping data}
\usage{
LoopsToRanges(..., loop_names = NULL, custom_cols = 0,
  custom_mcols = NULL)
}
\arguments{
\item{...}{Any number of tab delimited loop data files in BEDPE format}

\item{loop_names}{A character vector of names for the loop datasets (optional)}

\item{custom_cols}{An integer indicating the number of extra columns in the BEDPE file (default = 0)}

\item{custom_mcols}{An integer or vector of integers indicating which columns are used for metadata (optional)}
}
\value{
A \emph{LoopRanges} class object: list of GRanges looping data objects
}
\description{
Uses tab delimited looping data in the form of BEDPE files to create custom GRanges loop objects
}
\examples{
# Load loops into LoopRanges object 
ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", 
package = "LoopRig", mustWork = TRUE) 
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", 
package = "LoopRig", mustWork = TRUE)
LoopsToRanges(ovary_loops, spleen_loops, custom_cols = 0)
}
