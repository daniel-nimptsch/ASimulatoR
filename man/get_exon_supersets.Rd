% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_exon_supersets.R
\name{get_exon_supersets}
\alias{get_exon_supersets}
\title{Create or load exon_supersets from gtf/gff file}
\usage{
get_exon_supersets(gtf_path, ncores = 1L, save = TRUE)
}
\arguments{
\item{gtf_path}{path to the gtf/gff file from which exon supersets are created}

\item{ncores}{number of cores used to generate exon supersets in parallel.
Default \code{1}: no parallel computing.
\strong{This will spawn one process per core! Be aware that a lot of memory might be required for many processes.}}

\item{save}{should the exon_supersets be saved to .rda file?
Default \code{TRUE}}
}
\value{
exon supersets generated from gtf/gff file
}
\description{
Create or load exon_supersets from gtf/gff file
}
