% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jambio-gtf.R
\name{readGtf}
\alias{readGtf}
\title{read gtf/gff3 file}
\usage{
readGtf(GTF, nrows = -1, zcat_command = "zcat", verbose = FALSE, ...)
}
\arguments{
\item{GTF}{\code{character} file name sent to \code{data.table::fread()}. When the
file ends with ".gz", the \code{R.utils} package is recommended, otherwise
the fallback option is to make a system call to \code{gzcat}
to gunzip the file during the import step. Note this process fails
when \code{gzcat} is not available in the path of the user environment.
In general, the \code{R.utils} package is the best solution.}

\item{nrows}{\code{integer} number of rows to read from the GTF file, by default
-1 means all rows are imported. This parameter is useful to check the
results of a large GTF file using only a subset portion of the file.}

\item{zcat_command}{\code{character} name or path to zcat or gzcat executable,
only used when input \code{GTF} is a file with \code{".gz"} extension, and when
R package \code{R.utils} is not available.}

\item{verbose}{\code{logical} whether to print verbose output during processing.}

\item{...}{additional arguments are ignored.}
}
\value{
\code{data.frame} representing tab-delimited data stored in the
gtf or gff3 file.
}
\description{
read gtf/gff3 file
}
