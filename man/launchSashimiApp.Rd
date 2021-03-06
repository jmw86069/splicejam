% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicejam-shiny.R
\name{launchSashimiApp}
\alias{launchSashimiApp}
\title{Launch Sashimi R-shiny application}
\usage{
launchSashimiApp(..., options = list(width = 1200))
}
\arguments{
\item{...}{additional arguments are passed to \code{shiny::shinyApp()}.}
}
\description{
Launch Sashimi R-shiny application
}
\details{
This function launches the Sashimi visualization
R-shiny app.

The R objects required to prepare sashimi plots are
defined by the function \code{sashimiAppConstants()}, which
documents each R object required and how it is used.

The most straightforward way to run a new Sashimi R-shiny
app is to define \code{filesDF} and \code{gtf} in the global environment.
The \code{gtf} is a path or URL to a GTF file (which can be gzipped).
This GTF file will be used to derive all related annotation
data.
\itemize{
\item \code{txdb} -- TranscriptDb from which other objects are derived
\item \code{tx2geneDF} -- \code{data.frame} with transcript-to-gene relationship
\item \code{detectedTx} -- if not already defined, all transcripts are
used. \emph{Much better to use only a subset of detected transcripts.}
\item \code{detectedGenes} -- inferred from \code{detectedTx}, using \code{tx2geneDF}.
\item \code{flatExonsByGene}, \code{flatExonsByTx} -- these objects will combine
CDS exons and non-CDS exons to represent CDS and UTR regions.
}

Note that if \code{detectedTx} is not defined, it will use all transcripts
at this stage, which can be substantially slower than using only
the subset of "observed/detected" transcripts.

The first time running \code{launchSashimiApp()}
will populate several R objects in the global environment,
and these objects will be re-used during subsequent calls to
this function. To make changes in the content, these objects
can be edited or deleted so the object is created again.
For example, if \code{detectedTx} is edited, the object
\code{detectedGenes} should be removed so \code{detectedGenes}
will be created again during the next call to \code{launchSashimiApp()}.

The \code{filesDF} object should be a \code{data.frame} with colnames
\code{"sample_id"}, \code{"type"} (with values either \code{"bw"} or \code{"junction"}),
and \code{"url"} (a URL or file path to each file.) If coverage
or junctions are available in separate files, use the same
\code{sample_id} value for each file. Files with the same \code{sample_id}
value are combined using the sum, after multiplying each file
by a value in the optional \code{"scale_factor"} column.

This function calls \code{sashimiAppConstants()} which does the
heavy work of defining or deriving all necessary data objects,
then assigns the result to the relevant environment. The default
environment is \code{globalenv()} (also known as \code{.GlobalEnv}).
}
\examples{
# Note: disabled for web page examples
# launchSashimiApp();

}
\seealso{
Other splicejam R-shiny functions: 
\code{\link{sashimiAppConstants}()},
\code{\link{sashimiAppServer}()},
\code{\link{sashimiAppUI}()}
}
\concept{splicejam R-shiny functions}
