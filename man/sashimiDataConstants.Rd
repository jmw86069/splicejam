% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicejam-shiny-data.R
\name{sashimiDataConstants}
\alias{sashimiDataConstants}
\title{Prepare sashimi plot required data}
\usage{
sashimiDataConstants(
  gtf = NULL,
  txdb = NULL,
  tx2geneDF = NULL,
  exonsByTx = NULL,
  cdsByTx = NULL,
  detectedTx = NULL,
  detectedGenes = NULL,
  flatExonsByGene = NULL,
  flatExonsByTx = NULL,
  envir = NULL,
  empty_uses_farrisdata = TRUE,
  use_memoise = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{gtf, txdb, tx2geneDF, exonsByTx, cdsByTx}{objects used to define
the overall set of genes, transcripts, and associated exons and
CDS exons. See this function description for more detail.
Notes:
\itemize{
\item \code{gtf} can be a local file, in which case it will be loaded
from its current path without copying to the current directory,
as long as \code{file.exists(gtf)} is \code{TRUE}, which also means
it would not have the prefix \code{"file://"}.
However, when it has the prefix \code{"file://"} it will use
\code{curl::curl_download()} which will copy it to the current
directory.
In either case, when other files are derived, such as
\code{tx2gene} or \code{txdb}, those files are stored in the current
directory.
\item \code{tx2geneDF} is a \code{data.frame} with at minimum two columns:
\code{gene_name} and \code{transcript_id}. It is possible to customize
these column names, however it is easiest to use these defaults.
\item \code{exonsByTx} and \code{cdsByTx} are \code{GRangesList} objects, and
they are used to derive \code{flatExonsByTx} and \code{flatExonsByGene}
when those objects are not already provided.
Also, \code{exonsByTx} and \code{cdsByTx} are derived using \code{gtf} or \code{txdb}
when necessary to derive \code{flatExonsByTx} and \code{flatExonsByGene}.
}}

\item{detectedTx, detectedGenes, flatExonsByGene, flatExonsByTx}{objects used to derive a specific subset of gene-exon models
using only detected transcripts or genes. See this function
description for more detail.
\itemize{
\item \code{detectedTx} and \code{detectedGenes} are \code{character} vectors.
\item \code{flatExonsByTx} and \code{flatExonsByGene} are \code{GRangesList},
where each \code{GRanges} element contains disjoint (non-overlapping)
ranges. When not provided, they are derived from \code{exonsByTx}
and \code{cdsByTx}, which also requires \code{tx2geneDF} and either
\code{gtf} or \code{txdb}.
}}

\item{envir}{\code{environment} where data will be prepared, or when
\code{envir=NULL} a new environment will be created and returned.}

\item{empty_uses_farrisdata}{\code{logical} indicating whether to
use data from the Github R package \code{"jmw86069/farrisdata"}
if no data is supplied to this function. This behavior is
intended to make it easy to use farrisdata to recreate
the Sashimi plots in that publication.}

\item{use_memoise}{\code{logical} indicating whether to use \code{memoise}
to cache intermediate data files for exons, flattened exons,
transcript-gene data, and so on. This mechanism reduces
time to render sashimi plots that re-use the same gene.
All memoise cache folders are named with \code{"_memoise"}.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}

\item{default_gene}{\code{character} string indicating the default
gene to use for the initial R-shiny figure.}
}
\value{
\code{environment} that contains the required data objects
for splicejam sashimi plots. Note that the environment itself
is updated during processing, so the environment does not
need to be returned for the data contained inside it to
be updated by this function.
}
\description{
Prepare sashimi plot required data, deriving data objects as needed
}
\details{
This function performs a subset of steps performed by
\code{sashimiAppConstants()}, focusing only on data required
for gene-exon structure. The \code{sashimiAppConstants()} defines
\code{color_sub} and validates \code{filesDF}, then calls this function
\code{sashimiDataConstants()} to prepare and validate the gene-exon
data.

Data derived by this function \code{sashimiDataConstants()}:
\itemize{
\item \strong{txdb}: \code{TranscriptDb} object used to derive \code{exonsByTx}
and \code{cdsByTx} if either object does not already exist. If \code{txdb}
is not supplied, it is derived from \code{gtf} using
\code{GenomicFeatures::makeTxDbFromGFF()}.
\item \strong{tx2geneDF}: \code{data.frame} with colnames: \code{"transcript_id"} and
\code{"gene_name"}.
\item \strong{gtf}: \code{character} path to a GTF/GFF/GFF3 file, suitable for
\code{GenomicFeatures::makeTxDbFromGFF()}. The \code{gtf} is only used
if \code{tx2geneDF} or \code{exonsByTx} are not supplied. Note that
when \code{gtf} points to a remote server, the file is copied to
the current working directory for more rapid use.
If the file already exists in the local directory, it is re-used.
\item \strong{exonsByTx}: \code{GRangesList} object, named by \code{"transcript_id"},
containing all exons for each transcript. It is derived from \code{txdb}
if not supplied; and names should match \code{tx2geneDF$transcript_id}.
\item \strong{cdsByTx}: \code{GRangesList} object, named by \code{"transcript_id"},
containing only CDS (protein-coding) exons for each transcript.
It is derived from \code{txdb} if not supplied;
and names should match \code{tx2geneDF$transcript_id}.
\item \strong{detectedTx}: \code{character} vector of \code{tx2geneDF$transcript_id} values,
representing a subset of transcripts detected above background.
See \code{definedDetectedTx()} for one strategy to define detected transcripts.
If \code{detectedTx} does not exist, it is defined by all transcripts
present in \code{tx2geneDF$transcript_id}. Note this step can be the
rate-limiting step in the preparation of \code{flatExonsByTx}.
\item \strong{detectedGenes}: \code{character} vector of values that match
\code{tx2geneDF$gene_name}. If it is not supplied, it is inferred
from \code{detectedTx} and \code{tx2geneDF$transcript_id}.
\item \strong{flatExonsByGene}: \code{GRangesList} object containing non-overlapping
exons for each gene, whose names match \code{tx2geneDF$gene_name}. If not
supplied, it is derived using \code{flattenExonsBy()} and objects
\code{exonsByTx}, \code{cdsByTx}, \code{detectedTx}, and \code{tx2geneDF}. This step is
the key step for using a subset of detected transcripts, in order
to produce a clean gene-exon model.
\item \strong{flatExonsByTx}: \code{GRangesList} object containing non-overlapping
exons for each transcript. If not
supplied, it is derived using \code{flattenExonsBy()} and objects
\code{exonsByTx}, \code{cdsByTx}, \code{detectedTx}, and \code{tx2geneDF}. This step is
the key step for using a subset of detected transcripts, in order
to produce a clean transcript-exon model.
}

When \code{use_memoise=TRUE} several R objects are cached using
\code{memoise::memoise()}, to help re-use of prepared R objects,
and to help speed the re-use of data within the R-shiny app:
}
\seealso{
Other splicejam R-shiny functions: 
\code{\link{launchSashimiApp}()},
\code{\link{sashimiAppConstants}()},
\code{\link{sashimiAppUI}()}
}
\concept{splicejam R-shiny functions}
