% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jambio-plot.R
\name{prepareSashimi}
\alias{prepareSashimi}
\title{Prepare Sashimi plot data}
\usage{
prepareSashimi(
  flatExonsByGene = NULL,
  filesDF = NULL,
  gene,
  sample_id = NULL,
  minJunctionScore = 10,
  gapWidth = 200,
  addGaps = TRUE,
  baseline = 0,
  compressGR = TRUE,
  compress_introns = TRUE,
  ref2c = NULL,
  gap_feature_type = "intron",
  default_feature_type = "exon",
  feature_type_colname = "feature_type",
  exon_label_type = c("none", "repel", "mark"),
  junc_label_type = c("repel", "mark", "none"),
  return_data = c("df", "ref2c"),
  include_strand = c("both", "+", "-"),
  junc_color = jamba::alpha2col("goldenrod3", 0.7),
  junc_fill = jamba::alpha2col("goldenrod1", 0.4),
  doStackJunctions = TRUE,
  coord_method = c("coord", "scale", "none"),
  scoreFactor = 1,
  scoreArcFactor = 0.2,
  scoreArcMinimum = 100,
  covGR = NULL,
  juncGR = NULL,
  use_memoise = FALSE,
  memoise_coverage_path = "coverage_memoise",
  memoise_junction_path = "junctions_memoise",
  do_shiny_progress = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{flatExonsByGene}{GRangesList named by gene, whose GRanges
elements are flattened, disjoint, non-overlapping genomic ranges
per gene.}

\item{filesDF}{data.frame with columns \code{url}, \code{sample_id}, \code{type},
where: \code{url} is any valid file path or URL compatible with
\code{base::read.table()}; \code{sample_id} is an identified representing
a biological sample, used to group common files together;
\code{type} is one of \code{"bw"} for bigWig coverage, \code{"junction"} for
BED12 format splice junctions.}

\item{gene}{character string of the gene to prepare, which must be
present in \code{names(flatExonsByGene)}.}

\item{gapWidth}{numeric value of the fixed width to use for
gaps (introns) between exon features. If \code{NULL} then
\code{getGRgaps()} will use the default based upon the median exon
width.}

\item{addGaps}{logical indicating whether to include gap regions
in the coverage plot, for example including introns or intergenic
regions. When \code{compressGR=TRUE} then gaps regions are
down-sampled using running maximum signal with roughly the same
x-axis resolution as uncompressed regions.}

\item{baseline}{numeric vector named by \code{names(flatExonsByGene)}
where baseline is used to adjust the y-axis baseline position
above or below zero.}

\item{compressGR}{logical indicating whether to compress GRanges
coordinates in the output data, where gaps/introns are set
to a fixed width. When \code{ref2c} is not supplied, and
\code{compressGR=TRUE}, then \code{ref2c} is created using
\code{make_ref2compressed()}.}

\item{compress_introns}{logical indicating whether to compress
the coverage polygon coordinates to approximately the same
number of pixels per inch as the exon polygons. This option
greatly reduces the size of the polygon, since introns are
already about 50 to 100 times wider than exons, and when
\code{compressGR} is \code{TRUE}, the introns are visibly compressed
to a fixed width on the x-axis. The data has many more
x-axis coordinates than the data visualization, this argument
is intended to reduce the intron coordinates accordingly.}

\item{ref2c}{list object output from \code{make_ref2compressed()} used
to compress axis coordinates, to compress polygon coverage
data in compressed regions, and to adjust splice junction arcs
using compressed coordinates.}

\item{gap_feature_type}{the default feature_type value to use for
gaps when \code{addGaps=TRUE}.}

\item{include_strand}{character value, one of \code{"both"}, \code{"+"},
\code{"-"} indicating the strandedness of coverage and junctions
to display. The default \code{"both"} shows coverage on both strands,
otherwise coverage is filtered either by filename (presence of
\code{"pos"}, \code{"+"}, or \code{"plus"} indicates positive strand), or
by detecting strandedness by positive/negative coverage scores.
Detecting by filename is intended to avoid retrieving coverage
in the R-shiny app, to help efficiency.}

\item{doStackJunctions}{logical indicating whether to stack
junction arcs at each end, this argument is passed to
\code{grl2df()} which calls \code{stackJunctions()}.}

\item{covGR}{GRanges object containing coverage data in columns
stored as NumericList class, where \code{colnames(GenomicRanges::values(covGR))}
are present in \code{filesDF$url} when \code{files$type \%in\% "coverage_gr"}.}

\item{juncGR}{GRanges object containing splice junctions, where
\code{"score"} is used for the abundance of splice junction reads,
and \code{"sample_id"} is used to define the biological \code{sample_id}.}

\item{do_shiny_progress}{logical indicating whether to send
progress updates to a running shiny app, using the
\code{shiny::withProgress()} and \code{shiny::setProgress()} methods.
This function only calls \code{shiny::setProgress()} and
assumes the \code{shiny::withProgress()} has already been
initialized.}

\item{verbose}{logical indicating whether to print verbose output.}

\item{...}{additional arguments are passed to \code{make_ref2compressed()},
\code{getGRcoverageFromBw()}, \code{exoncov2polygon()}.}
}
\value{
list containing \code{ggSashimi} a ggplot2 graphical object
containing a full Sashimi plot; \code{ggCov} the RNA-seq coverage
subset of the Sashimi plot; \code{ggJunc} the splice junction
subsset of the Sashimi plot; \code{ref2c} the output of
\code{make_ref2compressed()} used for ggplot2 coordinate
visualization; \code{covDF}, \code{juncDF} data.frame objects
with the raw data used to create ggplot2 objects;
\code{covGR}, \code{juncGR} the GRanges objects used to create
the data.frames; \code{gr} the GRanges object representing the
exons for the gene of interest; \code{juncLabelDF} the data.frame
containing exon label coordinates used to add labels to
the splice junction arcs.
}
\description{
Prepare Sashimi plot data
}
\details{
This function is the workhorse function used to produce
Sashimi plots, and is intended to be a convenient wrapper
function for several other individual functions.

At a minimum, a Sashimi plot requires three things:
\enumerate{
\item Exons, usually from a gene of interest.
\item RNA-seq coverage data.
\item Splice junction data.
}

There is some required pre-processing before running
\code{prepareSashimi()}:
\itemize{
\item Prepare flattened exons by gene using \code{flattenExonsByGene()}
and corresponding data, including \code{exonsByGene}, \code{cdsByGene},
and \code{tx2geneDF}. Verify the gene exon model data using
\code{gene2gg()}.
\item Find file paths, or web URLs, for a set of bigWig coverage
files, representing RNA-seq coverage for each strand, for
the samples of interest. Test the coverage data using
\code{getGRcoverageFromBw()} for a small set of GRanges data.
\item Find file paths, or web URLs, for a set of BED6 or BED12
format files, note that it cannot currently use bigBed format
due to limitations in the \code{rtracklayer} package.
Test the splice junction data using \code{rtracklayer::import()}
for a small range of GRanges features, then send the data
to \code{spliceGR2junctionDF()} to prepare a data.frame summary.
}

The basic input for coverage and junction data is a data.frame,
which defines each file path or url, the type of data
\code{"bw"} or \code{"junction"}, and the biological sample \code{"sample_id"}.
Any file path compatible with \code{rtracklayer::import()} will
work, including web URLs and local files. When using a web URL
you may need to use \code{"https://"} format to force the use
of secure web requests, but this requirement varies by country.
}
\examples{
# The active example below uses sample data
suppressPackageStartupMessages(library(GenomicRanges));

data(test_exon_gr);
data(test_junc_gr);
data(test_cov_gr);
filesDF <- data.frame(url="sample_A",
   type="coverage_gr",
   sample_id="sample_A");
sh1 <- prepareSashimi(GRangesList(TestGene1=test_exon_gr),
   filesDF=filesDF,
   gene="TestGene1",
   covGR=test_cov_gr,
   juncGR=test_junc_gr);
plotSashimi(sh1);

}
\seealso{
Other jam RNA-seq functions: 
\code{\link{assignGRLexonNames}()},
\code{\link{closestExonToJunctions}()},
\code{\link{combineGRcoverage}()},
\code{\link{defineDetectedTx}()},
\code{\link{detectedTxInfo}()},
\code{\link{exoncov2polygon}()},
\code{\link{flattenExonsBy}()},
\code{\link{getGRcoverageFromBw}()},
\code{\link{groups2contrasts}()},
\code{\link{internal_junc_score}()},
\code{\link{makeTx2geneFromGtf}()},
\code{\link{make_ref2compressed}()},
\code{\link{runDiffSplice}()},
\code{\link{sortSamples}()},
\code{\link{spliceGR2junctionDF}()}

Other jam plot functions: 
\code{\link{bgaPlotly3d}()},
\code{\link{factor2label}()},
\code{\link{gene2gg}()},
\code{\link{grl2df}()},
\code{\link{jitter_norm}()},
\code{\link{plotSashimi}()},
\code{\link{stackJunctions}()}

Other splicejam core functions: 
\code{\link{exoncov2polygon}()},
\code{\link{gene2gg}()},
\code{\link{grl2df}()},
\code{\link{make_ref2compressed}()},
\code{\link{plotSashimi}()}
}
\concept{jam RNA-seq functions}
\concept{jam plot functions}
\concept{splicejam core functions}
