% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jambio-tx.R
\name{annotateGRfromGR}
\alias{annotateGRfromGR}
\title{Annotate GRanges using another GRanges object}
\usage{
annotateGRfromGR(
  GR1,
  GR2,
  grOL = NULL,
  numAsStrings = FALSE,
  stringShrinkFunc = function(...) {
     jamba::cPasteUnique(..., doSort = TRUE)
 },
  numShrinkFunc = sum,
  addStringCols = NULL,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE,
  select = "all",
  sep = ",",
  verbose = FALSE,
  DEBUG = FALSE,
  ...
)
}
\arguments{
\item{GR1}{GRanges object, the reference object to which annotations
are added.}

\item{GR2}{GRanges object used for annotations}

\item{grOL}{overlaps object, optionally used when the
\code{GenomicRanges::findOverlaps()} function has already been run, mostly
intended to save re-processing the overlaps for large objects.}

\item{numAsStrings}{logical indicating whether numerical values should
be treated as strings for the purpose of added annotations. When
\code{TRUE}, numerical values will be converted to character strings and
retained as if they were labels. This argument treats all numeric
columns as string, to treat only a subset use the \code{addStringCols}
argument.}

\item{stringShrinkFunc}{function used to shrink string annotations
by overlapping feature. This function should take a list as input,
and \code{sep} as an argument indicating the delimiter if applicable,
and return a character vector of the same length as the list. By
default, \code{jamba::cPasteUnique()} is used. Control over the sort
and uniqueness should be applied with this function.
Note: \code{stringShrinkFunc} can be a single function, or can be a
named list of function calls, whose names match the colnames
of \code{values}. If a named list is provided, the first entry is used
for any any names in \code{values} which are not in \code{names(stringShrinkFunc)}.}

\item{numShrinkFunc}{function used to shrink numerical values
by overlapping feature. For numeric values, the \code{data.table} package
is used to apply the function, which enables custom functions not
otherwise available via the \code{GenomicRanges} methods. Therefore, the
function does not take a list as input, instead takes a numeric
vector as input.
Note: \code{numShrinkFunc} can be a single function, or can be a
named list of function calls, whose names match the colnames
of \code{values}. If a named list is provided, the first entry is used
for any any names in \code{values} which are not in \code{names(numShrinkFunc)}.}

\item{addStringCols}{character vector, optional, of numeric colnames that
should be treated as string values. This option is a subset of the
\code{numAsStrings}.}

\item{type, ignore.strand, select}{arguments sent to
\code{GenomicRanges::findOverlaps()}. Note these options are only used
when \code{grOL} is not supplied.}

\item{sep}{character value indicating the delimiter to use between
string values.}

\item{verbose}{logical indicating whether to print verbose output.}

\item{DEBUG}{logical indicating whether to print debugging output.}

\item{...}{additional arguments are ignored.}
}
\value{
GRanges object with colnames added to \code{values}, with length
and order equal to the input \code{GR1} GRanges object.
}
\description{
Annotate GRanges using another GRanges object
}
\details{
This function adds annotations to features in the given GRanges
object, from overlapping features in a second GRanges object.
It is especially useful after performing a manipulation that
results in a GRanges object without any \code{values} annotation,
for example \code{GenomicRanges::reduce()} or \code{GenomicRanges::intersect()}.

In theory this function is relatively simple, it applies annotations
to overlapping entries. In practice, it gets complicated when multiple
annotations overlap the same GRange entry. In that case, numerical
values by default return the \code{base::mean()} while string values call
\code{jamba::cPasteUnique()} by default, and return the unique, sorted,
comma-delimited set of string values. For example, overlapping several
exons from the same gene, the resulting annotation might include
just one gene symbol.

The numeric values can be aggregated using another function, controlled
by \code{numShrinkFunction} named using the colname. For example:
\code{numShrinkFunction="min"} would return the \code{base::min()} value for all
numeric columns. But \code{numShrinkFunction=list(default=mean, score=max)}
would return the \code{base::mean()} for all numeric columns except the
\code{"score"} column which would report the \code{base::max()}.

Numeric values currently use the \code{data.table} package workflow, which
provides extremely efficient calculations for numeric values. However,
vectorized functions have the potential to be notably faster, as is
true with \code{jamba::cPasteUnique()} especially when applying
\code{jamba::mixedSort()} to sort alphanumeric values. In future, the
implementation may change to accomodate vectorized numeric functions,
instead of using \code{data.table}.

TODO: This function is a specific case where another function
\code{shrinkDataFrame()} may be a more general purpose use. That function
is yet to be ported to the Jam package suite.
}
\examples{
gr12 <- GenomicRanges::GRanges(
   seqnames=rep(c("chr1", "chr2", "chr1"), c(3,3,3)),
   ranges=IRanges::IRanges(
      start=c(100, 200, 400, 500, 300, 100, 200, 400, 600),
      width=c(100,150,50, 50,50,100, 50,200,50)
   ),
   strand=rep(c("+", "-", "+"), c(3,3,3)),
   gene_name=rep(c("GeneA", "GeneB", "GeneC"), each=3)
)
gr1 <- gr12[,0];

# Say for example you have a GRanges object with no annotation
gr1;

# And had another GRanges object with annotations
gr12;

# To add annotations
annotateGRfromGR(gr1, gr12);

# Notice certain features overlap multiple annotations,
# which may be acceptable.

# If you want to keep annotations distinct,
# use annotateGRLfromGRL()
grl1 <- GenomicRanges::split(gr12[,0],
   GenomicRanges::values(gr12)$gene_name);
grl2 <- GenomicRanges::split(gr12,
   GenomicRanges::values(gr12)$gene_name);

# The first object is a GRangesList with no annotations
grl1;

# The second object is a GRangesList with annotation,
# assumed to be in the same order
grl2;

annotateGRLfromGRL(grl1, grl2);

}
\seealso{
Other jam GRanges functions: 
\code{\link{addGRLgaps}()},
\code{\link{addGRgaps}()},
\code{\link{annotateGRLfromGRL}()},
\code{\link{assignGRLexonNames}()},
\code{\link{closestExonToJunctions}()},
\code{\link{combineGRcoverage}()},
\code{\link{exoncov2polygon}()},
\code{\link{findOverlapsGRL}()},
\code{\link{flattenExonsBy}()},
\code{\link{getFirstStrandedFromGRL}()},
\code{\link{getGRLgaps}()},
\code{\link{getGRcoverageFromBw}()},
\code{\link{getGRgaps}()},
\code{\link{grl2df}()},
\code{\link{jam_isDisjoint}()},
\code{\link{make_ref2compressed}()},
\code{\link{sortGRL}()},
\code{\link{spliceGR2junctionDF}()},
\code{\link{stackJunctions}()}
}
\concept{jam GRanges functions}
