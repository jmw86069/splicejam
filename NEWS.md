# splicejam version 0.0.34.900

## changes

* R-shiny new options: minimum junctions; show by strand.
* R-shiny now defines each to the global environment if they do
not already exists: exonsByTx, cdsByTx, flatExonsByGene, flatExonsByTx,
tx2geneDF, detectedTx, detectedGenes. This mechanism is used for now,
both to help define custom input data, but also to help define the
data values needed to produce plots manually outside R-shiny.

# splicejam version 0.0.33.900

## changes

* `stackJunctions()` now also returns "junction_rank"
scored from 1 to 3, where 1=minor, 2=partial, 3=major,
based upon the rank of junctions entering and leaving
exons.
* Junctions are drawn in order from dominant to minor,
which has the effect of ensuring smaller junctions are
displayed. The `junc_fill` is converted to gradient so
a darker color is used for dominant junctions, to try
to highlight the major isoform splice junctions per sample.
Unfortunately, plotly does not honor the polygon render
order (yet).

# splicejam version 0.0.32.900

## changes/fixes

* Added `flatExonsByTx` to R-shiny app, to be able to display
transcript isoform exon structures alongside the flattened
gene-exon model.
* Added R-shiny option to display transcript isoforms alongside
the gene-exon model. Also option to show all or detected
transcripts.

# splicejam version 0.0.31.900

## changes/fixes

* Splice junctions now require at least one junction end within
the range being displayed, to avoid junctions that do not connect
with the visible gene model.
* The x-axis range of multiple plots are more consistently controlled,
to avoid one panel autoscaling to display a wider junction than
other panels.

## additions

* Plotly views contain custom tooltip text, using the
sample_id, the name of the exon or feature, and the track
(referring to the name of the coverage).

# splicejam version 0.0.30.900

## additions

* New ggplot2 geom `geom_diagonal_wide_arc()`, extends
`ggforce::geom_diagonal_wide()` by connecting the left and right halves
of a diagonal into one continuous ribbon arc.

## changes

* `plotSashimi()` now uses `geom_diagonal_wide_arc()` to display
junction arcs instead of using two `ggforce::geom_diagonal_wide()`.
For now, `prepareSashimi()` still provides coordinates to create
two polygons, which are combined in the geom. Still todo: figure
out how to add a label; and make the middle x position immune to
x-axis rescaling.

# splicejam version 0.0.29.900

* R-shiny app now properly keeps interactive plot x-axis ranges
in sync, when zooming the plot.
* Fixed bug in `prepareSashimi()` that occurred when
no junctions overlapped another, resulting in error
"invalid 'type' (S4) of argument".

# splicejam version 0.0.28.900

## changes to existing functions

* `defineDetectedTx()` has a new argument `zeroAsNA`
to handle the special case where some expression values
reported as zero should be treated as `NA` (no data obtained)
and therefore will not be included in group mean calculations.
As transcript-exon models get more "comprehensive", kmer
quantitation tools such as Salmon and Kallisto sometimes
need to assign abundance to one of several nearly identical
isoforms, and in low count scenarios all counts may be
assigned to one or another isoform, leaving zero in the
alternate position. Excluding zero ensures that group
mean values represent only the assigned quantitation.

# splicejam version 0.0.27.900

## bug fixes

* Packages now imported for `launchSashimiApp()`: `shiny`,
`shinydashboard`, `shinydashboardPlus`, `shinyWidgets`.
* `getGRcoverageFromBw()` now returns NULL when a BigWig file
is not accessible. Previously any error caused
`"Error in seqinfo(con) : UCSC library operation failed"` which
could mean the file does not exist, or any number of other validation
checks failed. Since `getGRcoverageFromBw()` used a vector of
BigWig files, any error caused the entire set to fail.
* `getFirstStrandedFromGRL()` added package prefix to the use
of `IRanges::heads()` which was not imported directly.

## changes

* R-shiny plots now set the plot height for ggplot2 or plotly,
depending upon the number of samples, and presence of
gene-exon model. Future versions will allow R-shiny user
customization.
* `bgaPlotly3d()` removed some unnecessary print functions.

# splicejam version 0.0.26.900

## changes

* This version is an incremental update, while we evaluate
some display options of Sashimi plots in the R-shiny app.
* Minor aesthetic changes to R-shiny app, including plotly
interactive options, and evaluating some UI elements in sidebar.
* Added option to display gene-exon model in R-shiny.

# splicejam version 0.0.25.900

## changes

* `makeTx2geneFromGtf()` help docs recommend installing the `R.utils`
package when `.gz` files are required for import, since this process
uses `data.table::fread()`.
* Package dependencies were added for shiny, shinydashboard.

## new functions

* `to_basic.GeomShape()` is a hidden function intended only to
provide basic support of `ggforce::geom_shape()` when converting
ggplot objects to plotly.

# splicejam version 0.0.24.900

## changes

* R-shiny Sashimi app now allows zooming into coordinate ranges,
defined per gene.

# splicejam version 0.0.23.900

## additions

* First working Sashimi R-shiny app, still has issues with
certain genes that needs debugging. Defaults to using
Farris et al data from `farrisdata` package, but can
be overridden with global environment variables.

# splicejam version 0.0.22.900

## additions

* Initial R-shiny `launchSashimiApp()` function.
* Suggests the `farrisdata` package for example data used
for the Farris et all RNA-seq mouse hippocampus manuscript.

# splicejam version 0.0.21.900

## bug fixes

* Fixed examples that did not explicitly call library() for
required libraries.
* Fixed error in `codonUsage2df()` example data file path.

## changes

* Added specific TODO items.
* Added two function categories: `"jam ALE-specific RNA-seq functions"`,
and `"jam codon-usage RNA-seq functions"` to help organize the
large list of accessory functions.


# splicejam version 0.0.20.900

## changes

* Minor fixes to examples.

## additions

* Added data `test_exon_gr`, `test_junc_gr`, `test_cov_gr` for
easy example data for exons, junctions, and coverage, respectively.
* Added "wide" variants of the above test data, with introns
about 100x larger than exons, consistent with mammalian gene
structures.
* Added examples to each data, showing how one would use the
raw data to generate different visualizations used in
Sashimi plots.
* Added sample data to the vignettes and some examples.
* Renamed `flattenExonsByGene()` to `flattenExonsBy()` to reflect
that the function handles `by="gene"` and `by="tx"`.

# splicejam version 0.0.19.900

## additions

* Added `create-a-sashimi-plot.Rmd` which walks through a full example
showing how to create a Sashimi plot.
* Added examples to `stackJunctions()` with schematics, including
example on plotting junctions by themselves.

## changes to existing functions

* Removed `compressGRgaps()` for now, since the methods now try to
keep GRanges intact, and instead transform the x-axis scale to compress
visible gaps.

# splicejam version 0.0.18.900

## changes to existing functions

* `bgaPlotly3d()` now properly handles ellipsoid colors,
previously the colors were assigned but not honored by
`plotly::add_trace()`.
* `combineGRcoverage()` determines strandedness by requiring all
values to be at or below zero, with at least one negative value.
Otherwise, data can have position and negative values and will
be considered positive stranded.
* `plotSashimi()` replaces `jamSashimi()` because it is just more
intuitive... Ah well.
* Fixed small typo in plotSashimi.Rd help that included an unmatched quote.

# splicejam version 0.0.17.900

## new functions

* `jamSashimi()` is used to plot Sashimi data prepared by
`prepareSashimi()` in order to separate the download
and preparation of Sashimi data from the visualization.

## changes to existing functions

* `prepareSashimi()` is refactored to remove the plot functions,
moving them to the new `jamSashimi()`.
* `flattenExonsByGene()` can now handle Tx data, mainly useful
to add CDS regions to existing exon models.
* `gene2gg()` is more robust to edge input cases.

# splicejam version 0.0.16.900

## new functions

* `prepareSashimi()` is a workhorse function that prepares several
types of Sashimi-like data output, including ggplot2-based
Sashimi plots. This plot uses vanilla ggplot2 with custom x-axis
scaling to compress genomic coordinates, while keeping data in
proper genomic coordinate space.
* `gene2gg()` is a lightweight function that creates gene and
transcript exons models, and returns a ggplot2 object for
plotting. It can optionally return the data.frame used by ggplot2.
* `grl2df()` is similar to `fortify()` from ggplot2, or the
`broom::tidy()` functions, that take a custom object and return
a data.frame for downstream use. Here the function currently works
with "rectangle" data (like exons, introns, peaks, etc.), and
"junction" (like junction regions to be represented by
`ggforce::geom_diagonal_wide()`). In future it may also handle
sequence coverage polygons.
* `exoncov2polygon()` converts a GRanges object containing coverage
in the form of NumericList for each exon or intron, into a data.frame
suitable for use by `geom_polygon()` or `ggforce::geom_shape()`.
* `getGRcoverageFromBw()` loads bigWig coverage data for a GRanges of
exons and introns, returning a GRanges object with columns representing
each sample_id. It calls `combineGRcoverage()` to combine coverages
by strand within each `sample_id`.
* `flattenExonsByGene()` is intended to provide a set of unique exons
per gene, using all or only detected transcript exon models. It numbers
exons by contiguous segment, then labels subsections of each exon
with a letter suffix, for example "exon1", "exon2a", "exon2b", "exon3".
Lastly, it can annotate exons as CDS or non-CDS, if given a set
of `cdsByTx` data.
* `make_ref2compressed()` takes a GRanges object of exons (or any
useful feature like ChIP-seq peaks, etc) and returns functions needed
to compress x-axis coordinates, in order to shrink gaps/introns to
a fixed width, suitable for use by ggplot2.
* `spliceGR2junctionDF()` takes junctions GRanges data, and flatExonsByGene,
and summarizes and annotates each junction by the gene_exon connected,
and combines scores when multiple junctions are within "spliceBuffer"
distance of an exon boundary, by default 3 bases. It can accept a
GRanges object that was derived from multiple sample_id, and will return
data summarized within each sample_id. It calls `stackJunctions()` to
combine junctions per replicate.
* `simplifyXY()` takes a vector of coordinates, and simplifies them to
the minimal set of points to represent the polygon. For sequence
coverage data, that process can reduce individual points by 90%,
but it works with any coordinate data. When two or more consecutive line
segments have identical angle (with non-zero distance), they are combined
using only the first and last points.
* `getGRgaps()`, `getGRLgaps()`, `addGRgaps()`, `addGRLgaps()` are functions
that take GRanges input, and return either just the gaps between
non-overlapping regions, or the original features with gaps inserted between
them. Useful to convert a set of exons, to a set of exons and introns.

## other changes

* Added package dependencies to ggplot2, ggforce, ggrepel

## TODO

* Need a way of annotating GRangesList gaps by the parent feature name.
* In future, allow samples to have replicates, optionally allow each
replicate to be independently plotted, while being part of the
parent `sample_id`.
* Allow alternative input for `prepareSashimi()`, namely
BAM alignment files, from which coverage and junctions can be
derived.
* Consider making `ref2c` chromosome-wide, however it would require
x-axis coordinates to know their context, in terms of chromosome/seqname.
* Multiple vignettes to demonstrate workflows: Sashimi plots; preparing
exon GRangesList; using `annotateGRLfromGRL()` to add annotations, etc.

# splicejam version 0.0.15.900

## new functions

* `curateVtoDF()`, `curateDFtoDF()` are data curation functions
to curate a vector, or a data.frame, into a data.frame with
consistent, usable nomenclature for downstream analysis. They use
a flexible yaml format that should help automate analysis pipelines
that start with raw data file import.
* `exoncov2polygon()` and `compressPolygonM()` are basic functions
for sashimi plots, efficiently converting exon/intron coverages to
multi-polygons, then optionally compressing introns and reducing
the information content to roughly similar resolution as uncompressed
regions.
* `spliceGR2junctionDF()` is the summary function to convert a set
of splice junction read counts to gene-annotated junctions, grouped
and summed where needed.
* `closestExonToJunctions()` called by `spliceGR2junctionDF()` is used
to annotate junction ends near compatible exon boundaries, with
some buffer distance allowed to "snap" junctions to the edge.


# splicejam version 0.0.14.900

## bug fixes

* Fixed issue where `numTxs` was not getting populated in `runDiffSplice()`,
otherwise the stats summary is not changed.

# splicejam version 0.0.13.900

## changes

* Added "openxlsx" to Suggests, for exporting Rmarkdown stats tables
to Excel format.
* `runDiffSplice()` includes examples using `groups2contrasts()`, also
allows `txColname,geneColname` to be defined and therefore custom.
* Changed all `verbose=TRUE` to `verbose=FALSE` by default.
* Changed `makeTx2geneFromGtf()` to use `data.table::fread()` ability
to uncompress .gz files, hopefully making it cross-platform.
* `groups2contrasts()` was updated to handle single-factor experiments,
and be more robust to two-factor experiments with missing groups
in the full design table.

## additions

* Added basic RNA-seq workflow to vignettes.

# splicejam version 0.0.12.900

## new functions

* `runDiffSplice()` is a wrapped around `limma::diffSplice()` intended to
capture several steps of pre- and post-processing, optionally applying
`limma::voom()`.
* `groups2contrasts()` takes a vector or data.frame of sample groups,
and determines the relevant pairwise and two-way contrasts, returning
a design matrix and contrast matrix.
* `sortSamples()` sorts biological samples so that known patterns of control
group terms are returned before non-control groups, in order to help
provide useful defaults when setting up sample group contrasts.
* `strsplitOrdered()` provides `base::strsplit()` but returns factor output
whose factor levels are influenced either by factor input, or by other
arguments.

## changes

* `limma` package was added to Imports.
* `runDiffSplice()` has an argument `spliceTest` which defines `test`
when calling `limma::topSplice()`. Only the "t" (t-test) returns fold change,
therefore hits are not filtered by fold change for "F" or "simes".


# splicejam version 0.0.11.900

## new functions

* `ale2violin()` takes output from `tx2ale()` with some arguments, and
produces a ggplot2 violin plot object, as well as the underlying data.
It allows a custom filtering function, which allows filtering gene lists
for relevant regions of expression.

# splicejam version 0.0.10.900

## enhancements

* `shrinkMatrix()` minor fix to remove the shadow `x` in the resulting
colnames.
* `tx2ale()` modified to be tolerant of NA values in the expression matrix.

# splicejam version 0.0.9.900

## new functions

* `sortGRL()` sorts GRangesList objects by chromosome and position.
* `getFirstStrandedFromGRL()` return the first GRanges feature per
GRangesList, ordered properly by strand.
* `annotateGRfromGR()` which annotates one GRanges object based upon
overlaps with another GRanges object.
* `annotateGRLfromGRL()` which annotates one GRangesList object based upon
overlaps only with the same index entries in a second GRangesList object.
* `findOverlapsGRL()` runs `GenomicRanges::findOverlaps()` for the case of
two GRangesList objects, matching at the GRanges level but restricting
overlaps to those matching the same GRangesList index.
* `assignGRLexonNames()` assigns exon names and numbers to a disjoint
set of exons per gene model.

## updates

* Added "S4Vectors" to package imports, since some functions like
`lengths()` have been moved there. Another reason using a package prefix
is clunky at best. If methods have dispatch, and if they did not
conflict between S3 and S4 methods, that would be the better strategy.
It is not a solution to hardcode package names into function bodies,
since every package maintainer has to stay updated on the source package
of all other dependent functions. Those details should be irrelevant to
other package maintainers.


# splicejam version 0.0.8.900

## enhancements

* The "arules" package was moved to "Imports" since it is required
for the `list2im()` function. A slower workaround could be written,
but ultimately the arules package is preferred.
* Configured the package to use pkgdown for function references.
* Added "#' @imports jamba" to import all jamba R functions.

# splicejam version 0.0.7.900

## new functions

* `factor2label()` will convert a factor to a factor label, with the
same order of levels as the input factor, but including summary stats
like the number of items for each factor level. Useful for ggplot2
visualizations, to include counts in the color legend for example.

# splicejam version 0.0.6.900

## enhancements

* Added reference file "Mouse_codon_usage.txt" to the "extdata"
folder, as example input both for farrisSeq.Rmd, and for other
species. Access the file path with this command:
`system.file("extdata", "Mouse_codon_usage.txt", package="farrisdata")`

## new functions

* `codonUsage2df()` imports a text codon usage file, and returns a
data.frame. Support function `dna2codon()` simply takes a character
vector and returns vector whose elements all have three characters,
e.g. all(nchar(x) == 3).
* Two geometric mean functions: `jamGeomean()` is preferred by
the Jam packages, but `geomean()` is provided for direct comparison
to the classical approach.
* `detectedTxInfo()` summarizes the data used to define detected
transcripts for a given gene, or for a given set of transcripts.

# splicejam version 0.0.5.900

## bug fixes

* Minor fix to `makeTx2geneFromGtf()` to remove requirement for
`col.names` during import.
* Several minor documentation updates.

# splicejam version 0.0.4.900

## bug fixes

* Updated the DESCRIPTION file to include proper "Remotes" entries
and depencies for jamba, colorjam, and jamma packages.
* Fixed issue where the `data.table` package required a specific
`#' @import data.table` entry in the roxygen2 entry for the
`shrinkMatrix` function. This issue prevented `defineDetectedTx()`
from working within an R package.

## new RNA-seq functions

* `defineDetectedTx()` is a core RNA-seq function to determine
which transcript isoforms are considered "detected" based upon
counts (or pseudocounts), TPM values if available, and the
relative abundance of isoforms per gene. When "everything with
over 10 counts" is not sufficient.
* `tx2ale()` takes a set of transcript exon models, and returns
a set of alternate 3'UTR regions, similar to ALE (alternative
last exon.) It differs slightly from other methods in that it
aggregates transcript quantities by shared ALE regions,
producing a matrix based upon unique ALEs. It does not use
counts in the ALE region itself, but the aggregate abundance
of all isoforms that contain each ALE. This method allows
tools like Salmon or Kallisto to quantify isoforms using the
best available kmers per isoform, without being restricted
to the ALE regions which have a very wide distribution of
lengths.
* `makeTx2geneFromGtf()` takes a GTF file, and produces a
`data.frame` with tx (transcript) to gene associations.

## new functions

* `list2im()` converts a list to an incidence matrix. It uses
the `arules` package fast methods for creating
`transactions` objects, which are sparse binary matrices with
linked data.frames used to describe rows and columns. Similar
to `SummarizedExperiment` but with optimization specifically for
list-to-matrix and matrix-to-list conversion. Their implementation
is memory efficient as well.
* `bgaPlotly3d()` takes a `"bga"` class, as produced by the
`made4::bga()` function, and produces a 3-D plotly visualization.
It adds ability to group centroids into "supergroups" which are
connected by a spline 3-D curve.
* `spline3d()` is the supporting function to draw interpolated
3-D curves through a set of coordinates.
* `shrinkMatrix()` shrinks a numeric matrix by groups of rows,
essentially a wrapper function for the `data.table` package
functions.
