# splicejam version 0.0.17.900

## new functions

* `jamSashimi()` is used to plot Sashimi data prepared by
`prepareSashimi()` in order to separate the download
and preparation of Sashimi data from the visualization.

## changes to existing functions

* `prepareSashimi()` is refactored to remove the plot functions,
moving them to the new `jamSashimi()`.

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
