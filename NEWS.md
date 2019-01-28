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
