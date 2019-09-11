## TODO for splicejam

## Bug fixes

* `prepareSashimi()` throws an error when `flatExonsByGene`
contains different seqlevels than the `BigWigFile`, usually
when `flatExonsByGene` contains more seqlevels than present
in the `BigWigFile`. It happens even when the exon features
of interest involve seqlevels present in `BigWigFile`.
The fix is to reduce the seqlevels in `flatExonsByGene`
to match those in the `BigWigFile`, or to reduce seqlevels
to actual features in `flatExonsByGene`. At that point
an error would accurately reflect that `BigWigFile` does
not contain seqlevels for the requested seqlevels, and
thus would be helpful.


### Sashimi usability enhancements

* Separate the function `sashimiAppConstants()` to its own
proper standalone function, that prepares all the dependency
data objects like `flatExonsByGene`, `flatExonsByTx`,
`tx2geneDF`.
* New function to assemble `filesDF`. Essentially appends
a file to an existing `data.frame`.
* New function `plotlySashimi()` that optionally includes
the gene-transcript-exon model in the visualization.
* Zoom by exon range, using exon names. For example,
`zoom_by=c("Gria1_exon11", "Gria1_exon16")`.

### Sashimi longer term enhancements

* Enable display of sample replicates. Actually, enable
the grouping of `sample_id` entries. Default is to overlay
replicates in the same facet panel, facet by `sample_group`.
Alternative is to offset the y-axis similar to the `"ggridges"`
package. Unclear if the layers coverage polygons, and
junction arcs, will be out of sync; that is sample_1
should draw coverage and junctions, before sample_2 is
drawn.
* Enable BAM file for coverage data. Defer filtering rules
to the `Rsamtools` package, but show recommended examples
for properly paired reads.
* Enable BAM file for junction data. Determine recommended
filters for properly paired reads with minimum overlap width
at each end of the junction.
* Enable Sashimi plots over a region, not just one gene range.
Refactor of `prepareSashimi()` to take GRanges range
an optional coordinate range, or multiple genes (which would
imply the coordinate range if on the same chromosome).
* Create an example showing two genes with overlapping
exons.
* Mock up an example that combines gene exons, ChIP-seq peaks,
and ATAC-seq peaks, together into one schematic of a
genome region. It should show before and after, compressing
the "intervening sequences" between features.
* Mock up example showing `detectedTxInfo()` output in
the form of plots. For example, show percent max expression
as a heatmap, each call labeled to show the counts, TPM,
percentage, and whether it was called "detected" by
the thresholds.

### Farris R-shiny feedback

* Recalculate scale factors, use reference genes from in situ
hybridizations: Calm1, Camk2a, Camk2b, Actb, which show little to
no change.
* Plotly two-column layout does not keep all columns of panels
in sync when the x-axis is zoomed, the desired effect is to have the
x-axis of each sashimi panel zoomed consistently, across however
many columns are used for layout.
* Gria1 takes a while to display the first time, during which the app
is not usable. Consider making plot rendering asynchronous to decouple
the visualization from the use of the R-shiny app.
* The plotly exon coordinates show the overall coordinate range and
not the specific exon.

### Vignettes

* New vignette describing how to start a new R-shiny Sashimi App.

### Testing

* Implement methods in the testthat package to ensure a consistent
suite of tests to confirm full functionality with updates.
Highest priority functions, which would break parts of the workflow:
    * `flattenExonsBy()`
    * `assignGRLexonNames()` - dependent upon `jamba::makeNames()`
    * `annotateGRfromGR()` - dependent upon `shrinkMatrix()`
    * `shrinkMatrix()`
    * `combineGRcoverage()`
    * `df2colorSub()`

### data.frame versus tibble

* Evaluated tibble as possible replacement for data.frame in
`prepareSashimi()`. It shrinks data volume to about 7 times smaller,
reducing duplicated annotations per row for polygon coordinates.
However, ggplot2 still requires unnesting the tibble into a tall
format for plotting, it is unclear whether that step will incur
its own performance hit.
* ggplot2 polygon rendering should be faster -- possibly needs to
reduce the polygon detail prior to plotting, but unclear whether that
can happen since ggplot2 already chops polygons into multiple small
segments in order to handle transformed axis coordinates.

### R-shiny caching

* Review flatExonsByTx caching, which gets invalidated too frequently.
Consider using a fixed caching strategy, using defined keys instead
of the memoise default.

### R-shiny launchSashimiApp()

* Better error message for gene symbol not recognized.
* Plotly tooltips are so finicky, but try anyway:

    * coverage tooltip: y-position (coverage score), exon name (gr_name),
    track name (cov_name), preferably the coordinate range of the feature.
    * junction tooltip: junction score, name "gene_exonFrom gene_exonTo",
    and start-end coordinates of the junction.
    * plotly_onclick: show junction score or exon name in the middle
    of the respective polygon?
    
* Make font sizes more configurable.
* Add gene symbol somewhere, otherwise without the gene model shown,
there is no indication what gene is being displayed.
* Add a vignette describing how to set up an R-shiny splicejam.

### Sashimi plot functionality

* refactor prepareSashimi() to create full merged data.frame, thus
removing it from plotSashimi(); it should make plot faster.
* geom_diagonal_wide_arc(): new argument calc_midpoint=TRUE will
calculate the appropriate middle x-position using the active
x-axis transformation so the mid-point is the visible middle.
Optionally specify the fraction, 0.5 is middle, 0.25 is the first
quartile near the left edge, etc.
* Transformed axis:

    * customize breaks, minor_breaks, and labels functions so
    ggplot2 labels include as many non-overlapping x-axis labels
    at exon boundaries as possible. Make minor_breaks include
    tick marks and optional smaller labels when zoomed in, or
    when an exon is wide enough to have an intermediate label.
    * Transformed axis line: x-axis line indicating compresssed
    status, solid line = uncompressed, dotted line = compressed.
    It might need a new own geom.

### Portability/Usability

* Verify "Plan B for Windows users." Test and confirm whether
using Windows Subsystem for Linux (WSL) will enable bigWig import.
* Expand input types to allow BAM files, from which coverage will be
calculated. Will need an appropriate caching strategy, probably
at the level of covGR GRanges with RleList to store coverage
in value columns.
* Expand junction input to allow BAM files, from which junctions
can be inferred. (Check existing Bioconductor packages for
any that convert BAM to junctions.)

### Vignettes

* More Sashimi examples:

    * Using baseline to "shift" a set of exons up on the y-axis,
    use Ntrk2 as an example.
    * Use baseline with Gria1 to shift up/down the two alternate
    exons. (The original vision was something like a Sankey plot.)

* ALE-specific RNA-seq analysis workflow.
* Codon usage analysis workflow.

