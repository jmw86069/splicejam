# TODO for splicejam

## 12jul2021

* New feature (minor): Allow resizing the gene-transcript panel

   * Use case: Sometimes there are way too many transcripts, the
   panel should be much taller. Sometimes the exon labels are clipped
   off the bottom of the figure.
   * Note that adjusting "per panel size" affects gene-transcript
   panel, making it much too crowded.

* New feature (minor): Allow down-sampling the coverage profile

   * Use case: exons are rarely displayed one base per pixel, more
   often an exon is roughly 100-500 bases wide, but displayed in
   the equivalent of 25-50 pixels. The adjustment would be to
   down-sample data to help speed up the rendering step.

* New feature (minor): Optionally display gene model per column.

   * Use case: Currently with two-column display the gene model
   does not visually align to the coverages.

* New feature (inquiry): Test porting track type to ggbio or gviz.

   * Use case: To integrate other track types, coverage, peaks, genes,
   etc. the ggbio and gviz R packages are more feature-rich.
   * General idea is to provide two things: track "geom", and
   x-axis compressed axis which effectively compresses the intron ranges.
   * For Gviz, understand `GdObject-class` and create CustomTrack,
   which should mimic the workflow used by `AlignmentsTrack()`:
   
      ```R
      Gviz::CustomTrack(plottingFunction=function(GdObject, prepare=FALSE, ...){},
         variables=list(), name="CustomTrack", ...)
      ```

* New feature (minor): Color picker beside Shiny app Sample Selection

   * Use case is to display, and allow selection of colors per sample
   * Not sure if "easy" color picker widgets will work with
   the selectable table widget.


* New feature (major): multiple features in transcript window

   * Use case is to display more than one gene, or gene(s) alongside
   another track such as ChIP-seq peaks, or enhancer regions.
   * The main change would be to define a region by coordinates
   rather than by the gene symbol - which itself only implies
   coordinates.
   * All features displayed in the region would be used to define
   compressed exon ranges.
   * Major change: splice junction strandedness should be allowed
   to differ from the feature of interest, which may cause problems
   when using STAR junctions where strandedness is not always
   reliable.

* New features (major): coverage and junctions from BAM files

   * Use case is to allow BAM input, both for coverage and for junctions.
   * Major enabling feature: This step is a precursor to handling
   scRNA-seq, which could dynamically split the BAM reads by
   tSNE/UMAP clusters.

* New data (major): pre-defined hg19 transcriptome data required for sashimi plots

   * Gencode comprehensive with all genes, flat exons pre-computed.
   * Obvious enhancement: hg38, mm10, other common organisms.
   * Describe steps used so others can use their own GTF as needed.


* New feature (nice to have): ability to select/deselect transcripts
displayed in the R-shiny app for a given gene.

   Use case is when displaying all transcripts for a given gene,
   but there are too many un-expressed transcripts, e.g. human ACTB.
   Could allow a separate table to select transcripts to include/hide.
   When performing "Update" the compressed exons would be re-calculated
   for that gene using the selected exons.

## 19mar2021 todo items

* COMPLETE: Fix broken Sample Selection - change to table selection/ordering

   * include `sessionInfo()` in a hidden section on the R-shiny app
   to be able to confirm all packages and version numbers being used.

## 18mar2021 todo items

* COMPLETE: Fix the ggrepel exon labeling which is now above instead of below
the gene-exon model.
* for interactive plots, highlighting a splice junction also highlights
all splice junctions for the same exon span, which is good

   * for flat gene-exon models, name the gap using the format
   `geneSymbolexon_name1-geneSymbolexon_name2` which will allow
   the gap to be highlighted as well.
   * for flat transcript-exon models, not sure if this is
   appropriate.



## 23apr2020 todo items

* Consider adding transcript selection to the left menu items
* Workflow:

   * Idea is to allow manually selecting subset or superset of
   transcripts to include in the flatExonsByGene, partly to 
   allow removing transcripts that mess up the overall gene model,
   partly to allow showing a highlighted subset of transcripts
   for example for diffSplice/DEXseq splice hits. Sometimes a
   gene has 8 detected isoforms, but the predominant change
   involves only 2 or 3 of those isoforms. Would be great to
   be able to create the ideal subset.
   * (Note this workflow is already possible with manual
   steps outside the R-shiny app, but who wants to do that?)

* Allow expanding the x-axis genomic region being displayed, beyond
the gene body. Note this step requires defining the compression
to use for upstream region -- for example should it be 10:1
compressed in visual space?
* Note the genomic coordinates are not being displayed below
certain genes (human GAS5) -- why not? It might be because
this gene has one contiguous exon due to some noisy unspliced
isoforms. The coordinate labeling function might be trying to
use only the outside coordinate. Could adjust that function to
use either the edge, or if two labels are separated by more then
5% the total exon width, display that label.

   * Genomic axis label logic for compressed GRanges coordinates:
   Start by labeling outer edge of each exon, then internal exon
   boundaries. Add each label if it will be at least 5% distant
   from another label, based upon the total exon width, and total
   number of gaps.
   * Review new ggplot2 axis labeling rules and whether we can rely
   upon ggplot2 to hide axis labels that would otherwise overlap.

* Allow manually setting one common y-axis range. Currently only
possible to adjust the y-axis by using interactive plotly, which
loses some important junction labeling.
* Allow displaying any gene-exon model that sits inside the viewing
range, in case there is an interesting internal feature such as a miRNA.
This enhancement is the predecessor for supplying a genomic region
and having the sashimi plot generated in that region regardless of
genes involved. Note that junction strandedness will need to be
more fluid -- currently the junction strand is fixed consistent with
the target gene, mostly because STAR junction strandedness is
sometimes mis-annotated (i.e. 100% coverage on "+" strand, but one
junction might be labeled "-" strand.)
* Pie in the sky idea: some way to save per-gene settings, like
coordinate range, axis ranges, etc. So when someone goes to "Gria1"
they get the best default Gria1 view they can get. When they go to
"Smarca4" they get the best default Smarca4 view they can get.


## 21apr2020 todo items

* COMPLETE: It looks like junction stacking on negative strand is applied
to the wrong edge of each junction (left-right instead of right-left).
* COMPLETE: Allow the gene-transcript panel to be adjusted taller.

   * By default, gene panel height should be auto-adjusted based
   upon the number of transcripts shown.
   * Ideally the panel height should be adjusted based upon label height.

* Enable more default settings such as number of columns, panel height,
font sizes, etc. The app should be fully customizable upon startup,
so that new users get the intended default experience.

   * Available:
   
      * min_junction_reads - minimum junction read threshold
      * use_exon_names - one of `c("coordinates", "exon names")`
      * exon_range_selected - two gene_nameExon values
      * exon_range_choices_default - exon labels available for default gene
      * layout_ncol - number of columns in layout
      * include_strand - one of `c("+","-","both")` strands to display
      * share_y_axis - one of `c(TRUE, FALSE)` whether to use shared y-axis range
      
   * Todo:
   
      * panel height
      * show gene model
      * show transcripts
      * detected only
      * font size,
      * exon font size
      * interactive plot
      * minimum arc height? new option

* Consider minimum junction threshold based upon % coverage range
instead of absolute number. Not urgent, can be manually adjusted.
* Consider option to label gene exons by major exon number (exon1, exon2, etc.)
not using the disjoint exon numbers (exon1a, exon1b, exon2a, exon2b, etc.)
* `detectedTxInfo()`

   * option to return only detected results.
   * change labels to be user-friendly
   
* new function to convert `detectedTxInfo()` into ComplexHeatmap:

   * box around cells which meet the "detected" thresholds overall (cell_fun)
   * color cells by percent max expression per sample (default), option to use count, TPM
   * label cells: pseudo-counts (percent max expr, TPM abundance); 67.9 (13%, 1.2 TPM)
   * option to display "detected" or "all" Tx

* Consider optional per-panel label to supplement/replace the text
label used by ggplot2 facet wrap. Label could be positioned topright,
topleft, bottomright, bottomleft, or auto.
* Two-column displays are not visually synchronized with the gene
view, consider making a two-column gene view?
* When using plotly interactive plots with two columns, consider
syncing the two columns when zooming in one side, so the other column
is also zoomed into the same region.
* Consider adding two-column output selection to the main page.
* When viewing very large genes in plotly interactive mode, it
can take ages to render (after the ggplot object is already created,
it appears the slow step is creating all the plotly-specific javascript
to render every coordinate position, even though the region is zoomed
out.) Needs some optimization step to reduce the number of actual
x-axis points being displayed. SMARCA4 is an example gene, spans 150k
bases, yet every exon pixel is stored as a coordinate for plotly.
Consider reducing the resolution of x-axis points being rendered by
plotly somehow.

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

