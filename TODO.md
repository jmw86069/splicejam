## TODO for jambio

### R-shiny launchSashimiApp()

* Better error message for gene symbol not recognized.
* Optionally show transcript exon models.
* Optionally (new default) get coverage only for the gene strand.
    * Requires bigWig files are named with indication
    of strand. Current method of reading coverage to
    test for negative coverage values is still a time step.
    * Another option is to add strandedness to filesDF, but
    would only apply to bigWig files, not junctions.
* Color splice junctions using alternating lighter/darker
shading, to try to help visibility for crowded ribbons.

### Sashimi plot functionality

* Transformed axis:

    * customize breaks, minor_breaks, and labels functions so
    ggplot2 labels include as many non-overlapping x-axis labels
    at exon boundaries as possible. Make minor_breaks include
    tick marks and optional smaller labels when zoomed in, or
    when an exon is wide enough to have an intermediate label.
    * Transformed axis line: x-axis line indicating compresssed
    status, solid line = uncompressed, dotted line = compressed.
    It might need a new own geom.

### Portability

* Verify handling of gzip files across architectures, specifically for
makeTx2geneFromGtf(). Fallback plan is to exit gracefully with an informative
message. It appears that data.table::fread() can accept a command as input,
in which case it runs the command and imports the data stream.
* Coverage data on Windows platform:

    * Can Windows users load coverage data? Package "rtracklayer"
    cannot import bigWig format.
    * What is Plan B for coverage data?

### Testing

* include testthat into the package development workflow to help
test and confirm specific input and output criteria.

### Vignettes

* More Sashimi examples:

    * Using baseline to "shift" a set of exons up on the y-axis,
    use Ntrk2 as an example.
    * Use baseline with Gria1 to shift up/down the two alternate
    exons. (The original vision was something like a Sankey plot.)

* ALE-specific RNA-seq analysis workflow.
* Codon usage analysis workflow.

### Feature creep (testing new ideas):

* Test other labeling methods:

    * geom_mark_hull() worked for exon polygons, not ideal label placement
    * geom_text() direct label placement, no overlap management
    * geom_text_repel() need method to hide labels outside zoom range

* Sashimi panels:

    * Add flattened exon model at the baseline position for each exon,
    perhaps using the GRanges used by `make_ref2compresssed()` as
    an x-axis alternative?
