## TODO for jambio

### Sashimi plot functionality

* Need ability to zoom without polygon clipping issues in ggplot2. Ideas:

    * Provide zoom in `prepareSashimi()` based upon GRanges, or range of
    exons. However in testing it still displays junctions whose
   ends are in range, extending outside of the intended range.
   * Exploit ggforce::facet_zoom() that somehow manages zoom without
   clipping issues. Does require the facet stage of rendering? It
   removes all previous faceting, which is problematic.
   * `ggforce::facet_zoom()` keeps all labels that use
   `ggrepel::geom_text_repel()`, which makes no sense. Hide
   labels outside the plot range.
   
* Transformed axis:

    * customize breaks, minor_breaks, and labels functions so
    ggplot2 labels include as many non-overlapping x-axis labels
    at exon boundaries as possible. Make minor_breaks include
    tick marks and optional smaller labels when zoomed in, or
    when an exon is wide enough to have an intermediate label.
    * Transformed axis line: x-axis line indicating compresssed
    status, solid line = uncompressed, dotted line = compressed.
    It might need a new own geom.

* Junction arc:

    * Remove the seam in the middle of each junction arc.
    * Make arcs correctly center during x-axis transforms.
    * Consider a custom geom "geom_diagonal_arc_wide()",
    effectively a two-sided "ggforce::geom_diagonal_wide()".

### Portability

* Verify handling of gzip files across architectures, specifically for
makeTx2geneFromGtf(). Fallback plan is to exit gracefully with an informative
message. It appears that data.table::fread() can accept a command as input,
in which case it runs the command and imports the data stream.
* Coverage data on Windows platform:

    * Can Windows users load coverage data? Package "rtracklayer"
    cannot import bigWig format.
    * What is Plan B for coverage data?

### Example data

* Define a small subset of transcript gene models to use in examples
for each function for improved documentation. **COMPLETE**
* Add "mouse_codon_usage.txt" raw text file as example
input for codon usage analysis. **COMPLETE**
* Small bigWig file?

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
