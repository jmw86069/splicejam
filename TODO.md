## TODO for jambio

### Sashimi plot functionality

* Needs ability to zoom without polygon clipping issues in ggplot2.
Suggestions: prepare only subset of exons or a fixed chromosome region
during `prepSashimi()`; however it still displays any junction whose
ends are in this range, sometimes extending far outside of the intended
range.
* Customize the transformed axis breaks, minor_breaks, and labels
functions so ggplot2 labels include as many non-overlapping x-axis
labels at exon boundaries as possible. Make minor_breaks include
tick marks and optional smaller labels when zoomed in, or when an
exon is wide enough to have an intermediate label.

### Portability

* Verify handling of gzip files across architectures, specifically for
makeTx2geneFromGtf(). Fallback plan is to exit gracefully with an informative
message. It appears that data.table::fread() can accept a command as input,
in which case it runs the command and imports the data stream.

### Example data

* Define a small subset of transcript gene models to use in examples
for each function for improved documentation.
* Add "mouse_codon_usage.txt" raw text file as example
input for codon usage analysis.

### Testing

* include testthat into the package development workflow to help
test and confirm specific input and output criteria.
