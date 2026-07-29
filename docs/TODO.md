# TODO for splicejam

## 28jul2026

- Consider ability to use custom transcript labels, for example “Gria1
  flip” and “Gria1 flop”.

- Consider ability to provide `detectedTx` to an existing splicejam
  environment. It would rebuild `flatExonsByGene` using the new subset
  of transcripts.

  - It is already possible by re-creating the `environment` with the new
    `detectedTx` data.

- Bigger todo ideas:

  - Give differential isoform results, create the figure to help
    visualize the supporting data. E.g.
    [`limma::diffSplice()`](https://rdrr.io/pkg/limma/man/diffSplice.html)
    gives a set of significant ‘transcript_id’, use those to create the
    figure.
  - Given coverage data, determine whether it is feasible to determine
    “detected transcripts”, outside of using Salmon transcript quant
    output for example.
  - Port
    [`defineDetectedTx()`](https://jmw86069.github.io/splicejam/reference/defineDetectedTx.md)
    to kallisto input.

## 27jul2026

- Provide R-shiny method to customize y-axis ranges per panel. Probably
  re-use the Samples table since it saves adding a new UI element which
  would essentially be the same table.

- Provide R-shiny option to edit `scale_factor`, same approach used for
  y-axis limits.

- DONE. Provide a clear way to adjust the y-axis label font size, for
  example the counts, and separately the gene/transcript labels.
  Controlled with `base_size` while gene axis labels use `geneAxisSize`.

- FIXED. The gene/transcript labels are intended to be angled slightly
  down.

- Consider option to customize the x-axis tick marks, labels, label
  density, etc. Custom xlab?

- Re-enable the option to hide the gene-exon model. To be fair, it was
  rarely used.

- Write “importer” functions to use Bioconductor packages as input.

  - Alternative to create sashimi environment, which does not require a
    GTF file.

  - Organism transcript packages such as knownGenes `'TxDb'`, or EnsEMBL
    data via `'EnsDb'`.

  - Note ‘TxDb’ requires organism annotation such as `'org.Mm.eg.db'`
    for the relevant organism. The knownGenes packages use ENTREZID as
    the ‘gene_id’ by convention. It is required to create tx2geneDF.

  - Note ‘EnsDb’ packages have gene-to-transcript relationships stored
    inside, no need for ‘org’ annotation data to create tx2geneDF.

  - ‘TxDb’ seem to use ENTREZID for gene_id.

  - ‘tx2geneDF’ can be created with ‘TxDb’ and anno.

  - ‘EnsDb’ has complete tx2gene and exonsBy support.

- INCOMPLETE. Allow editing ‘scale_factor’ in the Sample tab.

  - The Shiny widget exists but is clunky. Values are not yet utilized
    in
    [`sashimiAppServer()`](https://jmw86069.github.io/splicejam/reference/sashimiAppServer.md).

- TODO. Move a lot of functions to ‘keywords internal’ to minimize the
  user-facing function space.

## 21jul2026

- DONE. Adjust junction count label to be below negative strand ribbon,
  not above.

- DONE. Add testthat tests for
  [`splicejamFigure()`](https://jmw86069.github.io/splicejam/reference/splicejamFigure.md).

- DONE. Convert `exoncov2polygons()` to keep `NumericList`.

  - DONE. Also modify the custom ggplot2 stat to accept `NumericList`
    and add baseline y=0 to begin and end of the `geom_shape()` for each
    polygon.

- TODO.
  [`sashimiDataConstants()`](https://jmw86069.github.io/splicejam/reference/sashimiDataConstants.md)
  option to filter for `detectedGenes` and/or `detectedTx` to reduce
  data volume.

- DEFER. Create a strategy to handle overlapping genes, or one plot with
  one or more genes, e.g. a coordinate range.

## 14jul2026

- PARTIAL. Enhance plotly output from
  [`splicejamFigure()`](https://jmw86069.github.io/splicejam/reference/splicejamFigure.md)

  - DEFER. Enable highlighting, synchronize panels by ‘feature’.
    Previously called
    [`plotly::highlight_key()`](https://rdrr.io/pkg/plotly/man/highlight_key.html)
    which is actually `crosstalk::SharedData$new()` and which alters the
    data from `data.frame` to R6 object. Then
    [`plotly::highlight()`](https://rdrr.io/pkg/plotly/man/highlight.html).
    It fails via downstream issue.
  - Improve the padding at the top.
  - Improve the hover text for each type of feature. See plotly docs on
    how to customize ggplotly output, specifically how to customize
    hover text after creation. Junction should show the score, from-to
    exon label. Coverage should show running coverage at that position,
    exon name. Gene model should show the feature name and type.

## 06jul2026

- DONE. Update ‘README.Rmd’ and vignettes with more simplified workflow,
  using
  [`sashimiDataConstants()`](https://jmw86069.github.io/splicejam/reference/sashimiDataConstants.md).

- DONE. Consider replacing logic in
  [`launchSashimiApp()`](https://jmw86069.github.io/splicejam/reference/launchSashimiApp.md)
  with
  [`splicejamFigure()`](https://jmw86069.github.io/splicejam/reference/splicejamFigure.md).

- DONE. Transition to `progressr` for Shiny and CLI progress indicators.

- DONE. Add optional progress bar to
  [`splicejamFigure()`](https://jmw86069.github.io/splicejam/reference/splicejamFigure.md).

- DEFER. When applying ‘junc_color’ only for fill_scheme=‘exon’, also
  apply the light-to-dark shading as used with fill_scheme=‘sample_id’.

- TODO.
  [`gene2gg()`](https://jmw86069.github.io/splicejam/reference/gene2gg.md)
  option to set gene/transcript order.

- TODO. Method to validate the splicejam environment. flat exons,
  tx2geneDF, detectedGenes/detectedTx as needed.

- DONE. Update
  [`splicejamFigure()`](https://jmw86069.github.io/splicejam/reference/splicejamFigure.md)
  with ‘use_memoise=TRUE’ the same as with
  [`launchSashimiApp()`](https://jmw86069.github.io/splicejam/reference/launchSashimiApp.md).

- DONE. Refactor prepareSashimi(), plotSashimi() to accept ‘envir’ as
  input as alternative to multiple objects.

- Consider replacing ‘envir’ with proper S4/S7 object.

## 24jun2026

- DONE. Add some wrapper function to do the steps wrapped inside the
  Splicejam Shiny app, specifically range of exons or range of genome
  coordinates to display. See 14apr2026 note below.
- DONE.
  [`prepareSashimi()`](https://jmw86069.github.io/splicejam/reference/prepareSashimi.md)
  update to permit coordinate range input.
- DONE. Consider accepting `SJ.out` format directly, validate using
  junctionBed bed12-format from the same file.
- PARTIAL. Optimize
  [`exoncov2polygon()`](https://jmw86069.github.io/splicejam/reference/exoncov2polygon.md),
  which has become rate-limiting.
- PARTIAL. Consider refactoring splicejam ‘df’ format for polycon
  coordinates as list on each row. Edit: Coverage now stores the entire
  sequence on one row using AsIs with `NumericList` data direct from
  cpp11bigwig. Woot.
- DONE. Add visual unit tests (vdiffr) for sashimi plot output, to cover
  the range of intermediate steps altogether.
- DEFER. Add unit tests for all the various core processing steps.

## 16jun2026

- DONE. Fix bigwig load failures/slowness with remote URLs.

## 30apr2026

- DEFER. Make `GRanges` utility functions more prominent for wider
  usefulness.

  - [`annotateGRfromGR()`](https://jmw86069.github.io/splicejam/reference/annotateGRfromGR.md) -
    annotate `GRanges` from another.
  - [`annotateGRLfromGRL()`](https://jmw86069.github.io/splicejam/reference/annotateGRLfromGRL.md) -
    annotate `GRangesList` from another, keeping list elements distinct.
    E.g. exons within a transcript, exons within a gene. Good when
    operating on a sub-grouping, then keeping the original parent
    annotations.

## 14apr2026

- DONE. Cleaner workflow:

  - `assembleSashimiData()`? - does the work of
    [`sashimiDataConstants()`](https://jmw86069.github.io/splicejam/reference/sashimiDataConstants.md),
    returns an environment to be used in other functions.

  - `plotSashimiGene()`

    - does work of
      [`prepareSashimi()`](https://jmw86069.github.io/splicejam/reference/prepareSashimi.md),
      [`plotSashimi()`](https://jmw86069.github.io/splicejam/reference/plotSashimi.md)
    - adds gene panel using
      [`gene2gg()`](https://jmw86069.github.io/splicejam/reference/gene2gg.md)
    - uses patchwork as needed
    - zoom by coordinate or by exon name, as with the Shiny app.

- DEFER. Consider moving BGA plot functions to another smaller R
  package, accept other PCA objects. Consider tSNE/UMAP.

  - Dedicated R package for PCA/BGA-like features?
  - Or add to something like: jamma, jamses, or platjam
  - Consider dynamic group centroid calculations, with optional
    supergroups as with BGA.

- DEFER. Support non-gene/transcript features/tracks within the “gene
  panel”. The gene panel could include suitable “tracks”.

- DEFER. Ideal world: Insert Splicejam into something like plotgardner.
  Needs research to determine feasibility.

## 07apr2026

- Consider a way to specify gene/transcript order in
  [`gene2gg()`](https://jmw86069.github.io/splicejam/reference/gene2gg.md).

  - One example is with neighboring genes: MT1L, MT1E, MT1A. It “works”
    by showing all genes and transcripts, but all genes are shown
    together, then all transcripts are shown together, they are not
    ordered to keep transcripts with gene.
  - Transcripts are ordered alphanumerically, it might be useful to
    specify a particular order. (Can limits be used here?)

## 31mar2026

- [`launchSashimiApp()`](https://jmw86069.github.io/splicejam/reference/launchSashimiApp.md)
  Consider methods for user authentication.
- Revisit app layout, shiny UI components, packages, etc.

## 19mar2026

- **Make a proper vignette describing the key steps!**

- Add `testthat` unit tests.

- [`launchSashimiApp()`](https://jmw86069.github.io/splicejam/reference/launchSashimiApp.md)

  - DONE. Enable multi-column gene-exon plot.

  - Add options:

    - DONE. Hide junction label
    - DEFER. Hide exon label. (See gene: TTN, Ttn). ggrepel already
      hides labels when there are too many, but junction counts
      sometimes aren’t useful.
    - DEFER. Option to zoom y-axis range? Defer, since interactive plot
      should provide similar relief.
    - DONE. Focus y-axis range based on display coordinates.

  - PARTIAL. Enable more custom defaults: panel_height, share_y_axis,
    junction_arc_factor, default_gene (already works, use as a model).

  - Consider option to flatten exons by gene as needed, when not
    provided up-front. It only requires ‘txdb’.

  - Add “bookmarks” to store per-gene settings, like exon range.

  - Consider option to have **no gene** load at startup.

- DEFER. Debug slow processing on Windows, and
  `simpleError in seqinfo(con): UCSC library operation failed` even
  though `seqinfo()` isn’t called directly.

  - The only workaround is to have bigwig files local on Windows.
  - Apparently known issue on Windows, cannot load remote bigwig files.
    `import.bw()` calls `seqinfo(bw)` which calls `expandPath()` then
    `expandURL()` for remote files, and that process calls `GET()` with
    `nobody=1` which retrieves the entire file then discards the content
    body, which for bigwig files defeats the purpose of using a
    position-indexed file. That process is used simply to verify the
    URI. `GET(uri, config)$url` pulls the entire bigwig to get the url.
    Perhaps linux hosts correctly use only HEAD? Windows gets
    everything. Workaround is to use local files, though that’s just
    silly.
  - Speed is reasonably good with local bigwig files.

- Check memoise cache steps, when no coverage, it tries to repair. What
  if there is no coverage, it spends time trying to repair coverage.

- Simplify the workflow.

  - Consider SplicejamData object to contain the elements needed, making
    it clear how to prepare data, and confirm it is ready to use.
