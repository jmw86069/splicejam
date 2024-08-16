
# @import shiny
# @import shinydashboard
# @import htmltools

#' Sashimi Shiny app constants
#'
#' Sashimi Shiny app constants
#'
#' This function defines several constant values
#' used by the R-shiny Splicejam Sashimi viewer.
#' The required coverage and junction data is prepared
#' and defined by `sashimiDataConstants()`. The remaining
#' items used in the R-shiny app are defined for inline documentation
#' in the R-shiny app, including `aboutExtra` which is
#' included in the `"About"` tab, intended to describe
#' the source of data included in the R-shiny app.
#' Data is returned in an `environment` which by default
#' is the global environment `globalenv()`. However it
#' is recommended to use a custom environment, for
#' example: `shiny_envir <- new.env()`.
#'
#' When the R-shiny app is defined `launchSashimiApp()`,
#' it calls `shiny::shinyApp()` using arguments `server`, `ui`,
#' and `options`. This function `sashimiAppConstants()`
#' prepares `environment` which are assigned to the `ui`
#' and `server` objects. The process therefore makes data
#' inside these environments available to the `ui` and `server`
#' functions.
#'
#' The following values will be used from the environment,
#' searching up the environment parent chain until it finds
#' a match, until searching the global environment. Similarly,
#' this function also defines variables in the environment using
#' the `<<-` operator, which by default also searches up the
#' environment chain until it finds a match, otherwise
#' populating the global environment.
#'
#' If a variable is not found, the corresponding data will be
#' derived from relevant source data. If no data is provided,
#' the default argument `empty_uses_farrisdata=TRUE` means the
#' data and `filesDF` will use data from the publication Farris et al,
#' from Github package `"jmw86069/farrisdata"`.
#'
#' #' The `filesDF` object should be a `data.frame` with at least three colnames:
#'
#' * `"sample_id"`
#' * `"type"` (with values either `"bw"` or `"junction"`)
#' * `"url"` (a URL or file path to each file.)
#'
#' It can optionally include colname `"scale_factor"` with `numeric`
#' values used to multiply the coverage or junction values, the default
#' `scale_factor=1`.
#'
#' Other data derived by this function or by `sashimiDataConstants()`:
#'
#' * **color_sub**: `character` vector of R colors, whose names are
#'    used to match `filesDF$sample_id`. When not supplied,
#'    colors are defined by `colorjam::group2colors()` and
#'    `unique(filesDF$sample_id)`.
#' * **txdb**: `TranscriptDb` object used to derive `exonsByTx`
#'    and `cdsByTx` if either object does not already exist. If `txdb`
#'    is not supplied, it is derived from `gtf` using
#'    `GenomicFeatures::makeTxDbFromGFF()`.
#' * **tx2geneDF**: `data.frame` with colnames: `"transcript_id"` and
#'    `"gene_name"`.
#' * **gtf**: `character` path to a GTF/GFF/GFF3 file, suitable for
#'    `GenomicFeatures::makeTxDbFromGFF()`. The `gtf` is only used
#'    if `tx2geneDF` or `exonsByTx` are not supplied. Note that
#'    when `gtf` points to a remote server, the file is copied to
#'    the current working directory for more rapid use.
#'    If the file already exists in the local directory, it is re-used.
#' * **exonsByTx**: `GRangesList` object, named by `"transcript_id"`,
#'    containing all exons for each transcript. It is derived from `txdb`
#'    if not supplied; and names should match `tx2geneDF$transcript_id`.
#' * **cdsByTx**: `GRangesList` object, named by `"transcript_id"`,
#'    containing only CDS (protein-coding) exons for each transcript.
#'    It is derived from `txdb` if not supplied;
#'    and names should match `tx2geneDF$transcript_id`.
#' * **detectedTx**: `character` vector of `tx2geneDF$transcript_id` values,
#'    representing a subset of transcripts detected above background.
#'    See `definedDetectedTx()` for one strategy to define detected transcripts.
#'    If `detectedTx` does not exist, it is defined by all transcripts
#'    present in `tx2geneDF$transcript_id`. Note this step can be the
#'    rate-limiting step in the preparation of `flatExonsByTx`.
#' * **detectedGenes**: `character` vector of values that match
#'    `tx2geneDF$gene_name`. If it is not supplied, it is inferred
#'    from `detectedTx` and `tx2geneDF$transcript_id`.
#' * **flatExonsByGene**: `GRangesList` object containing non-overlapping
#'    exons for each gene, whose names match `tx2geneDF$gene_name`. If not
#'    supplied, it is derived using `flattenExonsBy()` and objects
#'    `exonsByTx`, `cdsByTx`, `detectedTx`, and `tx2geneDF`. This step is
#'    the key step for using a subset of detected transcripts, in order
#'    to produce a clean gene-exon model.
#' * **flatExonsByTx**: `GRangesList` object containing non-overlapping
#'    exons for each transcript. If not
#'    supplied, it is derived using `flattenExonsBy()` and objects
#'    `exonsByTx`, `cdsByTx`, `detectedTx`, and `tx2geneDF`. This step is
#'    the key step for using a subset of detected transcripts, in order
#'    to produce a clean transcript-exon model.
#'
#' When `use_memoise=TRUE` several R objects are cached using
#' `memoise::memoise()`, to help re-use of prepared R objects,
#' and to help speed the re-use of data within the R-shiny app:
#'
#' * **flatExonsByGene**
#' * **flatExonsByTx**
#' * **exonsByTx**
#' * **cdsByTx**
#'
#' To include a description of data used in your R-shiny app,
#' define the variable `aboutExtra` either using `character` text,
#' or as `htmltools::tags()` sufficient to be displayed in the
#' R-shiny UI. The content is displayed in the tab
#' `"About Sashimi Plots"` at the top of the app.
#'
#' @family splicejam R-shiny functions
#'
#' @return `environment` that contains the data required for the
#'    splicejam R-shiny app. It also includes data returned by
#'    `sashimiDataConstants()`. Note that if `envir` is supplied,
#'    the data will be updated inside that environment.
#'
#'
#' @param ... additional arguments are passed to `sashimiDataConstants()`
#' @param filesDF `data.frame` that contains at least these colnames:
#'    `"sample_id"`, `"url"`, `"type"`. This `data.frame` defines the
#'    source data used to create sashimi plots.
#' @param color_sub `character` vector of R colors, whose names
#'    match values in `filesDF$sample_id`. If not supplied,
#'    or if not all names are present in `color_sub`, the
#'    remaining names are converted to colors using
#'    `colorjam::group2colors()`.
#' @param aboutExtra character string or html tag from `"htmltools"`
#'    suitable for use in a R-shiny app. This text is displayed
#'    in the Help tab, and is intended to describe the data
#'    content shown in the R-shiny app.
#' @param envir `environment` in which the data should be loaded,
#'    which takes priority over argument `assign_global`.
#'    When `envir=NULL` and `assign_global=TRUE` the default
#'    environment is `globalenv()`. When `assign_global=FALSE`
#'    and `envir=NULL` a new environment is created using
#'    `new.env(parent=emptyenv())` so there is no parent environment,
#'    thereby preventing it from searching `globalenv()` for
#'    variables not defined in its own environment.
#' @param assign_global `logical` indicating whether the default
#'    environment should be `globalenv()`. Note this is not
#'    typically recommended, however it can be convenient to
#'    operate using only the user global environment, and is
#'    the default approach.
#' @param use_memoise `logical` indicating whether to use `memoise`
#'    to cache intermediate data files for exons, flattened exons,
#'    transcript-gene data, and so on. This mechanism reduces
#'    time to render sashimi plots that re-use the same gene.
#'    All memoise cache folders are named with `"_memoise"`.
#' @param empty_uses_farrisdata `logical` indicating whether to
#'    use data from the Github R package `"jmw86069/farrisdata"`
#'    if no data is supplied to this function. This behavior is
#'    intended to make it easy to use farrisdata to recreate
#'    the Sashimi plots in that publication.
#' @param gtf,txdb,tx2geneDF,exonsByTx,cdsByTx arguments passed to
#'    `sashimiDataConstants()`.
#' @param detectedTx,detectedGenes,flatExonsByGene,flatExonsByTx
#'    arguments passed to `sashimiDataConstants()`.
#' @param verbose `logical` indicating whether to print verbose output.
#'
#' @export
sashimiAppConstants <- function
(...,
 filesDF=NULL,
 color_sub=NULL,
 aboutExtra=NULL,
 envir=NULL,
 assign_global=TRUE,
 use_memoise=TRUE,
 empty_uses_farrisdata=TRUE,
 gtf=NULL,
 txdb=NULL,
 tx2geneDF=NULL,
 exonsByTx=NULL,
 cdsByTx=NULL,
 detectedTx=NULL,
 detectedGenes=NULL,
 flatExonsByGene=NULL,
 flatExonsByTx=NULL,
 verbose=FALSE)
{
   # define environment in which to store resulting data
   if (is.environment(envir)) {
      env_name <- deparse(substitute(envir));
      if (nchar(env_name) == 0) {
         env_name <- "envir"
      }
      if (verbose) {
         jamba::printDebug("sashimiAppConstants(): ",
            "Using ", "envir", " '", env_name, "' as provided.");
      }
   } else if (assign_global) {
      if (verbose) {
         jamba::printDebug("sashimiAppConstants(): ",
            "Using '", "globalenv()", "' due to assign_global=TRUE.");
      }
      envir <- globalenv();
      env_name <- "globalenv";
   } else {
      envir <- new.env(parent=emptyenv());
      env_name <- "new.env";
      if (verbose) {
         jamba::printDebug("sashimiAppConstants(): ",
            "Using '", "new.env()", "' due to assign_global=FALSE.");
      }
   }
   if (verbose) {
      jamba::printDebug("sashimiAppConstants(): ",
         "Using environment ", env_name);
      print(ls(envir=envir));
   }
   ## Quietly load an otherwise loud package dependency
   # 15jul2021: commented out for testing
   #suppressPackageStartupMessages(require(GenomicFeatures));

   envir$aboutExtra <- get_fn_envir("aboutExtra",
      envir=envir);
   # if present, include aboutExtra
   if (length(envir$aboutExtra) > 0) {
      if (!inherits(envir$aboutExtra, c("shiny.tag", "shiny.tag.list"))) {
         envir$aboutExtra <- htmltools::tags$p(envir$aboutExtra);
      }
      if (verbose) {
         jamba::printDebug("sashimiAppConstants(): ",
            "Using provided ", "aboutExtra");
      }
   }

   # update gene-transcript data as required
   envir <- sashimiDataConstants(gtf=gtf,
      txdb=txdb,
      tx2geneDF=tx2geneDF,
      exonsByTx=exonsByTx,
      cdsByTx=cdsByTx,
      detectedTx=detectedTx,
      detectedGenes=detectedGenes,
      flatExonsByGene=flatExonsByGene,
      flatExonsByTx=flatExonsByTx,
      empty_uses_farrisdata=empty_uses_farrisdata,
      use_memoise=use_memoise,
      verbose=verbose,
      envir=envir,
      ...);
   params <- c("filesDF",
      "color_sub");
   for (i in params) {
      assign(i,
         value=get_fn_envir(i,
            envir=envir,
            verbose=verbose - 1),
         envir=envir);
   }

   # assert filesDF is available with proper colnames
   if (length(envir$filesDF) == 0 ||
         nrow(envir$filesDF) == 0) {
      if (envir$empty_uses_farrisdata && nchar(system.file(package="farrisdata")) > 0) {
         if (verbose) {
            jamba::printDebug("sashimiAppConstants(): ",
               "Using filesDF from ",
               "farrisdata::farris_sashimi_files_df");
         }
         envir$filesDF <- farrisdata::farris_sashimi_files_df;
      } else {
         envir$empty_uses_farrisdata <- FALSE;
         stop("filesDF must be provided, or set empty_uses_farrisdata=TRUE to use Farris et al. data from Github R package 'jmw86069/farrisdata'");
      }
   }
   if (!all(c("sample_id", "url", "type") %in% colnames(envir$filesDF))) {
      stop("filesDF must contain colnames: 'sample_id', 'url', 'type'");
   }

   if (length(envir$aboutExtra) == 0 &&
         jamba::igrepHas("farris", envir$filesDF$url) &&
         envir$empty_uses_farrisdata) {
      if (length(envir$aboutExtra) == 0) {
         if (verbose) {
            jamba::printDebug("sashimiAppConstants(): ",
               "Using farrisdata ",
               "'aboutExtra'",
               " text.");
         }
         envir$aboutExtra <- htmltools::tags$p("Data is provided by the ",
            htmltools::strong("farrisdata"),
            " package, which provides mouse hippocampal subregion-
            and compartment-specific RNA-seq data described in
            Farris et al 2019. Each 'sample_id' represents the
            normalized RNA-seq aligned sequence coverage after
            combining three biological replicates per sample group.
            The purpose of this resource is to provide the community
            with a user-friendly interface to mine the data for
            isoform-specific differences across hippocampal subregions
            (CA1, CA2, CA3, DG) and subcellular compartments
            (CB = Cell Body, DE = Dendrites).");
      }
   }

   ## Define color_sub
   if (length(envir$color_sub) == 0 &&
         envir$empty_uses_farrisdata &&
         nchar(system.file(package="farrisdata")) > 0) {
      envir$color_sub <- farrisdata::colorSub;
      if (!all(envir$filesDF$sample_id %in% names(envir$color_sub))) {
         color_sub_new <- colorjam::group2colors(unique(envir$filesDF$sample_id));
         is_new <- setdiff(names(color_sub_new), names(envir$color_sub))
         envir$color_sub[is_new] <- color_sub_new[is_new];
      }
   }
   if (length(envir$color_sub) == 0 && nrow(envir$filesDF) > 0) {
      envir$color_sub <- colorjam::group2colors(unique(envir$filesDF$sample_id));
      if (verbose) {
         jamba::printDebug("sashimiAppConstants(): ",
            "Defined new color_sub for unique sample_id:\n",
            names(color_sub),
            fgText=list("darkorange", "dodgerblue",
               jamba::setTextContrastColor(color_sub)),
            bgText=list(NA, NA, color_sub));
      }
   }


   # guides
   # define guides tab
   envir$guidesTab <- shiny::fluidPage(
      htmltools::tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
      shiny::sidebarLayout(
         shiny::mainPanel(
            width=7,
            shinydashboard::tabBox(
               width=12,
               shiny::tabPanel(
                  title="About Sashimi Plots",
                  shiny::uiOutput("sashimiplot_guide")
               ),
               shiny::tabPanel(
                  title="Creating a Sashimi Plot",
                  shiny::uiOutput("sashimiplotviz_guide")
               )
            )
         ),
         shiny::sidebarPanel(
            width=5,
            "Sashimi viewer visualized sequence coverage
            data alongside splice junction-spanning sequence reads,
            using compressed intron genomic coordinates.",
            htmltools::tags$ul(
               htmltools::tags$li(
                  htmltools::strong(style="color:firebrick",
                     "The methods were developed in support of this publication"),
                  htmltools::br(),
                  htmltools::a("S. Farris, J. M. Ward, K.E. Carstens, M. Samadi, Y. Wang and S. M. Dudek. ",
                     "Cell Reports 2019. ",
                     htmltools::em("Hippocampal subregions express distinct dendritic transcriptomes that reveal unexpected differences in mitochondrial function in CA2."),
                     href="https://github.com/jmw86069/jampack")
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:firebrick",
                     "Sashimi plots were originally envisioned by MISO:"),
                  htmltools::br(),
                  htmltools::a("Katz, Y, Wang ET, Silterra J, Schwartz S, Wong B, Thorvaldsdóttir H, Robinson JT, Mesirov JP, Airoldi EM, Burge, CB.:",
                     htmltools::em("Sashimi plots: Quantitative visualization of alternative isoform expression from RNA-seq data."),
                     href="http://biorxiv.org/content/early/2014/02/11/002576")
               )
            ),
            htmltools::tags$p("Relevant R version info:"),
            htmltools::tags$ul(
               htmltools::tags$li(
                  htmltools::strong(style="color:black", R.version.string)
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "jampack:"),
                  as.character(packageVersion("jampack"))
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "splicejam:"),
                  as.character(packageVersion("splicejam"))
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "jamba:"),
                  as.character(packageVersion("jamba"))
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "colorjam:"),
                  as.character(packageVersion("colorjam"))
               ),
               if (jamba::check_pkg_installed("farrisdata")) {
                  htmltools::tags$li(
                     htmltools::strong(style="color:black", "farrisdata:"),
                     as.character(packageVersion("farrisdata"))
                  )
               } else {
                  htmltools::tags$li(
                     htmltools::strong(style="color:black", "farrisdata:"),
                     as.character("not installed")
                  )
               },
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "ggplot2:"),
                  as.character(packageVersion("ggplot2"))
               ),
               htmltools::tags$li(
                  htmltools::strong(style="color:black", "plotly:"),
                  as.character(packageVersion("plotly"))
               )
            ),
            shiny::fluidRow(
               shiny::column(
                  width=12,
                  style="padding:0px",
                  shinydashboardPlus::box(
                     title="Full R sessionInfo():",
                     #status="warning",
                     solidHeader=TRUE,
                     collapsible=TRUE,
                     collapsed=TRUE,
                     width=12,
                     shiny::pre(shiny::htmlOutput("sessionInfo"))
                  )
               )
            )
         )
      )
   );

   # sashimiplot_guide
   envir$sashimiplot_guide <- shiny::fluidPage(
      htmltools::h1("About Sashimi Plots",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         aboutExtra,
         htmltools::tags$h3("The Sashimi Plot Tab"),
         htmltools::tags$p(
            "The Sashimi Plot tab",
            " provides a visualization of
            gene transcript data from one or more biological samples,
            specific for RNA-seq (RNA sequencing of gene transcripts.)
            It combines several types of data important for insightful
            interpretation of the raw results:"),
         htmltools::tags$ul(
            htmltools::tags$li(
               htmltools::strong("Sequence coverage", style="color:navy"),
               " - a visual indication of the number of sequencing
               reads that overlap each position along the genome.
               Coverage data is presented as filled polygon, where
               the y-axis score is proportional to the number of
               normalized RNA-seq reads aligned."
            ),
            htmltools::tags$li(
               htmltools::strong("Splice junction reads", style="color:navy"),
               " - sequence reads that have a special gapped alignment
               across two separate genome locations, often defined during
               alignment to transcript (exon) sequences. Junctions are
               represented as wide arcs, where the y-axis height of the
               arc is proportional to the number of normalized RNA-seq
               reads support the gapped alignment across this junction."
            ),
            htmltools::tags$li(
               htmltools::strong("Gene-transcript-exon model", style="color:navy"),
               htmltools::tags$ul(
                  htmltools::tags$li(
                     "The gene-exon model is defined for each gene locus, and is
                     visualized as a series of rectangles joined by thin lines,
                     where each rectangle is an 'exon' and each thin line is
                     an 'intron' between exons. Note that all transcript
                     models are flattened to produce the gene-exon model.
                     Exons are numbered using contiguous exons, and each
                     sub-section of an exon is given a letter suffix, for example:
                     exon1a, exon1b, exon2, exon3, exon4a, exon4b."
                  ),
                  htmltools::tags$li(
                     "The transcript-exon models are defined for each
                     transcript isoform annotated to the gene locus. A
                     subset are annotated 'detected' in order to hide
                     transcript isoforms for which there is little or no
                     supporting RNA-seq data evidence."
                  )
               )
            ),
            htmltools::tags$li(
               htmltools::strong("Compressed intron coordinates", style="color:navy"),
               " - the intron regions are visually compressed along the x-axis,
               because the purpose of this Sashimi plot is to visualize
               gene transcript data relative to transcript-scale features,
               and not relative to genome-scale features. Introns are
               typically 10x to 100x larger than exons, and their magnitude
               would otherwise obscure all exon features."
            )
         )
      )
   );

   envir$sashimiplotviz_guide <- shiny::fluidPage(
      htmltools::h1("Creating a Sashimi Plot",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         htmltools::tags$p("The typical workflow for viewing a Sashimi plot is described
            below:"),
         htmltools::tags$h3("Select the Sashimi Plot tab"),
         htmltools::tags$ul(
            htmltools::tags$li(
               htmltools::strong("Select a gene", style="color:navy"),
               "- all detected genes are searchable",
               htmltools::tags$ul(
                  htmltools::tags$li(
                     htmltools::strong("Slider bar measurement"),
                     " - optionally define a smaller region for the x-axis range,
                     where 'coordinates' uses genome coordianates, and
                     'exon names' uses fixed gene-exon boundaries."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Show coverage by strand"),
                     " - optionally restrict the display of RNA-seq coverage
                     to a specific DNA strand, (+) refers to the top strand
                     for genes transcribed left-to-right, and (-) refers to
                     the bottom strand for genes transcribed right-to-left."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Minimum junction reads"),
                     " - optionally define a minimum threshold for junction
                     reads. This option is helpful to reduce the display of
                     spurious junctions that have very few supporting RNA-seq
                     reads."
                  )
               )
            ),
            htmltools::tags$li(
               htmltools::strong("Click 'Update Sashimi Plots'", style="color:navy"),
               htmltools::tags$ul(
                  htmltools::tags$li(
                     " - when this button is enabled, clicking will begin
                     downloading and processing a new Sashimi plot. Data is
                     cached as it is produced, so the R-shiny app should
                     become more responsive over time."
                  )
               )
            ),
            htmltools::tags$li(
               htmltools::strong("Click the icon '",
                  style="color:navy"),
               shiny::icon("info"),
               htmltools::strong("' for more visual options", style="color:navy"),
               htmltools::tags$ul(
                  htmltools::tags$li(
                     htmltools::strong("Height per panel"),
                     " - define the pixel height of each Sashimi plot panel,
                     with one panel per biological sample."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Font sizing"),
                     " - optionally scale up or down the overall font size"
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Interactive plot"),
                     " - when checked, render using plotly with interactive
                     features. This feature is under active development to
                     enable as many useful features as possible. Uncheck
                     to view a static plot as created using ggplot."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Show gene-exon model"),
                     " - when checked, the flattened gene-exon model is
                     displayed below the Sashimi panels. When using
                     interactive plotting, and one column of panels, the
                     x-axis range can be zoomed by clicking and dragging."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Show transcript-exon model"),
                     " - when checked, the transcript-exon model is
                     displayed below the Sashimi panels, including all
                     transcripts, or a subset of 'detected' transcripts.
                     Viewing the transcript models can be helpful when
                     interpreting which transcript isoform may be
                     differentially regulated across biological samples."
                  ),
                  htmltools::tags$li(
                     htmltools::strong("Shared y-axis range"),
                     " - when checked, all Sashimi panels share the same
                     y-axis range, which helps visualize differences in
                     absolute gene expression levels. When unchecked, each
                     panel is independently scaled, which helps interpret
                     changes in transcript isoforms across samples."
                  )
               )
            )
         )
      )
   );

   envir$dataresources_guide <- shiny::fluidPage(
      htmltools::h1("Data Resources",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         htmltools::tags$ul(
            htmltools::tags$li(
               htmltools::strong("List of data resources",
                  style="color:dimgrey"),
               " go here.")
         )
      )
   );
   envir$nbsp <- htmltools::HTML("&nbsp;");
   envir$nbsp3 <- htmltools::HTML("&nbsp;&nbsp;&nbsp;");

   return(invisible(envir))
}
