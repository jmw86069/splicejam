
#' Sashimi Shiny app constants
#'
#' Sashimi Shiny app constants
#'
#' This function defines several constant values
#' used by the R-shiny Splicejam Sashimi viewer,
#' typically taking from the global environment
#' where possible.
#'
#' The R-shiny app is started by `launchSashimiApp()`, which
#' calls `shiny::shinyApp()`, using arguments `server`, `ui`,
#' `onStart`, and `options`. This function fulfills the
#' argument `onStart`.
#'
#' This function is intended to be called only from within
#' an R-shiny app, since it defines several variables in the
#' parent (global) environment.
#'
#' The following values will be used from the environment,
#' searching up the environment parent chain until it finds
#' a match, until searching the global environment. Similarly,
#' this function also defines variables in the environment using
#' the `<<-` operator, which by default also searches up the
#' environment chain until it finds a match, otherwise
#' populating the global environment.
#'
#' If a variable is not found, the corresponding default
#' from the `farrisdata` package will be used.
#'
#' * **filesDF**: `data.frame` with colnames: `"sample_id"`,
#' `"type"`, `"url"`, `"scale_factor"`. When `"filesDF"` does
#' not exist, it is taken from `farrisdata::farris_sashimi_files_df`,
#' and the default `"aboutExtra"` text is used.
#' * **color_sub**: character vector of R colors, whose names are
#' used to match `"sample_id"` or other values used to colorize
#' ggplot visualizations. It is either populated from
#' `farrisdata::colorSub`, or `colorjam::group2colors()` is called
#' using the unique `sample_id` values in `filesDF`.
#' * **txdb**: TxDb object used to derive `exonsByTx`, and `cdsByTx`
#' if either does not exist. If the TxDb object does not exist,
#' it is created from the `gtf` file using
#' `GenomicFeatures::makeTxDbFromGFF()`.
#' * **tx2geneDF**: `data.frame` with colnames: `"transcript_id"`,
#' `"gene_name"`. This data is used to convert `"gene_name"` to
#' a vector of `"transcript_id"` values.
#' * **gtf**: character path to a GTF/GFF/GFF3 file, suitable for
#' `GenomicFeatures::makeTxDbFromGFF()`. The `gtf` is only used
#' if `tx2geneDF` or `exonsByTx` are not available. Note that
#' when `gtf` points to a remote server, the file is copied to
#' the current working directory. If the file already exists in
#' the local directory, it is re-used.
#' * **exonsByTx**: `GRangesList` object, named by `"transcript_id"`,
#' containing all exons for each transcript.
#' * **cdsByTx**: `GRangesList` object, named by `"transcript_id"`,
#' containing only CDS (protein-coding) exons for each transcript.
#' * **detectedTx**: character vector of `"transcript_id"` values,
#' representing the subset of transcripts detected above background.
#' If it does not exist, it is derived from `farrisdata::farrisTxSE`
#' `"TxDetectedByTPM"`. If those values are not all present in
#' `tx2geneDF` then `detectedTx` is defined by all transcripts
#' present in `tx2geneDF$transcript_id`.
#' * **detectedGenes**: inferred using tx2geneDF and detectedTx.
#' * **flatExonsByGene**: `GRangesList` object containing non-overlapping
#' exons for each gene, or it is derived from `exonsByTx`,
#' `cdsByTx`, `detectedTx`, and `tx2geneDF`, using
#' `flattenExonsBy()`.
#' * **flatExonsByTx**: `GRangesList` object containing non-overlapping
#' exons for each gene, or it is derived from `exonsByTx`,
#' `cdsByTx`, `detectedTx`, and `tx2geneDF`, using
#' `flattenExonsBy()`.
#'
#' Several R objects are cached using `memoise::memoise()`, to
#' avoid having to create the R object each time the R-shiny app
#' is started:
#'
#' * **flatExonsByGene**
#' * **flatExonsByTx**
#' * **exonsByTx**
#' * **cdsByTx**
#'
#' A special custom variable is used to describe the specific
#' data included with the R-shiny app, `"aboutExtra"`. The
#' `"aboutExtra"` variable is expected to contain HTML tags,
#' for example from `htmltools::tags()` to format the text.
#' The `"aboutExtra"` content is displayed on the tab
#' `"About Sashimi Plots"` at the top.
#'
#' @family splicejam R-shiny functions
#'
#' @import shiny
#' @import shinydashboard
#' @import htmltools
#'
#' @param ... additional arguments are ignored.
#' @param empty_uses_farrisdata logical indicating whether to
#'    use data from the `"farrisdata"` R package if no default
#'    values are provided, and if that package is available.
#' @param aboutExtra character string or html tag from `"htmltools"`
#'    suitable for use in a shiny app. This text is intended
#'    to display a description of data visualized in the R-shiny
#'    app, for example when run with `launchSashimiApp()`.
#'
#' @export
sashimiAppConstants <- function
(...,
 assign_global=TRUE,
 empty_uses_farrisdata=TRUE,
 use_memoise=TRUE,
 aboutExtra=NULL,
 filesDF=NULL,
 color_sub=NULL,
 gtf=NULL,
 txdb=NULL,
 tx2geneDF=NULL,
 exonsByTx=NULL,
 cdsByTx=NULL,
 detectedTx=NULL,
 detectedGenes=NULL,
 flatExonsByGene=NULL,
 flatExonsByTx=NULL)
{
   ## Define filesDF here
   # dots <- list(...);
   if (assign_global) {
      sashimi_env <- globalenv();
   } else {
      sashimi_env <- new.env();
   }
   ## Quietly load an otherwise loud package dependency
   suppressPackageStartupMessages(require(GenomicFeatures));

   if (length(aboutExtra) == 0) {
      if (exists("aboutExtra", envir=globalenv(), inherits=FALSE)) {
         aboutExtra <- get("aboutExtra", globalenv());
      }
   }
   if (length(aboutExtra) > 0) {
      if (!inherits(aboutExtra, c("shiny.tag", "shiny.tag.list"))) {
         aboutExtra <- htmltools::tags$p(aboutExtra);
      }
      printDebug("Using existing ",
         "'aboutExtra'",
         " as provided.")
      assign("aboutExtra",
         value=aboutExtra,
         envir=sashimi_env);
   }
   params <- c("filesDF",
      "color_sub",
      "gtf",
      "txdb",
      "tx2geneDF",
      "exonsByTx",
      "cdsByTx",
      "detectedTx",
      "detectedGenes",
      "flatExonsByGene",
      "flatExonsByTx");
   for (i in params) {
      if (length(get(i)) == 0 && exists(i, envir=globalenv(), inherits=FALSE)) {
         ival <- get(i,
            envir=globalenv(),
            inherits=FALSE);
         ## Assign to local function space
         assign(i,
            ival);
         ## Assign to the specified environment
         assign(i,
            value=ival,
            envir=sashimi_env);
      }
   }
   if (length(filesDF) == 0 ||
         nrow(filesDF) == 0) {
      if (empty_uses_farrisdata && suppressPackageStartupMessages(require(farrisdata))) {
         printDebug("Using filesDF from ",
            "farrisdata::farris_sashimi_files_df");
         data(farris_sashimi_files_df);
         filesDF <- farris_sashimi_files_df;
         assign("filesDF",
            value=filesDF,
            envir=sashimi_env);
      } else {
         empty_uses_farrisdata <- FALSE;
      }
   }
   if (length(aboutExtra) == 0 &&
         jamba::igrepHas("farris", filesDF$url) &&
         empty_uses_farrisdata) {
      if (length(aboutExtra) == 0) {
         printDebug("Using farrisdata ",
            "'aboutExtra'",
            " text.");
         aboutExtra <- htmltools::tags$p("Data is provided by the ",
            strong("farrisdata"),
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
         assign("aboutExtra",
            value=aboutExtra,
            envir=sashimi_env);
      }
   }
   if (!all(c("sample_id", "url", "type") %in% colnames(filesDF))) {
      stop("filesDF must contain colnames: 'sample_id', 'url', and 'type'.");
   }

   ## Define color_sub
   if (length(color_sub) == 0 &&
         empty_uses_farrisdata &&
         suppressPackageStartupMessages(require(farrisdata))) {
      color_sub <- farrisdata::colorSub;
      if (!all(filesDF$sample_id %in% names(color_sub))) {
         color_sub_new <- colorjam::group2colors(unique(filesDF$sample_id));
         is_new <- setdiff(names(color_sub_new), names(color_sub))
         color_sub[is_new] <- color_sub_new[is_new];
      }
   }
   if (length(color_sub) == 0 && nrow(filesDF) > 0) {
      color_sub <- colorjam::group2colors(filesDF$sample_id);
   }
   assign("color_sub",
      value=farrisdata::colorSub,
      envir=sashimi_env);

   ## Define flat exons by gene
   ## One-time setup cost when using GTF input
   if (length(tx2geneDF) == 0 || length(exonsByTx) == 0) {
      if (length(gtf) == 0 && length(txdb) == 0 && empty_uses_farrisdata) {
         # use default GTF file if not defined
         gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz";
         printDebug("Defining default gtf file:",
            gtf,
            " from which ",
            c("tx2geneDF, exonsByTx, and cdsByTx"),
            " will be derived.");
      }
      if (length(gtf) == 0 && length(txdb) == 0) {
         stop(paste0("The 'gtf' or 'txdb' argument are required ",
            "when either 'tx2geneDF' or 'exonsByTx' are not provided."));
      }
      if (length(gtf) > 0) {
         gtfBase <- basename(gtf);
         if (!file.exists(gtfBase)) {
            printDebug("Downloading gtf:", gtf,
               " to:", gtfBase);
            curl::curl_download(url=gtf,
               destfile=gtfBase);
         }
         if (length(tx2geneDF) == 0) {
            tx2geneFile <- gsub("[.](gff|gff3|gtf).*$",
               ".tx2geneDF.txt",
               gtfBase,
               ignore.case=TRUE);
            if (!file.exists(tx2geneFile)) {
               printDebug("Deriving tx2geneDF from gtf:",
                  gtfBase,
                  " then storing:",
                  tx2geneFile);
               tx2geneDF <- makeTx2geneFromGtf(GTF=gtfBase,
                  verbose=FALSE);
               write.table(file=tx2geneFile,
                  x=tx2geneDF,
                  sep="\t",
                  quote=FALSE,
                  na="",
                  col.names=TRUE);
            } else {
               printDebug("Reloading stored tx2geneDF:",
                  tx2geneFile);
               tx2geneDF <- read.table(file=tx2geneFile,
                  sep="\t",
                  header=TRUE,
                  quote="\"",
                  comment.char="");
            }
         } else {
            printDebug("Using tx2geneDF from environment.");
         }
      }
      if (length(tx2geneDF) == 0) {
         stop("The 'tx2geneDF' argument is required when 'gtf' is not supplied.");
      }
      # Assign in the parent environment
      assign("tx2geneDF",
         value=tx2geneDF,
         envir=sashimi_env);
      ## Now make TxDb in order to derive exonsByTx and cdsByTx
      #if (length(exonsByTx) == 0 || length(cdsByTx) == 0) {
      if (length(exonsByTx) == 0) {
         localDb <- gsub("[.](gff|gff3|gtf).*$",
            ".txdb",
            gtfBase,
            ignore.case=TRUE);
         if (!file.exists(localDb)) {
            printDebug("Deriving txdb from gtf:",
               gtfBase,
               " to store as:",
               localDb);
            txdb <- GenomicFeatures::makeTxDbFromGFF(gtfBase);
            AnnotationDbi::saveDb(x=txdb, file=localDb);
         } else {
            printDebug("Reloading txdb:", localDb);
            txdb <- AnnotationDbi::loadDb(file=localDb);
         }
         if (!DBI::dbIsValid(AnnotationDbi::dbconn(txdb))) {
            printDebug("Refreshing db connection:", localDb);
            txdb <- AnnotationDbi::loadDb(file=localDb);
         }
         # Assign in the parent environment
         assign("txdb",
            value=txdb,
            envir=sashimi_env);

         # First obtain exons by transcript
         printDebug("Deriving exonsByTx from txdb");
         suppressPackageStartupMessages(require(GenomicFeatures));
         if (use_memoise) {
            exonsBy_m <- memoise::memoise(exonsBy,
               cache=memoise::cache_filesystem("exonsBy_memoise"));
            exonsBy_m_cached <- memoise::has_cache(exonsBy_m)(
               txdb,
               by="tx",
               use.names=TRUE);
            printDebug("exonsBy_m_cached:",
               exonsBy_m_cached);
         } else {
            exonsBy_m <- exonsBy;
         }
         exonsByTx <- exonsBy_m(
            txdb,
            by="tx",
            use.names=TRUE);
         values(exonsByTx@unlistData)$feature_type <- "exon";
         values(exonsByTx@unlistData)$subclass <- "exon";
      } else {
         if (!"feature_type" %in% names(values(exonsByTx@unlistData))) {
            values(exonsByTx@unlistData)$feature_type <- "exon";
         }
         if (!"subclass" %in% names(values(exonsByTx@unlistData))) {
            values(exonsByTx@unlistData)$subclass <- "exon";
         }
      }
      # Assign in the parent environment
      assign("exonsByTx",
         value=exonsByTx,
         envir=sashimi_env);

      if (length(cdsByTx) == 0 && exists("txdb")) {
         printDebug("Deriving cdsByTx from txdb");
         suppressPackageStartupMessages(require(GenomicFeatures));
         if (use_memoise) {
            cdsBy_m <- memoise::memoise(cdsBy,
               cache=memoise::cache_filesystem("cdsBy_memoise"));
            cdsBy_m_cached <- memoise::has_cache(cdsBy_m)(
               txdb,
               by="tx",
               use.names=TRUE);
            printDebug("cdsBy_m_cached:",
               cdsBy_m_cached);
         } else {
            cdsBy_m <- cdsBy;
         }

         cdsByTx <- cdsBy_m(
            txdb,
            by="tx",
            use.names=TRUE);
         values(cdsByTx@unlistData)$feature_type <- "cds";
         values(cdsByTx@unlistData)$subclass <- "cds";
      } else {
         if (!"feature_type" %in% names(values(cdsByTx@unlistData))) {
            values(cdsByTx@unlistData)$feature_type <- "exon";
         }
         if (!"subclass" %in% names(values(cdsByTx@unlistData))) {
            values(cdsByTx@unlistData)$subclass <- "exon";
         }
      }
      # Assign in the parent environment
      assign("cdsByTx",
         value=cdsByTx,
         envir=sashimi_env);
   }

   ## Define detectedTx
   if (!exists("detectedTx") || length(detectedTx) == 0) {
      if (empty_uses_farrisdata && suppressPackageStartupMessages(require(farrisdata))) {
         data(farrisTxSE);
         printDebug("Using detectedTx from farrisTxSE");
         detectedTx <- subset(SummarizedExperiment::rowData(farrisTxSE),
            TxDetectedByTPM)$transcript_id;
      }
      if (!all(detectedTx %in% tx2geneDF$transcript_id)) {
         printDebug("Using detectedTx <- unique(tx2geneDF$transcript_id)");
         detectedTx <- unique(tx2geneDF$transcript_id);
      }
   } else {
      detectedTx <- unique(tx2geneDF$transcript_id);
   }
   # Assign in the parent environment
   assign("detectedTx",
      value=detectedTx,
      envir=sashimi_env);

   ## Infer available genes
   if (length(detectedGenes) == 0) {
      printDebug("Inferring detectedGenes from ",
         format(big.mark=",", length(detectedTx)),
         " detectedTx entries.");
      #printDebug("head(tx2geneDF, 3):");
      #print(head(tx2geneDF, 3));
      detectedGenes <- jamba::mixedSort(
         unique(
            subset(tx2geneDF,
               transcript_id %in% detectedTx)$gene_name));
   }
   # Assign in the parent environment
   assign("detectedGenes",
      value=detectedGenes,
      envir=sashimi_env);

   ## Define flatExonsByGene
   ## define memoised function
   if (use_memoise) {
      flattenExonsBy_m <- memoise::memoise(flattenExonsBy,
         cache=memoise::cache_filesystem("flattenExonsBy_memoise"));
   } else {
      flattenExonsBy_m <- flattenExonsBy;
   }

   if (length(flatExonsByGene) == 0) {
      printDebug("Deriving flatExonsByGene from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      if (use_memoise) {
         flattenExonsByGene_m_cached <- memoise::has_cache(flattenExonsBy_m)(
            exonsByTx=exonsByTx,
            cdsByTx=cdsByTx,
            detectedTx=detectedTx,
            by="gene",
            tx2geneDF=tx2geneDF,
            verbose=FALSE);
         printDebug("flattenExonsByGene_m_cached:",
            flattenExonsByGene_m_cached);
      }
      flatExonsByGene <- flattenExonsBy_m(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         detectedTx=detectedTx,
         by="gene",
         tx2geneDF=tx2geneDF,
         verbose=FALSE);
   }
   # Assign in the parent environment
   assign("flatExonsByGene",
      value=flatExonsByGene,
      envir=sashimi_env);

   ## Define flatExonsByTx
   if (length(flatExonsByTx) == 0) {
      printDebug("Deriving flatExonsByTx from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      if (use_memoise) {
         flattenExonsByTx_m_cached <- memoise::has_cache(flattenExonsBy_m)(
            exonsByTx=exonsByTx,
            cdsByTx=cdsByTx,
            detectedTx=detectedTx,
            tx2geneDF=tx2geneDF,
            by="tx",
            verbose=FALSE);
         printDebug("flattenExonsByTx_m_cached:",
            flattenExonsByTx_m_cached);
      }
      flatExonsByTx <- flattenExonsBy_m(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         detectedTx=detectedTx,
         tx2geneDF=tx2geneDF,
         by="tx",
         verbose=FALSE)
   }
   # Assign in the parent environment
   assign("flatExonsByTx",
      value=flatExonsByTx,
      envir=sashimi_env);

   if (!exists("verbose")) {
      verbose <- FALSE;
   }

   # guides
   # define guides tab
   guidesTab <- fluidPage(
      tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
      sidebarLayout(
         mainPanel(
            width=7,
            tabBox(
               width=12,
               tabPanel(
                  title="About Sashimi Plots",
                  uiOutput("sashimiplot_guide")
               ),
               tabPanel(
                  title="Creating a Sashimi Plot",
                  uiOutput("sashimiplotviz_guide")
               )
            )
         ),
         sidebarPanel(
            width=5,
            "Sashimi viewer visualized transcriptome RNA-seq coverage
            data alongside splice junction-spanning sequence reads,
            using compressed intron genomic coordinates.",
            tags$ul(
               tags$li(
                  strong(style="color:firebrick",
                     "The methods were developed in support of this manuscript"),
                  br(),
                  a("S. Farris, J. M. Ward, K.E. Carstens, M. Samadi, Y. Wang and S. M. Dudek. ",
                     "Cell Reports 2019 (Accepted). ",
                     em("Hippocampal subregions express distinct dendritic transcriptomes that reveal unexpected differences in mitochondrial function in CA2."),
                     href="https://github.com/jmw86069/jampack")
               ),
               tags$li(
                  strong(style="color:firebrick",
                     "Sashimi plots were originally envisioned by MISO:"),
                  br(),
                  a("Katz, Y, Wang ET, Silterra J, Schwartz S, Wong B, Thorvaldsdóttir H, Robinson JT, Mesirov JP, Airoldi EM, Burge, CB.:",
                     em("Sashimi plots: Quantitative visualization of alternative isoform expression from RNA-seq data."),
                     href="http://biorxiv.org/content/early/2014/02/11/002576")
               )
            ),
            tags$p("Relevant R version info:"),
            tags$ul(
               tags$li(
                  strong(style="color:black", R.version.string)
               ),
               tags$li(
                  strong(style="color:black", "jampack:"),
                  as.character(packageVersion("jampack"))
               ),
               tags$li(
                  strong(style="color:black", "splicejam:"),
                  as.character(packageVersion("splicejam"))
               ),
               tags$li(
                  strong(style="color:black", "jamba:"),
                  as.character(packageVersion("jamba"))
               ),
               tags$li(
                  strong(style="color:black", "colorjam:"),
                  as.character(packageVersion("colorjam"))
               ),
               if(suppressPackageStartupMessages(require(farrisdata))) {
                  tags$li(
                     strong(style="color:black", "farrisdata:"),
                     as.character(packageVersion("farrisdata"))
                  )
               } else {
                  tags$li(
                     strong(style="color:black", "farrisdata:"),
                     as.character("not installed")
                  )
               },
               tags$li(
                  strong(style="color:black", "ggplot2:"),
                  as.character(packageVersion("ggplot2"))
               ),
               tags$li(
                  strong(style="color:black", "plotly:"),
                  as.character(packageVersion("plotly"))
               )
            )
         )
      )
   );
   assign("guidesTab",
      value=guidesTab,
      envir=sashimi_env);

   # sashimiplot_guide
   sashimiplot_guide <- fluidPage(
      h1("About Sashimi Plots",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         aboutExtra,
         tags$h3("The Sashimi Plot Tab"),
         tags$p(
            "The Sashimi Plot tab",
            " provides a visualization of
            gene transcript data from one or more biological samples,
            specific for RNA-seq (RNA sequencing of gene transcripts.)
            It combines several types of data important for insightful
            interpretation of the raw results:"),
         tags$ul(
            tags$li(
               strong("Sequence coverage", style="color:navy"),
               " - a visual indication of the number of sequencing
               reads that overlap each position along the genome.
               Coverage data is presented as filled polygon, where
               the y-axis score is proportional to the number of
               normalized RNA-seq reads aligned."
            ),
            tags$li(
               strong("Splice junction reads", style="color:navy"),
               " - sequence reads that have a special gapped alignment
               across two separate genome locations, often defined during
               alignment to transcript (exon) sequences. Junctions are
               represented as wide arcs, where the y-axis height of the
               arc is proportional to the number of normalized RNA-seq
               reads support the gapped alignment across this junction."
            ),
            tags$li(
               strong("Gene-transcript-exon model", style="color:navy"),
               tags$ul(
                  tags$li(
                     "The gene-exon model is defined for each gene locus, and is
                     visualized as a series of rectangles joined by thin lines,
                     where each rectangle is an 'exon' and each thin line is
                     an 'intron' between exons. Note that all transcript
                     models are flattened to produce the gene-exon model.
                     Exons are numbered using contiguous exons, and each
                     sub-section of an exon is given a letter suffix, for example:
                     exon1a, exon1b, exon2, exon3, exon4a, exon4b."
                  ),
                  tags$li(
                     "The transcript-exon models are defined for each
                     transcript isoform annotated to the gene locus. A
                     subset are annotated 'detected' in order to hide
                     transcript isoforms for which there is little or no
                     supporting RNA-seq data evidence."
                  )
               )
            ),
            tags$li(
               strong("Compressed intron coordinates", style="color:navy"),
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
   assign("sashimiplot_guide",
      value=sashimiplot_guide,
      envir=sashimi_env);

   sashimiplotviz_guide <- fluidPage(
      h1("Creating a Sashimi Plot",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$p("The typical workflow for viewing a Sashimi plot is described
            below:"),
         tags$h3("Select the Sashimi Plot tab"),
         tags$ul(
            tags$li(
               strong("Select a gene", style="color:navy"),
               "- all detected genes are searchable",
               tags$ul(
                  tags$li(
                     strong("Slider bar measurement"),
                     " - optionally define a smaller region for the x-axis range,
                     where 'coordinates' uses genome coordianates, and
                     'exon names' uses fixed gene-exon boundaries."
                  ),
                  tags$li(
                     strong("Show coverage by strand"),
                     " - optionally restrict the display of RNA-seq coverage
                     to a specific DNA strand, (+) refers to the top strand
                     for genes transcribed left-to-right, and (-) refers to
                     the bottom strand for genes transcribed right-to-left."
                  ),
                  tags$li(
                     strong("Minimum junction reads"),
                     " - optionally define a minimum threshold for junction
                     reads. This option is helpful to reduce the display of
                     spurious junctions that have very few supporting RNA-seq
                     reads."
                  )
               )
            ),
            tags$li(
               strong("Click 'Update Sashimi Plots'", style="color:navy"),
               tags$ul(
                  tags$li(
                     " - when this button is enabled, clicking will begin
                     downloading and processing a new Sashimi plot. Data is
                     cached as it is produced, so the R-shiny app should
                     become more responsive over time."
                  )
               )
            ),
            tags$li(strong("Click the icon '", style="color:navy"),
               icon("info"),
               strong("' for more visual options", style="color:navy"),
               tags$ul(
                  tags$li(
                     strong("Height per panel"),
                     " - define the pixel height of each Sashimi plot panel,
                     with one panel per biological sample."
                  ),
                  tags$li(
                     strong("Font sizing"),
                     " - optionally scale up or down the overall font size"
                  ),
                  tags$li(
                     strong("Interactive plot"),
                     " - when checked, render using plotly with interactive
                     features. This feature is under active development to
                     enable as many useful features as possible. Uncheck
                     to view a static plot as created using ggplot."
                  ),
                  tags$li(
                     strong("Show gene-exon model"),
                     " - when checked, the flattened gene-exon model is
                     displayed below the Sashimi panels. When using
                     interactive plotting, and one column of panels, the
                     x-axis range can be zoomed by clicking and dragging."
                  ),
                  tags$li(
                     strong("Show transcript-exon model"),
                     " - when checked, the transcript-exon model is
                     displayed below the Sashimi panels, including all
                     transcripts, or a subset of 'detected' transcripts.
                     Viewing the transcript models can be helpful when
                     interpreting which transcript isoform may be
                     differentially regulated across biological samples."
                  ),
                  tags$li(
                     strong("Shared y-axis range"),
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
   assign("sashimiplotviz_guide",
      value=sashimiplotviz_guide,
      envir=sashimi_env);

   dataresources_guide <- fluidPage(
      h1("Data Resources",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(strong("List of data resources", style="color:dimgrey"),
               " go here.")
         )
      )
   );
   assign("dataresources_guide",
      value=dataresources_guide,
      envir=sashimi_env);
   nbsp <- HTML("&nbsp;");
   nbsp3 <- HTML("&nbsp;&nbsp;&nbsp;");
   assign("nbsp",
      value=nbsp,
      envir=sashimi_env);
   assign("nbsp3",
      value=nbsp3,
      envir=sashimi_env);

   invisible(sashimi_env);
}
