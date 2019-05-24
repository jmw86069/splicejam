
#' SALSA Shiny app constants
#'
#' SALSA Shiny app constants
#'
#' This function is intended to define constant values used in
#' the creation of the SALSA shiny UI.
#'
#' This function is intended to be called only from within
#' an R-shiny app, since it defines several variables in the
#' parent (global) environment.
#'
#' Specifically, this function will attempt to re-use several
#' R data objects by name if they exist in any parent environment.
#' If they do not exist, defaults will be used based upon the
#' Farris et all data from the `farrisdata` package.
#'
#' * filesDF
#' * color_sub
#' * txdb (used to derive exonsByTx, cdsByTx if either does not exist)
#' * tx2geneDF
#' * gtf (used if tx2geneDF or exonsByTx are not available)
#' * exonsByTx
#' * cdsByTx
#' * detectedTx (derived either from farrisdata::farrisTxSE
#'    `"TxDetectedByTPM"` or using all entries in `tx2geneDF$transcript_id`)
#' * detectedGenes (inferred using tx2geneDF and detectedTx)
#' * flatExonsByGene (defined using exonsByTx, cdsByTx, detectedTx,
#'    tx2geneDF)
#' * flatExonsByTx (defined using exonsByTx, cdsByTx, detectedTx,
#'    tx2geneDF)
#'
#' @family SALSA Shiny functions
#'
#' @import shiny
#' @import shinydashboard
#' @import htmltools
#'
#' @export
sashimiAppConstants <- function
(...)
{
   ## Define filesDF here
   if (!exists("filesDF")) {
      if (suppressPackageStartupMessages(require(farrisdata))) {
         printDebug("Using filesDF from ",
            "farrisdata::farris_sashimi_files_df");
         data(farris_sashimi_files_df);
         filesDF <<- farris_sashimi_files_df;
      }
   }

   ## Define color_sub
   if (!exists("color_sub")) {
      color_sub <<- farrisdata::colorSub;
   }

   ## Define flat exons by gene
   ## One-time setup cost when using GTF input
   if (!exists("tx2geneDF") || !exists("exonsByTx")) {
      if (!exists("gtf")) {
         # use default GTF file if not defined
         gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz";
         printDebug("Defining default gtf file:", gtf,
            " from which ",
            c("tx2geneDF, exonsByTx, and cdsByTx"),
            " will be derived.");
      }
      gtfBase <- basename(gtf);
      if (!file.exists(gtfBase)) {
         printDebug("Downloading gtf:", gtf,
            " to:", gtfBase);
         curl::curl_download(url=gtf,
            destfile=gtfBase);
      }
      if (!exists("tx2geneDF")) {
         tx2geneFile <- gsub("[.]gtf.*$",
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
               sep="\t", quote=FALSE, na="", col.names=TRUE);
         } else {
            printDebug("Reloading stored tx2geneDF:",
               tx2geneFile);
            tx2geneDF <- read.table(file=tx2geneFile,
               sep="\t",
               header=TRUE,
               quote="\"",
               comment.char="");
         }
         # Assign in the parent environment
         tx2geneDF <<- tx2geneDF;
      } else {
         printDebug("Using tx2geneDF from environment.");
      }
      ## Now make TxDb
      if (!exists("exonsByTx") || !exists("cdsByTx")) {
         localDb <- gsub("[.]gtf.*$",
            ".txdb",
            gtfBase,
            ignore.case=TRUE);
         if (!file.exists(localDb)) {
            printDebug("Deriving txdb from gtf:", gtfBase);
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
         txdb <<- txdb;
      }
      if (!exists("exonsByTx")) {
         # First obtain exons by transcript
         printDebug("Deriving exonsByTx from txdb");
         require(GenomicFeatures);
         exonsBy_m <- memoise::memoise(exonsBy,
            cache=memoise::cache_filesystem("exonsBy_memoise"));
         exonsBy_m_cached <- memoise::has_cache(exonsBy_m)(
            txdb,
            by="tx",
            use.names=TRUE);
         printDebug("exonsBy_m_cached:",
            exonsBy_m_cached);
         if (!exonsBy_m_cached) {
            printDebug("exonsBy_m_cached digest(txdb):", digest::digest(txdb));
         }
         exonsByTx <- exonsBy_m(
            txdb,
            by="tx",
            use.names=TRUE);
         values(exonsByTx@unlistData)$feature_type <- "exon";
         values(exonsByTx@unlistData)$subclass <- "exon";
         # Assign in the parent environment
         exonsByTx <<- exonsByTx;
      }
      if (!exists("cdsByTx")) {
         printDebug("Deriving cdsByTx from txdb");
         require(GenomicFeatures);
         cdsBy_m <- memoise::memoise(cdsBy,
            cache=memoise::cache_filesystem("cdsBy_memoise"));
         cdsBy_m_cached <- memoise::has_cache(cdsBy_m)(
            txdb,
            by="tx",
            use.names=TRUE);
         printDebug("cdsBy_m_cached:",
            cdsBy_m_cached);
         if (!cdsBy_m_cached) {
            printDebug("cdsBy_m_cached digest(txdb):", digest::digest(txdb));
         }

         cdsByTx <- cdsBy_m(
            txdb,
            by="tx",
            use.names=TRUE);
         values(cdsByTx@unlistData)$feature_type <- "cds";
         values(cdsByTx@unlistData)$subclass <- "cds";
         # Assign in the parent environment
         cdsByTx <<- cdsByTx;
      }
   }

   ## Define detectedTx
   if (!exists("detectedTx") || length(detectedTx) == 0) {
      if (suppressPackageStartupMessages(require(farrisdata))) {
         data(farrisTxSE);
         printDebug("Using detectedTx from farrisTxSE");
         detectedTx <- subset(SummarizedExperiment::rowData(farrisTxSE),
            TxDetectedByTPM)$transcript_id;
      }
      if (!all(detectedTx %in% tx2geneDF$transcript_id)) {
         printDebug("Using detectedTx <- unique(tx2geneDF$transcript_id)");
         detectedTx <- unique(tx2geneDF$transcript_id);
      }
      # Assign in the parent environment
      detectedTx <<- detectedTx;
   }
   ## Infer available genes
   if (!exists("detectedGenes") || length(detectedGenes) == 0) {
      printDebug("Inferring detectedGenes from ",
         format(big.mark=",", length(detectedTx)),
         " detectedTx entries:");
      printDebug("head(tx2geneDF, 3):");
      print(head(tx2geneDF, 3));
      detectedGenes <- jamba::mixedSort(
         unique(
            subset(tx2geneDF,transcript_id %in% detectedTx)$gene_name));
      # Assign in the parent environment
      detectedGenes <<- detectedGenes;
   }

   ## Define flatExonsByGene
   ## define memoised function
   flattenExonsBy_m <- memoise::memoise(flattenExonsBy,
      cache=memoise::cache_filesystem("flattenExonsBy_memoise"));

   if (!exists("flatExonsByGene")) {
      printDebug("Deriving flatExonsByGene from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      flattenExonsByGene_m_cached <- memoise::has_cache(flattenExonsBy_m)(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         detectedTx=detectedTx,
         by="gene",
         tx2geneDF=tx2geneDF,
         verbose=FALSE);
      printDebug("flattenExonsByGene_m_cached:",
         flattenExonsByGene_m_cached);
      if (!flattenExonsByGene_m_cached) {
         printDebug(" digest(exonsByTx): ", digest::digest(exonsByTx));
         printDebug("   digest(cdsByTx): ", digest::digest(cdsByTx));
         printDebug("digest(detectedTx): ", digest::digest(detectedTx));
         printDebug(" digest(tx2geneDF): ", digest::digest(tx2geneDF));
      }
      flatExonsByGene <- flattenExonsBy_m(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         detectedTx=detectedTx,
         by="gene",
         tx2geneDF=tx2geneDF,
         verbose=FALSE);
      # Assign in the parent environment
      flatExonsByGene <<- flatExonsByGene;
   }
   ## Define flatExonsByTx
   if (!exists("flatExonsByTx")) {
      printDebug("Deriving flatExonsByTx from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      flattenExonsByTx_m_cached <- memoise::has_cache(flattenExonsBy_m)(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         tx2geneDF=tx2geneDF,
         by="tx",
         verbose=FALSE);
      printDebug("flattenExonsByTx_m_cached:",
         flattenExonsByTx_m_cached);
      if (!flattenExonsByTx_m_cached) {
         printDebug(" digest(exonsByTx): ", digest::digest(exonsByTx));
         printDebug("   digest(cdsByTx): ", digest::digest(cdsByTx));
         printDebug("digest(detectedTx): ", digest::digest(detectedTx));
         printDebug(" digest(tx2geneDF): ", digest::digest(tx2geneDF));
      }
      flatExonsByTx <- flattenExonsBy_m(
         exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         tx2geneDF=tx2geneDF,
         by="tx",
         verbose=FALSE)
      # Assign in the parent environment
      flatExonsByTx <<- flatExonsByTx;
   }
   if (!exists("verbose")) {
      verbose <- FALSE;
   }

   # guides
   # define guides tab
   guidesTab <<- fluidPage(
      tags$style(type="text/css", "a{color:steelblue; font-weight:bold}"),
      sidebarLayout(
         mainPanel(
            width=7,
            tabBox(
               width=12,
               tabPanel(
                  title="What is a Sashimi Plot?",
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
                     "The methods were developed in support of:"),
                  br(),
                  a("S. Farris, J. M. Ward, K.E. Carstens, M. Samadi, Y. Wang and S. M. Dudek, 2019:",
                     em("Hippocampal subregions express distinct dendritic transcriptomes that reveal an unexpected role for enhanced mitochondrial function in CA2"),
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
            )
         )
      )
   );

   # sashimiplot_guide
   sashimiplot_guide <<- fluidPage(
      h1("What is a Sashimi Plot?",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$p("The Sashimi Plot tab provides a visualization of
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
   sashimiplotviz_guide <<- fluidPage(
      h1("Creating a Sashimi Plot",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$p("The typical workflow for viewing a Sashimi plot is described
            below:"),
         tags$h2("Select the Sashimi Plot tab", style="color:navy"),
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

   dataresources_guide <<- fluidPage(
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
   nbsp <- HTML("&nbsp;");
   nbsp3 <- HTML("&nbsp;&nbsp;&nbsp;");

   invisible(NULL);
}
