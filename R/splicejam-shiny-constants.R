
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
#' * tx2geneDF
#' * gtf (used if tx2geneDF or exonsByTx are not available)
#' * exonsByTx
#' * cdsByTx
#' * detectedTx
#' * detectedGenes, inferred from tx2geneDF and detectedTx
#' * flatExonsByGene
#' * flatExonsByTx
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
         gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz";
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
            printDebug("Deriving tx2geneDF from gtf:", gtfBase);
            tx2geneDF <- makeTx2geneFromGtf(GTF=gtfBase,
               verbose=FALSE);
            write.table(file=tx2geneFile,
               x=tx2geneDF,
               sep="\t", quote=FALSE, na="", col.names=TRUE);
         } else {
            printDebug("Reloading tx2geneDF:", tx2geneFile);
            tx2geneDF <- read.table(file=tx2geneFile,
               sep="\t", header=TRUE, quote="\"", comment.char="");
         }
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
      }
      if (!exists("exonsByTx")) {
         # First obtain exons by transcript
         printDebug("Deriving exonsByTx from txdb");
         require(GenomicFeatures);
         exonsBy_m <- memoise::memoise(exonsBy,
            cache=memoise::cache_filesystem("exonsBy_memoise"));
         exonsByTx <- exonsBy_m(txdb,
            by="tx",
            use.names=TRUE);
         values(exonsByTx@unlistData)$feature_type <- "exon";
         values(exonsByTx@unlistData)$subclass <- "exon";
         exonsByTx <<- exonsByTx;
      }
      if (!exists("cdsByTx")) {
         printDebug("Deriving cdsByTx from txdb");
         require(GenomicFeatures);
         if (!exists("cdsBy_m")) {
            cdsBy_m <- memoise::memoise(cdsBy,
               cache=memoise::cache_filesystem("cdsBy_memoise"));
         }
         cdsByTx <- cdsBy_m(txdb,
            by="tx",
            use.names=TRUE);
         values(cdsByTx@unlistData)$feature_type <- "cds";
         values(cdsByTx@unlistData)$subclass <- "cds";
         cdsByTx <<- cdsByTx;
      }
   }

   ## Define detectedTx
   if (!exists("detectedTx") || length(detectedTx) == 0) {
      if (suppressPackageStartupMessages(require(farrisdata))) {
         data(farrisTxSE);
         printDebug("Using detectedTx from farrisTxSE");
         detectedTx <<- subset(SummarizedExperiment::rowData(farrisTxSE),
            TxDetectedByTPM)$transcript_id;
      }
      if (!all(detectedTx %in% tx2geneDF$transcript_id)) {
         printDebug("Using detectedTx <- unique(tx2geneDF$transcript_id)");
         detectedTx <<- unique(tx2geneDF$transcript_id);
      }
   }
   ## Infer available genes
   if (!exists("detectedGenes") || length(detectedGenes) == 0) {
      printDebug("Inferring detectedGenes from detectedTx:",
         head(detectedTx));
      printDebug("head(tx2geneDF):");
      print(head(tx2geneDF));
      detectedGenes <<- jamba::mixedSort(
         unique(
            subset(tx2geneDF, transcript_id %in% detectedTx)$gene_name));
   }

   ## Define flatExonsByGene
   if (!exists("flatExonsByGene") || !exists("flatExonsByTx")) {
      flattenExonsBy_m <- memoise::memoise(flattenExonsBy,
         cache=memoise::cache_filesystem("flattenExonsBy_memoise"));
   }
   if (!exists("flatExonsByGene")) {
      printDebug("Deriving flatExonsByGene from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      flatExonsByGene <<- flattenExonsBy_m(exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         detectedTx=detectedTx,
         by="gene",
         tx2geneDF=tx2geneDF,
         verbose=FALSE);
   }
   ## Define flatExonsByTx
   if (!exists("flatExonsByTx")) {
      printDebug("Deriving flatExonsByTx from:",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      flatExonsByTx <<- flattenExonsBy_m(exonsByTx=exonsByTx,
         cdsByTx=cdsByTx,
         tx2geneDF=tx2geneDF,
         by="tx",
         verbose=FALSE)
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
                  title="Sashimi Plot",
                  uiOutput("sashimiplot_guide")
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
                     "The methods were developed in support of "),
                  a("S. Farris, J. M. Ward, K.E. Carstens, M. Samadi, Y. Wang and S. M. Dudek, 2019:",
                     em("Hippocampal subregions express distinct dendritic transcriptomes that reveal an unexpected role for enhanced mitochondrial function in CA2"),
                     href="https://github.com/jmw86069/jampack")
               ),
               tags$li(
                  strong(style="color:firebrick",
                     "Sashimi plots were originally envisioned by MISO "),
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
      h1("Sashimi Plot",
         style="color:firebrick"),
      shinydashboard::box(
         width=12,
         status="primary",
         style="background-color:aliceblue",
         tags$ul(
            tags$li(
               strong("Define filesDF", style="color:navy"),
               "- data.frame with bigWig coverage, BED junctions, and sample_id",
               tags$ul(
                  tags$li(
                     strong("url"),
                     " url, file path, or R object name."
                  ),
                  tags$li(
                     strong("type"),
                     " one of 'coverage_gr', 'bw', or 'junction'."
                  ),
                  tags$li(
                     strong("sample_id"),
                     " unique identifier per biological sample."
                  )
               )
            ),
            tags$li(strong("Import gene-exon models", style="color:navy"),
               tags$ul(
                  tags$li(
                     strong("gencode GTF"),
                     " import gencode-compatible GTF file."
                  ),
                  tags$li(
                     strong("TxDb"),
                     " import Bioconductor 'TxDb' package."
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
