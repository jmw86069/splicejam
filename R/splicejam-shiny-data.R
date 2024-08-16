
#' Prepare sashimi plot required data
#'
#' Prepare sashimi plot required data, deriving data objects as needed
#'
#' This function performs a subset of steps performed by
#' `sashimiAppConstants()`, focusing only on data required
#' for gene-exon structure. The `sashimiAppConstants()` defines
#' `color_sub` and validates `filesDF`, then calls this function
#' `sashimiDataConstants()` to prepare and validate the gene-exon
#' data.
#'
#' Data derived by this function `sashimiDataConstants()`:
#'
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
#' @family splicejam R-shiny functions
#'
#' @return `environment` that contains the required data objects
#'    for splicejam sashimi plots. Note that the environment itself
#'    is updated during processing, so the environment does not
#'    need to be returned for the data contained inside it to
#'    be updated by this function.
#'
#' @param gtf,txdb,tx2geneDF,exonsByTx,cdsByTx objects used to define
#'    the overall set of genes, transcripts, and associated exons and
#'    CDS exons. See this function
#'    description for more detail.
#' @param detectedTx,detectedGenes,flatExonsByGene,flatExonsByTx
#'    objects used to derive a specific subset of gene-exon models
#'    using only detected transcripts or genes. See this function
#'    description for more detail.
#' @param default_gene `character` string indicating the default
#'    gene to use for the initial R-shiny figure.
#' @param envir `environment` where data will be prepared, or when
#'    `envir=NULL` a new environment will be created and returned.
#' @param empty_uses_farrisdata `logical` indicating whether to
#'    use data from the Github R package `"jmw86069/farrisdata"`
#'    if no data is supplied to this function. This behavior is
#'    intended to make it easy to use farrisdata to recreate
#'    the Sashimi plots in that publication.
#' @param use_memoise `logical` indicating whether to use `memoise`
#'    to cache intermediate data files for exons, flattened exons,
#'    transcript-gene data, and so on. This mechanism reduces
#'    time to render sashimi plots that re-use the same gene.
#'    All memoise cache folders are named with `"_memoise"`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @import data.table
#'
#' @export
sashimiDataConstants <- function
(gtf=NULL,
 txdb=NULL,
 tx2geneDF=NULL,
 exonsByTx=NULL,
 cdsByTx=NULL,
 detectedTx=NULL,
 detectedGenes=NULL,
 flatExonsByGene=NULL,
 flatExonsByTx=NULL,
 envir=NULL,
 empty_uses_farrisdata=TRUE,
 use_memoise=TRUE,
 verbose=FALSE,
 ...)
{
   if (!is.environment(envir)) {
      if (verbose) {
         jamba::printDebug("sashimiDataConstants(): ",
            "Assigning new.env() because envir was empty");
      }
      envir <- new.env(parent=emptyenv());
   }

   params <- setdiff(names(formals(sashimiDataConstants)),
      c("...", "envir"));
   for (param in params) {
      if (verbose > 1) {
         jamba::printDebug("sashimiDataConstants(): ",
            "Assigning param: ", param);
      }
      assign(x=param,
         value=get_fn_envir(param,
            envir=envir,
            verbose=verbose - 1),
         envir=envir)
   }
   rm(list=params);
   #envir <- attach(envir, name="sashimi_env");
   #on.exit(detach(name="sashimi_env"));


   # Logic flow:
   # - If tx2geneDF OR exonsByTx are empty
   #    - if gtf and txdb are empty
   #       - if empty_uses_farrisdata then use farrisdata gtf
   #       - otherwise STOP
   #    - if gtf exists
   #       - download gtf
   #       - if tx2geneDF is empty
   #          - create tx2geneDF, save to cache / or load from cache
   #          - assign tx2geneDF
   #    - if exonsByTx is empty
   #       - if txdb is empty
   #          - create txdb from gtf, save to cache / or load from cache
   #          - assign txdb
   #       - create exonsByTx from txdb, using memoise
   #       - assign exonsByTx
   #    - if cdsByTx is empty and txdb exists
   #       - create cdsByTx from txdb, using memoise
   #       - assign cdsByTx
   #
   # - if detectedTx is empty
   #    - if empty_uses_farrisdata, and farrisdata detectedTx matches tx2genedF, use farrisdata
   #    - define detectedTx as all transcripts
   #    - assign detectedTx
   #
   # - if detectedGenes is empty
   #    - define detectedGenes using tx2geneDF and detectedTx
   #    - assign detectedGenes
   #
   # - if flatExonsByGene is empty
   #    - define flatExonsByGene with flattenExonsBy() with exonsByTx, cdsBytx, detectedTx, tx2geneDF
   #    - assign flatExonsByGene
   #
   # - if flatExonsByTx is empty
   #    - define flatExonsByTx with flattenExonsBy() with exonsByTx, cdsBytx, detectedTx, tx2geneDF
   #    - assign flatExonsByTx

   msg <- function(...) {
      jamba::printDebug("sashimiDataConstants(): ",
         ...);
   }
   ## Define flat exons by gene
   ## One-time setup cost when using GTF input
   if (length(envir$tx2geneDF) == 0 || length(envir$exonsByTx) == 0) {
      if (verbose) msg("length(tx2geneDF) == 0 || length(exonsByTx) == 0");
      if ((!exists("gtf", envir=envir) || length(envir$gtf) == 0) && length(envir$txdb) == 0) {
         if (envir$empty_uses_farrisdata) {
            # use default GTF file if not defined
            envir$gtf <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz";
            jamba::printDebug("sashimiDataConstants(): ",
               "Defining default gtf file: '",
               envir$gtf,
               "' from which ",
               c("tx2geneDF", " exonsByTx", " and cdsByTx"),
               " will be derived.");
         } else {
            stop(paste0("The 'gtf' or 'txdb' argument are required ",
               "when either 'tx2geneDF' or 'exonsByTx' are not provided."));
         }
      } else {
         envir$empty_uses_farrisdata <- FALSE
      }
      if (length(envir$gtf) == 0 && length(envir$txdb) == 0) {
         stop(paste0("The 'gtf' or 'txdb' argument are required ",
            "when either 'tx2geneDF' or 'exonsByTx' are not provided."));
      }
      ## tx2geneDF
      if (length(envir$gtf) > 0) {
         gtfBase <- basename(envir$gtf);
         if (!file.exists(gtfBase)) {
            jamba::printDebug("sashimiDataConstants(): ",
               "Downloading gtf: '", envir$gtf,
               "' to: '", gtfBase, "'");
            curl::curl_download(url=envir$gtf,
               destfile=gtfBase);
         }
         if (length(envir$tx2geneDF) == 0) {
            tx2geneFile <- gsub("[.](gff|gff3|gtf).*$",
               ".tx2geneDF.txt",
               gtfBase,
               ignore.case=TRUE);
            if (!file.exists(tx2geneFile)) {
               if (verbose) {
                  jamba::printDebug("sashimiDataConstants(): ",
                     "Deriving tx2geneDF from gtf: '",
                     gtfBase,
                     "' then storing: '",
                     tx2geneFile, "'");
               }
               envir$tx2geneDF <- makeTx2geneFromGtf(GTF=gtfBase,
                  verbose=FALSE);
               data.table::fwrite(x=envir$tx2geneDF,
                  file=tx2geneFile,
                  sep="\t",
                  quote=FALSE,
                  na="",
                  row.names=FALSE,
                  col.names=TRUE);
            } else {
               jamba::printDebug("sashimiDataConstants(): ",
                  "Reloading stored tx2geneDF: '",
                  tx2geneFile, "'");
               envir$tx2geneDF <- data.table::fread(file=tx2geneFile,
                  sep="\t",
                  data.table=FALSE);
            }
            # Assign in the appropriate environment
            #assign("tx2geneDF",
            #   value=tx2geneDF,
            #   envir=envir);
         } else if (verbose) {
            jamba::printDebug("sashimiDataConstants(): ",
               "Using tx2geneDF as supplied.");
         }
      }
      if (length(envir$tx2geneDF) == 0) {
         stop("The 'tx2geneDF' argument is required when 'gtf' is not supplied.");
      }

      ## Now make TxDb in order to derive exonsByTx and cdsByTx
      #if (length(exonsByTx) == 0 || length(cdsByTx) == 0) {
      if (length(envir$exonsByTx) == 0) {
         if (length(envir$txdb) > 0) {
            jamba::printDebug("sashimiDataConstants(): ",
               "Using supplied txdb.");
         } else {
            localDb <- gsub("[.](gff|gff3|gtf).*$",
               ".txdb",
               gtfBase,
               ignore.case=TRUE);
            if (!file.exists(localDb)) {
               jamba::printDebug("sashimiDataConstants(): ",
                  "Deriving txdb from gtf: '",
                  gtfBase,
                  "' to store as: '",
                  localDb, "'");
               envir$txdb <- GenomicFeatures::makeTxDbFromGFF(gtfBase);
               AnnotationDbi::saveDb(x=envir$txdb,
                  file=localDb);
            } else {
               jamba::printDebug("sashimiDataConstants(): ",
                  "Reloading txdb from: '", localDb, "'");
               envir$txdb <- AnnotationDbi::loadDb(file=localDb);
            }
            if (!DBI::dbIsValid(AnnotationDbi::dbconn(envir$txdb))) {
               if (verbose) {
                  jamba::printDebug("sashimiDataConstants(): ",
                     "Refreshing db connection: '", localDb, "'");
               }
               envir$txdb <- AnnotationDbi::loadDb(file=localDb);
            }
         }

         # First obtain exons by transcript
         jamba::printDebug("sashimiDataConstants(): ",
            c("Deriving ","exonsByTx"," from txdb"), sep="");
         #suppressPackageStartupMessages(require(GenomicFeatures));
         if (use_memoise) {
            exonsBy_m <- memoise::memoise(GenomicFeatures::exonsBy,
               cache=memoise::cache_filesystem("exonsBy_memoise"));
            exonsBy_m_cached <- memoise::has_cache(exonsBy_m)(
               envir$txdb,
               by="tx",
               use.names=TRUE);
            jamba::printDebug("exonsBy_m_cached:",
               exonsBy_m_cached);
         } else {
            exonsBy_m <- GenomicFeatures::exonsBy;
         }
         envir$exonsByTx <- exonsBy_m(
            envir$txdb,
            by="tx",
            use.names=TRUE);
         GenomicRanges::values(envir$exonsByTx@unlistData)$feature_type <- "exon";
         GenomicRanges::values(envir$exonsByTx@unlistData)$subclass <- "exon";
      } else {
         if (!"feature_type" %in% names(GenomicRanges::values(envir$exonsByTx@unlistData))) {
            GenomicRanges::values(envir$exonsByTx@unlistData)$feature_type <- "exon";
         }
         if (!"subclass" %in% names(GenomicRanges::values(envir$exonsByTx@unlistData))) {
            GenomicRanges::values(envir$exonsByTx@unlistData)$subclass <- "exon";
         }
      }

      ## create cdsByTx
      if (length(envir$cdsByTx) == 0 && exists("txdb", envir=envir, inherits=FALSE)) {
         jamba::printDebug("sashimiDataConstants(): ",
            "Deriving cdsByTx from txdb.");
         #suppressPackageStartupMessages(require(GenomicFeatures));
         if (use_memoise) {
            cdsBy_m <- memoise::memoise(GenomicFeatures::cdsBy,
               cache=memoise::cache_filesystem("cdsBy_memoise"));
            cdsBy_m_cached <- memoise::has_cache(cdsBy_m)(
               envir$txdb,
               by="tx",
               use.names=TRUE);
            jamba::printDebug("cdsBy_m_cached:",
               cdsBy_m_cached);
         } else {
            cdsBy_m <- GenomicFeatures::cdsBy;
         }

         envir$cdsByTx <- cdsBy_m(
            envir$txdb,
            by="tx",
            use.names=TRUE);
         GenomicRanges::values(envir$cdsByTx@unlistData)$feature_type <- "cds";
         GenomicRanges::values(envir$cdsByTx@unlistData)$subclass <- "cds";
      } else {
         if (!"feature_type" %in% names(GenomicRanges::values(envir$cdsByTx@unlistData))) {
            GenomicRanges::values(envir$cdsByTx@unlistData)$feature_type <- "exon";
         }
         if (!"subclass" %in% names(GenomicRanges::values(envir$cdsByTx@unlistData))) {
            GenomicRanges::values(envir$cdsByTx@unlistData)$subclass <- "exon";
         }
      }
   } else {
      if (verbose) msg("tx2geneDF or exonsByTx were provided.");
   }

   ## Define detectedTx
   if (length(envir$detectedTx) == 0) {
      if (verbose) msg("length(detectedTx) == 0");
      if (envir$empty_uses_farrisdata && nchar(system.file(package="farrisdata")) > 0) {
         #data(farrisTxSE);
         envir$detectedTx <- subset(SummarizedExperiment::rowData(farrisdata::farrisTxSE),
            TxDetectedByTPM)$transcript_id;
         if (verbose) {
            jamba::printDebug("sashimiDataConstants(): ",
               c("Using ",
                  jamba::formatInt(length(envir$detectedTx)),
                  " detectedTx from '", "farrisdata::farrisTxSE", "'"),
               sep="");
         }
      } else {
         envir$detectedTx <- unique(envir$tx2geneDF$transcript_id);
         if (verbose) {
            jamba::printDebug("sashimiDataConstants(): ",
               c("Defined ",
                  jamba::formatInt(length(envir$detectedTx)),
                  " detectedTx using '",
                  "tx2geneDF$transcript_id",
                  "' since no detectedTx were supplied."),
               sep="");
         }
      }
   }
   ## Confirm detectedTx are present in tx2geneDF
   ## - if none are present, use all from tx2geneDF
   if (!all(envir$detectedTx %in% envir$tx2geneDF$transcript_id)) {
      detlen <- length(unique(envir$detectedTx));
      envir$detectedTx <- intersect(envir$detectedTx,
         envir$tx2geneDF$transcript_id);
      if (length(envir$detectedTx) == 0) {
         envir$detectedTx <- unique(envir$tx2geneDF$transcript_id);
         if (verbose) {
            jamba::printDebug("sashimiDataConstants(): ",
               c("None of the ",
                  jamba::formatInt(detlen),
                  " detectedTx entries were found in '",
                  "tx2geneDF$transcript_id",
                  " therefore ",
                  "detectedTx",
                  " will use all of ",
                  "tx2geneDF$transcript_id"),
               sep="");
         }
         if (length(envir$detectedTx) == 0) {
            jamba::printDebug("sashimiDataConstants(): ",
               c("No values were present in ",
                  "tx2geneDF$transcript_id",
                  ". Printing ",
                  "head(tx2geneDF, 20)", ":"),
               sep="");
            print(head(envir$tx2geneDF, 20));
            stop("There were no values in tx2geneDF$transcript_id");
         }
      } else {
         if (verbose) {
            jamba::printDebug("sashimiDataConstants(): ",
               c("Matched ",
                  jamba::formatInt(detlen),
                  " entries in ", "detectedTx", " to ",
                  jamba::formatInt(length(envir$detectedTx)),
                  " entries in ",
                  "tx2geneDF$transcript_id"),
               sep="");
         }
      }
   }

   ## Infer available genes
   if (length(envir$detectedGenes) == 0) {
      envir$detectedGenes <- jamba::mixedSort(
         unique(
            subset(envir$tx2geneDF,
               transcript_id %in% envir$detectedTx)$gene_name));
      if (verbose) {
         jamba::printDebug("sashimiDataConstants(): ",
            c("Inferred ",
               jamba::formatInt(length(envir$detectedGenes)),
               " detectedGenes from ",
               jamba::formatInt(length(envir$detectedTx)),
               " detectedTx."),
            sep="");
      }
      if (length(envir$detectedGenes) == 0) {
         stop("No detectedGenes were found using detectedTx and tx2geneDF$gene_name");
      }
   }

   ## Define flatExonsByGene
   ## define memoised function
   if (use_memoise) {
      flattenExonsBy_m <- memoise::memoise(flattenExonsBy,
         cache=memoise::cache_filesystem("flattenExonsBy_memoise"));
   } else {
      flattenExonsBy_m <- flattenExonsBy;
   }


   if (length(envir$flatExonsByGene) == 0) {
      jamba::printDebug("sashimiDataConstants(): ",
         "Deriving flatExonsByGene using: ",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      if (verbose && use_memoise) {
         jamba::printDebug("sdim(envir):");
         print(jamba::sdim(envir));
         flattenExonsByGene_m_cached <- memoise::has_cache(flattenExonsBy_m)(
            exonsByTx=envir$exonsByTx,
            cdsByTx=envir$cdsByTx,
            detectedTx=envir$detectedTx,
            by="gene",
            tx2geneDF=envir$tx2geneDF,
            verbose=FALSE);
         jamba::printDebug("sashimiDataConstants(): ",
            "flattenExonsByGene_m_cached:",
            flattenExonsByGene_m_cached);
      }
      envir$flatExonsByGene <- flattenExonsBy_m(
         exonsByTx=envir$exonsByTx,
         cdsByTx=envir$cdsByTx,
         detectedTx=envir$detectedTx,
         by="gene",
         tx2geneDF=envir$tx2geneDF,
         verbose=FALSE);
   }

   ## Define flatExonsByTx
   if (length(envir$flatExonsByTx) == 0) {
      jamba::printDebug("sashimiDataConstants(): ",
         "Derived flatExonsByTx using: ",
         c("exonsByTx", "cdsByTx", "detectedTx", "tx2geneDF"));
      if (verbose && use_memoise) {
         flattenExonsByTx_m_cached <- memoise::has_cache(flattenExonsBy_m)(
            exonsByTx=envir$exonsByTx,
            cdsByTx=envir$cdsByTx,
            detectedTx=envir$detectedTx,
            tx2geneDF=envir$tx2geneDF,
            by="tx",
            verbose=FALSE);
         jamba::printDebug("sashimiDataConstants(): ",
            "flattenExonsByTx_m_cached:",
            flattenExonsByTx_m_cached);
      }
      envir$flatExonsByTx <- flattenExonsBy_m(
         exonsByTx=envir$exonsByTx,
         cdsByTx=envir$cdsByTx,
         detectedTx=envir$detectedTx,
         tx2geneDF=envir$tx2geneDF,
         by="tx",
         verbose=FALSE)
   }

   ## Define default_gene
   if (length(envir$default_gene) == 0 || nchar(envir$default_gene) == 0) {
      envir$default_gene <- head(
         jamba::provigrep(c("Gria1",
            "Ntrk2",
            "Actb",
            "Gapd",
            "^[A-Z][a-z]{3}",
            "^[A-Za-z]{4}",
            "^[A-Za-z]{3}",
            "."),
            detectedGenes),
         1);
   }

   return(envir)
}

#' Get value from function call or specific environment
#'
#' Get value from function call or specific environment in that order
#'
#' This function is a helper function intended to return a
#' variable value, if it exists and is not NULL, by searching
#' these locations in order:
#'
#' 1. The calling function, which is the environment of the
#' function that called `get_fn_envir()`.
#' 2. The environment or environments provided in `envir`.
#' 3. It returns `NULL` if the previous steps do not find
#' the object named by `x`.
#'
#' @return object represented by variable name given in `x` from either
#'    the calling function, or the environment `envir`, or `NULL`
#'    if not defined in either case.
#'
#' @param x `character` string indicating the name of an R object.
#' @param envir `environment` or `list` of `environment` objects.
#' @param assign_to_envir `logical` indicating whether to assign values
#'    to environment `envir`, if `envir` is not NULL. This option is
#'    helpful to combine function arguments with environment values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- 10;
#' get_fn_envir("x")
#'
#' test_x <- function(x=NULL, envir=NULL, verbose=FALSE, ...) {
#'    get_fn_envir("x", envir, verbose=verbose)
#' }
#'
#' test_x()
#'
#' test_x(envir=globalenv())
#'
#' test_x(x=5)
#'
#' test_x(x=5, envir=globalenv())
#' test_x(x=NULL, envir=globalenv())
#'
#' test_x(envir=globalenv())
#'
#' # create new environment
#' testenv <- new.env();
#' testenv$x <- 100;
#'
#' test_x(envir=testenv, verbose=TRUE)
#' test_x(x=1000, envir=testenv, verbose=TRUE)
#'
#' # search testenv then globalenv()
#' test_x(x=12, envir=c(testenv, globalenv()), verbose=TRUE)
#' test_x(envir=c(testenv, globalenv()), verbose=TRUE)
#'
#' testenv$x <- NULL;
#' test_x(envir=c(testenv, globalenv()), verbose=TRUE)
#'
#' rm("x", envir=testenv);
#' test_x(envir=c(testenv, globalenv()), verbose=TRUE)
#'
#' @export
get_fn_envir <- function(x,
 envir=NULL,
 verbose=FALSE,
 ...)
{
   #jamba::printDebug("get_fn_envir_value(): ", "search():\n", search());
   if (exists(x, envir=parent.frame(1)) && length(get(x, envir=parent.frame(1))) > 0) {
      if (verbose) {
         jamba::printDebug("get_fn_envir(): ",
            "Variable '", x, "' found in ", "parent.frame(1)");
      }
      return(get(x, envir=parent.frame(1)))
   } else if (length(envir) > 0) {
      if (!is.list(envir) && is.environment(envir)) {
         envir <- list(envir);
      }
      for (i in seq_along(envir)) {
         if (is.environment(envir[[i]]) &&
               exists(x, envir=envir[[i]]) &&
               length(get(x, envir=envir[[i]])) > 0) {
            if (verbose) {
               jamba::printDebug("get_fn_envir(): ",
                  "Variable '", x, "' found in the provided ",
                  c("envir",
                     ifelse(length(envir) > 1, paste0("[[", i, "]]"), "")),
                  sep="");
            }
            return(get(x, envir=envir[[i]]))
         }
      }
   } else {
      if (verbose) {
         jamba::printDebug("get_fn_envir(): ",
            "Variable '", x, "' not found");
      }
   }
   return(NULL);
}
