
#' Get coverage for GRanges from bigWig files
#'
#' Get coverage for GRanges from bigWig files
#'
#' This function takes a GRanges object to define genomic regions,
#' for which coverage data is loaded for each `bwUrls` input file.
#'
#' Note that this function uses `rtracklayer::import.bw()` which
#' they describe does not work on the Windows platform.
#'
#' Update in version 0.0.68.900: This function was updated in two
#' subtle ways, to work around a bug in `rtracklater::import.bw()`,
#' which returns data sorted by chromosome in the order it is indexed
#' in the bigWig file, then within each chromosome entries are returned
#' in the order requested. This `getGRcoverageFromBw()` was updated to:
#'
#' 1. Confirm input `gr` GRanges contains names, or assigns names
#' as needed.
#' 2. The output coverage from `rtracklayer::import.gw()` is ordered
#' by `names(gr)` to confirm the output coverage is returned in the
#' identical order as requested.
#'
#' The updates above were done outside the scope of memoise file caching,
#' so that stored coverage cache files will still be valid, but the
#' order of named entries from the cache will be dependent upon the
#' order requested. In the event the cache coverage contains no names,
#' they will be returned in the same order as stored, however it is
#' possible the cache will be invalidated by the addition of names
#' to `gr`, though unclear exactly how deeply memoise checks such things.
#'
#' @return DataFrame object, whose colnames are defined using
#'    either `names(bwUrls)` or by `jamba::makeNames(basename(bwUrls))`
#'    then removing the `.bw` or `.bigWig` file extension,
#'    case-insensitively.
#'    Each column is type `IRanges::NumericList-class` which is a
#'    list of numeric coverage values.
#'
#'
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#'
#' @param gr GRanges object
#' @param bwUrls character vector of full file paths or web URLs
#'    to bigWig files, suitable for use by `rtracklayer::import()`.
#' @param addGaps logical indicating whether gaps between GRanges
#'    should be added to the query. Gaps are determined using
#'    `getGRgaps()`. Practically, when `addGaps=TRUE` loads the
#'    coverage data between exons, which can be a substantially
#'    larger region than exons. When `addGaps=FALSE` the coverage
#'    data is not loaded in intron/gap regions, and therefore is
#'    not displayed in downstream plots like sashimi plots.
#' @param feature_type_colname,gap_feature_type,default_feature_type
#'    When `addGaps=TRUE` a
#'    new column named using `feature_type_colname` is added to `values(gr)`,
#'    whose value for gap regions is `gap_feature_type`. When
#'    `feature_type_colname` is already present in `gr` it is not modified,
#'    otherwise the column is created with value `default_feature_type`.
#'    By default, this function adds a column `"feature_type"` with
#'    value `"gap"`.
#' @param use_memoise logical indicating whether to use `memoise::memoise()`
#'    to store coverage data in cache files, which can be re-used in
#'    subsequent R sessions, given consistent values for
#'    `memoise_coverage_path`.
#'    Note that the primary reason to use memoise during this step,
#'    is that cache files will be stored *for each bigWig file* and
#'    not for the set of bigWig files. For example, adding one bigWig
#'    file to `bwUrls` will cause creating of one new memoise cache file,
#'    but will re-use any pre-existing memoise cache files for the
#'    previously cached `bwUrls` entries.
#' @param memoise_coverage_path character path to file folder
#'    used to store coverage data in memoise cache files.
#'    By default, the folder is a subfolder of
#'    the current working directory (see `getwd()`) so it should be
#'    changed to an absolute path if needed for wider re-use in any
#'    working directory.
#' @param do_shiny_progress logical indicating whether to update
#'    shiny progress bar, using `shiny::setProgress()`. It assumes
#'    the progress bar is already initiated.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
getGRcoverageFromBw <- function
(gr,
 bwUrls,
 addGaps=FALSE,
 gap_feature_type="gap",
 default_feature_type="exon",
 feature_type_colname="feature_type",
 use_memoise=FALSE,
 memoise_coverage_path="coverage_memoise",
 do_shiny_progress=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to get coverage from bigWig files,
   ## returning GRanges with columns containing NumericList
   ## get coverage across the regions of interest.
   ##
   ## TODO: consider splitting bwUrls into bwUrlsPos, bwUrlsNeg
   ## in order to allow strand-specificity
   if (!jamba::igrepHas("GRanges", class(gr))) {
      stop("gr must be a GRanges object.");
   }
   if (length(names(bwUrls)) == 0) {
      names(bwUrls) <- jamba::makeNames(
         gsub("[.](bw|bigWig)$",
            "",
            ignore.case=TRUE,
            basename(bwUrls)));
   }
   default_feature_type <- head(c(default_feature_type, "exon"), 1);
   gap_feature_type <- head(c(gap_feature_type, "gap"), 1);
   if (addGaps) {
      if (verbose) {
         jamba::printDebug("getGRcoverageFromBw(): ",
            "addGRgaps()");
      }
      newValues <- list(feature_type=gap_feature_type);
      names(newValues)[1] <- feature_type_colname;
      if (!feature_type_colname %in% colnames(GenomicRanges::values(gr))) {
         GenomicRanges::values(gr)[[feature_type_colname]] <- default_feature_type;
      }
      gr <- addGRgaps(gr,
         newValues=newValues,
         ...);
   }
   import_or_null <- function
   (bwUrl,
      gr) {
      cov1 <- tryCatch({
         rtracklayer::import(bwUrl,
            selection=rtracklayer::BigWigSelection(gr),
            as="NumericList");
      }, error=function(e){
         ## Note: errors occur most commonly when the file is not available
         warnText <- paste0("getGRcoverageFromBw():",
            "import_or_null() error:",
            "BigWig file not accessible:'",
            bwUrl,
            "', returning NULL.");
         jamba::printDebug(warnText);
         # print the error
         print(e);
         # print any associated warnings
         print(warnings());
         warning(warnText);
         NULL;
      });
      cov1;
   }
   if (use_memoise) {
      import_or_null_m <- memoise::memoise(import_or_null,
         cache=memoise::cache_filesystem(memoise_coverage_path));
   }
   ## version 0.0.68.900
   ## - fix issue with rtracklayer::import.bw() returning in bigWig
   ##    chromosome index order and not the input order
   ## - confirm names(gr) exist and are unique
   if (length(names(gr)) == 0) {
      names(gr) <- paste0("GR", seq_along(gr));
   }
   grnames <- jamba::nameVector(names(gr));
   names(gr) <- names(grnames);

   ## Iterate bwUrls and get coverage from each
   if (verbose > 1) {
      jamba::printDebug("getGRcoverageFromBw(): ",
         "bwUrls:");
      print(bwUrls);
   }
   covL <- lapply(jamba::nameVectorN(bwUrls), function(iBw){
      bwUrl <- bwUrls[[iBw]];
      iBwNum <- match(iBw, names(bwUrls));
      iBwPct <- iBwNum / length(bwUrls);
      if (do_shiny_progress && !is.na(iBwPct)) {
         shiny::setProgress(
            value=0/2 + iBwPct/2,
            detail=paste0("Importing coverage (",
               iBwNum,
               " of ",
               length(bwUrls),
               ")"));
      }
      if (verbose) {
         jamba::printDebug("getGRcoverageFromBw(): ",
            "Importing bwUrl:",
            bwUrl);
      }
      if (use_memoise) {
         if (verbose) {
            jamba::printDebug("bwUrl:");
            print(bwUrl);
            cov_has_cache <- memoise::has_cache(import_or_null_m)(
               bwUrl,
               gr=gr);
            jamba::printDebug("   cov_has_cache:", cov_has_cache);
         }
         cov1 <- import_or_null_m(
            bwUrl,
            gr=gr);
         if (length(cov1) == 0) {
            # Repeat once for NULL cached results
            if (verbose) {
               jamba::printDebug("getGRcoverageFromBw(): ",
                  "Repairing coverage cache.",
                  fgText=c("darkorange", "seagreen2"));
            }
            if (do_shiny_progress && !is.na(iBwPct)) {
               shiny::setProgress(
                  value=0/2 + iBwPct/2,
                  detail=paste0("Repairing cov cache (",
                     iBwNum,
                     " of ",
                     length(bwUrls),
                     ")"));
            }
            cov_has_cache <- memoise::has_cache(import_or_null_m)(
               bwUrl,
               gr=gr);
            if (cov_has_cache) {
               cov1 <- tryCatch({
                  memoise::drop_cache(import_or_null_m)(
                     bwUrl,
                     gr=gr);
                  import_or_null_m(
                     bwUrl,
                     gr=gr);
               }, error=function(e){
                  jamba::printDebug("getGRcoverageFromBw(): ",
                     "Error calling memoise::drop_cache(), calling import_or_null() directly.");
                  import_or_null(
                     bwUrl,
                     gr=gr);
               });
            }
            if (length(cov1) == 0) {
               jamba::printDebug("getGRcoverageFromBw(): ",
                  "Failed to repair coverage cache for bwUrl: ",
                  c("'", bwUrl, "'"), sep="",
                  fgText=c("darkorange", "red"));
               if (do_shiny_progress && !is.na(iBwPct)) {
                  shiny::setProgress(
                     value=0/2 + iBwPct/2,
                     detail=paste0("Failed repair cov cache (",
                        iBwNum,
                        " of ",
                        length(bwUrls),
                        ")"));
               }
            }
         }
      } else {
         cov1 <- import_or_null(bwUrl,
            gr=gr);
      }
      if (length(cov1) > 0 && length(names(cov1)) > 0) {
         cov1 <- cov1[match(names(grnames), names(cov1))]
         names(cov1) <- grnames;
      }
      cov1;
   });
   cov1nonzero <- (lengths(covL) > 0);
   if (any(cov1nonzero)) {
      GenomicRanges::values(gr)[,names(covL)[cov1nonzero]] <- S4Vectors::DataFrame(covL[cov1nonzero]);
   }
   if (any(!cov1nonzero)) {
      attr(gr, "some_null") <- TRUE;
   }
   return(gr);
}

#' Combine GRanges coverage replicates
#'
#' Combine GRanges coverage replicates
#'
#' This function takes a GRanges object as output from
#' `getGRcoverageFromBw()` and combines the coverages into
#' one coverage per strand for each `covName` (equivalent
#' to `sample_id`). Each coverage value is multiplied by
#' its `scaleFactors` value, then the sum is returned for
#' each strand, for each `covName` (`sample_id`).
#'
#' The strand is inferred by the
#' presence of negative values, where any negative value
#' indicates the column is negative strand.
#'
#' @return GRanges object whose colnames contain the
#'    `covName` (`sample_id`) for each observed strand,
#'    with coverage combined taking the sum of individual
#'    coverages after multiplying each by `scaleFactors`.
#'
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#'
#' @param gr `GRanges` object containing coverage data in columns
#'    containing NumericList class data.
#' @param covNames `character` vector of `colnames(GenomicRanges::values(gr))`
#'    representing columns in `gr` that contain coverage data
#'    in NumericList format, for example data prepared
#'    with `getGRcoverageFromBw()`.
#' @param covName `character` vector with length equal to
#'    `length(covNames)` representing the `sample_id` for each
#'    `covNames` entry.
#' @param strands `character` vector, or NULL, indicating the strand
#'    for which the coverage data was obtained. When `NULL` the strand
#'    is inferred by the presence of any negative values.
#' @param scaleFactors `numeric` vector length equal to `length(covNames)`
#'    or expanded to that length. Values are multiplied by each
#'    coverage data result, intended to apply a normalization
#'    to each coverage value. A `-1` value can also be used
#'    to flip the score of negative strand data, in the event the
#'    source coverage data is scored only using positive values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
combineGRcoverage <- function
(gr,
 covNames=NULL,
 covName=NULL,
 strands=NULL,
 scaleFactors=1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to combine replicate coverage values
   ##
   ## covName is used to group together replicates of a sample_id
   ## covNames is a unique name per coverage

   ## define scaleFactors=1 if not supplied
   if (length(scaleFactors) == 0) {
      scaleFactors <- 1;
   }
   ## ensure scaleFactors has length=length(covNames)
   scaleFactors <- rep(scaleFactors,
      length.out=length(covNames));
   names(scaleFactors) <- covNames;

   ## define covName="cov" if not supplied
   if (length(covName) == 0) {
      covName <- "cov";
   }
   ## ensure covName has length=length(covNames)
   covName <- rep(covName,
      length.out=length(covNames));
   names(covName) <- covNames;

   ## determine which covNames are present in the supplied gr object
   ## to be tolerant of missing data
   has_coverage <- covNames %in% colnames(GenomicRanges::values(gr));
   some_null <- attr(gr, "some_null");
   if (!all(has_coverage)) {
      jamba::printDebug("covName:", covName);
      jamba::printDebug("covNames:", covNames);
      jamba::printDebug("head(gr):");
      print(head(gr));
      if (verbose) {
         print(data.frame(covName=covName,
            covNames=covNames,
            has_coverage=has_coverage));
      }
      some_null <- TRUE;
      ## re-order
      covName <- covName[has_coverage];
      covNames <- covNames[has_coverage];
      scaleFactors <- scaleFactors[has_coverage];
   }
   #if (length(covName) == 0) {
   #   stop("combineGRcoverage() found no covNames present in colnames(GenomicRanges::values(gr)).");
   #}
   if (length(covName) > 0) {
      if (verbose) {
         jamba::printDebug("combineGRcoverage(): ",
            "scaleFactors: ",
            paste0(names(scaleFactors),
               ":",
               format(digits=2, trim=TRUE, scaleFactors)),
            sep=", ");
      }

      ## define strands if not supplied
      if (length(strands) == 0) {
         strands <- factor(sapply(seq_along(covNames), function(i){
            iCov <- covNames[[i]];
            ifelse(any(
               any(GenomicRanges::values(gr)[[iCov]] * scaleFactors[i] < 0) &
                  all(GenomicRanges::values(gr)[[iCov]] * scaleFactors[i] <= 0)),
               "-",
               "+")
         }), levels=c("+", "-"));
      }
      ## ensure strands is a factor with proper order "+" then "-"
      if (!"factor" %in% class(strands)) {
         strands <- factor(strands,
            levels=jamba::provigrep(c("[+]|pos|plus", "[-]|neg|minus", "."),
               unique(strands)));
      }
      covNamesL <- split(covNames, paste0(covName, strands));
      covNameL <- split(covName, paste0(covName, strands));
      for (iCovN in names(covNamesL)) {
         iCovV <- covNamesL[[iCovN]];
         if (verbose) {
            jamba::printDebug("combineGRcoverage(): ",
               "iCovN:", iCovN);
            jamba::printDebug("combineGRcoverage(): ",
               "   iCovV:", iCovV);
         }
         iCovX <- Reduce("+",
            lapply(iCovV, function(i){
               GenomicRanges::values(gr)[[i]] * scaleFactors[i]
            }));
         if (verbose) {
            jamba::printDebug("combineGRcoverage(): ",
               "iCovN:",
               iCovN);
         }
         GenomicRanges::values(gr)[[iCovN]] <- iCovX;
      }
   }
   ## Remove the original coverage data columns from the output
   keep_gr_colnames <- setdiff(colnames(GenomicRanges::values(gr)),
      covNames);
   gr <- gr[, keep_gr_colnames];

   if (length(covName) > 0) {
      if (exists("covNamesL")) {
         attr(gr, "covNames") <- names(covNamesL);
      }
      attr(gr, "covName") <- jamba::cPasteU(covNameL);
   }
   attr(gr, "some_null") <- some_null;
   return(gr);
}
