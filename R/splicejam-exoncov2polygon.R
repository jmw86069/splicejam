#
# exoncov2polygon()

#' Convert exon coverage to polygons
#'
#' Convert exon coverage to polygons
#'
#' This function is a workhorse function that converts a GRanges
#' object containing column values with NumericList coverage data,
#' into a full data.frame sufficient to define ggplot2 and other
#' coverage polygon plots.
#'
#' An interesting argument is `baseline` which allows each exon
#' in the `gr` GRanges object to be offset from zero, in order to
#' make certain features visually easier to distinguish.
#'
#' This function also calls `simplifyXY()` which reduces the
#' stored polygon detail for regions whose coordinates are compressed
#' on the x-axis, taking roughly the max value for each point.
#'
#' The default output is roughly similar to `broom::tidy()` in
#' that it converts a custom R object into a tidy data.frame
#' suitable for use by ggplot2 and other tidy workflows.
#'
#' The function `getGRcoverageFromBw()` takes a set of bigWig files
#' and returns a GRanges object whose columns contain NumericList data,
#' which is the intended input for `exoncov2polygon()`.
#'
#' @family jam GRanges functions
#' @family jam RNA-seq functions
#' @family splicejam core functions
#'
#' @param gr GRanges where `colnames(GenomicRanges::values(gr))` is present in `covNames`,
#'    and contains data with class `NumericList`.
#' @param covNames character vector contained in `colnames(GenomicRanges::values(gr))`.
#' @param baseline numeric vector of length 0, 1 or `length(gr)`. If
#'    `baseline` has names matching `names(gr)` they will be used for
#'    each `gr` entry; if `baseline` is not named, it is extended
#'    to `length(gr)`. The `baseline` value is added to the coverage
#'    for each exon to offset the polygon as needed.
#' @param gapWidth numeric value sent to `make_ref2compressed()`.
#' @param coord_style character value to define the output style:
#'    `"base"` returns a matrix with polygons separated by a row of `NA`
#'    values; `"fortify"` returns a `data.frame` intended for ggplot2,
#'    with columns `"cov"` and `"gr"` indicating the values in `covNames` and
#'    `names(gr)` used to separate each polygon.
#' @param ref2c optional list containing output from `make_ref2compressed()`,
#'    used to compress the GRanges coordinates.
#' @param compress_introns logical indicating whether to compress
#'    the coverage polygon coordinates to approximately the same
#'    number of pixels per inch as the exon polygons. This option
#'    greatly reduces the size of the polygon, since introns are
#'    already about 50 to 100 times wider than exons, and when
#'    `ref2c` is supplied, the introns are visibly compressed
#'    to a fixed width on the x-axis. The data has many more
#'    x-axis coordinates than the data visualization, this argument
#'    is intended to reduce the intron coordinates accordingly.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @seealso `test_cov_wide_gr()` for examples
#'
#' @returns output dependent upon `coord_style`:
#'    * `"fortify"` returns `data.frame` in tall format, sufficient
#'    to use with `ggplot2` functions. Each polygon is separated
#'    by rows where `x,y` values are `NA`.
#'    * `"list"` returns a `list` with `numeric` matrix objects.
#'    * `"base"` returns a `list` with `numeric` matrix objects,
#'    each of which contains `NA` as the last coordinate, so that
#'    each `matrix` can be `rbind()` and used for vectorized plotting.
#'    * `"all"` returns a `list` with all the data formats.
#'
#' @examples
#' # use some test data
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' data(test_cov_gr);
#' # prepare polygon coordinates
#' exondf <- exoncov2polygon(test_cov_gr, covNames="sample_A");
#' # create a ggplot
#' gg3 <- ggplot(exondf,
#'       aes(x=x, y=y, group=gr, fill=gr, color=gr)) +
#'    ggforce::geom_shape(alpha=0.8) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    colorjam::scale_color_jam();
#' print(gg3);
#' 
#' @export
exoncov2polygon <- function
(gr,
 covNames=NULL,
 sample_id=NULL,
 baseline=NULL,
 gapWidth=250,
 doPlot=FALSE,
 coord_style=c("fortify", "base", "list", "all"),
 ref2c=NULL,
 compress_introns=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take exon coverage in the form of NumericList
   ## associated with GRanges, and produce a list of polygons
   ## with baseline=0
   ##
   ## Workflow:
   ## - input GRanges with coverages in the values columns
   ##   stored as NumericList
   ## create polygons
   suppressPackageStartupMessages(library(GenomicRanges));
   coord_style <- match.arg(coord_style);

   if (!jamba::igrepHas("granges", class(gr))) {
      stop("Input gr should be a GRanges object.");
   }
   if (length(covNames) == 0) {
      covNames <- jamba::provigrep(c("pos|[+]$", "neg|[-]$"),
         colnames(GenomicRanges::values(gr)));
   }
   if (length(covNames) == 0) {
      # Condition occurs when no coverage is available,
      # often when bigwig files are not accessible at all.
      stop_msg <- paste0("Empty covnames suggests coverage data ",
         "are not available or do not match ",
         "colnames(GenomicRanges::values(gr)): ",
         jamba::cPaste(colnames(GenomicRanges::values(gr))));
      if (verbose) {
         jamba::printDebug(stop_msg);
      }
      stop(stop_msg);
      ## 0.0.84.900: return NULL instead of stop()
      ## to allow subsequent steps to succeeed.
      ## Edge case: Let it display only junctions, without coverage.
      return(NULL);
   }
   retVals <- list();

   ## Define compressed coordinate space
   #ref2c <- make_ref2compressed(
   #   subset(gr, feature_type %in% "exon"),
   #   gapWidth=gapWidth);

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   baselineV <- jamba::nameVector(
      rep(0,
         length.out=length(gr)),
      names(gr));
   if (length(baseline) > 0) {
      if (length(names(baseline)) > 0) {
         baselineV[names(baseline)] <- baseline;
      } else {
         baselineV[] <- baseline;
      }
   }

   ## List of lists of polygons
   if (verbose) {
      jamba::printDebug("exoncov2polygon(): ",
         "covNames:",
         covNames);
   }
   covPolyL <- lapply(jamba::nameVector(covNames), function(iName){
      polyL <- lapply(jamba::nameVectorN(gr), function(iGRname){
         yVals1a <- unlist(GenomicRanges::values(gr[iGRname])[[iName]]);
         iBase <- baselineV[iGRname];
         yVals1 <- yVals1a + iBase;
         ## Note: we define xVals1 width using yVals1 width,
         ## because this GRanges might have compressed coordinates
         ## and therefore the width(gr) is not an accurate measure
         ## of the actual width of coverage data
         xVals1 <- seq(from=GenomicRanges::start(gr[iGRname]),
            to=GenomicRanges::end(gr[iGRname]),
            length.out=length(yVals1));
         xVals <- c(rep(head(xVals1, 1), 2) - 0.5,
            xVals1,
            rep(tail(xVals1, 1), 2) + 0.5,
            head(xVals1, 1) - 0.5);
         yVals <- c(iBase,
            head(yVals1, 1),
            yVals1,
            tail(yVals1, 1),
            iBase,
            iBase);
         if (length(xVals) != length(yVals) && verbose) {
            jamba::printDebug("length(xVals):", length(xVals));
            jamba::printDebug("length(yVals):", length(yVals));
         }
         ## TODO: compress coordinates where y-value doesn't change
         if (length(xVals) > 0 && length(yVals) > 0) {
            xy <- cbind(x=xVals, y=yVals);
            ## Simplify the coordinates, typically removing up to 90% rows
            xy <- simplifyXY(xy,
               restrictDegrees=c(0,180));
         } else {
            xy <- cbind(x=0, y=0)[0,,drop=FALSE];
         }
      });
   });
   if ("list" %in% coord_style) {
      return(covPolyL);
   }
   retVals$covPolyL <- covPolyL;

   if ("base" %in% coord_style) {
      ## coordinate matrix suitable for vectorized polygon() by separating
      ## each polygon with a row of NA values
      covPolyML <- lapply(covPolyL, function(iL){
         tail(jamba::rbindList(lapply(iL, function(iM){
            rbind(cbind(x=NA, y=NA),
               iM)
         })), -1);
      });
      return(covPolyML);
   }
   if ("fortify" %in% coord_style) {
      covPolyML <- lapply(jamba::nameVectorN(covPolyL), function(iN){
         iL <- covPolyL[[iN]];
         jamba::rbindList(lapply(jamba::nameVectorN(iL), function(iN2){
            iM <- iL[[iN2]];
            data.frame(check.names=FALSE,
               stringsAsFactors=FALSE,
               iM,
               cov=iN,
               gr=iN2);
         }));
      });
      retVals$covPolyML <- covPolyML;

      ## Compress polygon
      if (length(ref2c) > 0 && compress_introns) {
         t1 <- Sys.time();
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2c=ref2c,
               coord_style=coord_style,
               verbose=verbose);
         });
         t2 <- Sys.time();
         if (verbose > 1) {
            jamba::printDebug("exoncov2polygon(): ",
               "Completed polygon compression:",
               format(t2 - t1));
         }
         retVals$compCovPolyML <- compCovPolyML;
         #return(compCovPolyML);
         covPolyML <- compCovPolyML;
      }

      covPolyDF <- jamba::rbindList(covPolyML);
      covPolyDF$cov <- factor(covPolyDF$cov,
         levels=unique(c(covNames, covPolyDF$cov)));
      covPolyDF$gr <- factor(covPolyDF$gr,
         levels=unique(c(names(gr),
            covPolyDF$gr)));
      ## Optionally add sample_id
      if (length(sample_id) > 0) {
         cov2sample <- jamba::nameVector(sample_id, covNames);
         covPolyDF$sample_id <- cov2sample[covPolyDF$cov];
      }
      return(covPolyDF);
   }

      ## Compress polygon
      if (1 == 2 && length(ref2c) > 0 && compress_introns) {
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2compressed=ref2c,
               coord_style=coord_style);
         });
         retVals$compCovPolyML <- compCovPolyML;
      }

   ## Optionally plot coverage
   if (doPlot) {
      if ("base" %in% coord_style) {
         if (exists("compCovPolyML")) {
            opar <- par("mfrow"=c(2,1));
            on.exit(par(opar));
         }
         polyCol <- colorjam::rainbowJam(length(covPolyL[[1]]));
         plot(jamba::rbindList(covPolyML),
            pch=".",
            xaxt="n",
            col="transparent");
         polygon(jamba::rbindList(covPolyML),
            border=polyCol,
            col=polyCol);
         if (exists("compCovPolyML")) {
            plot(jamba::rbindList(compCovPolyML),
               pch=".",
               xaxt="n",
               col="transparent");
            polygon(compCovPolyML[[1]],
               border=polyCol,
               col=polyCol);
            polygon(compCovPolyML[[2]],
               border=polyCol,
               col=polyCol);
         }
      }
   }
   return(retVals);
}
