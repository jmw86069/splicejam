




#' Create a ref2compressed function to compress GR gaps
#'
#' Create a ref2compressed function to compress GR gaps
#'
#' This function takes a set of GRanges which are to be maintained
#' with fixed aspect ratio, and it defines a function to compress
#' coordinates of the gaps between GRanges features.
#'
#' @family jam GRanges functions
#' @family splicejam core functions
#' @family jam RNA-seq functions
#'
#' @return list with `trans_grc` which is class `"trans"` suitable
#'    for use in ggplot2 functions; `transform` a function that converts
#'    chromosome coordinates to compressed coordinates; `inverse` a function
#'    that converts compressed coordinates to chromosome coordinates;
#'    `scale_x_grc` a function used similar to `ggplot2::scale_x_continuous()`
#'    during ggplot2 creation; `gr` a function that compresses coordinates
#'    in a `GRanges` object; `grl` a function that compresses coordinates in
#'    a `GRangesList` object. Attributes `"lookupCoordDF"` is a two-column
#'    data.frame with chromosome coordinates and compressed coordinates,
#'    which is used to create the other transformation functions via
#'    `stats::approx()`; `"gapWidth"` the gap width used, since it can
#'    be programmatically defined; `"gr"` the `GRanges` input data used
#'    to train the transformation.
#'
#'
#' @param gr GRanges object containing regions not to compress. Regions
#'    which are unstranded gaps are compressed to fixed width.
#' @param gapWidth integer value used for fixed gap width, or when
#'    NULL the gap width is defined as 3 times the median feature width.
#' @param keepValues logical indicating whether to keep feature values
#'    in the GRanges data.
#' @param upstream,downstream,upstreamGapWidth,downstreamGapWidth used
#'    to define the compression of coordinates upstream and downstream
#'    the supplied GRanges. In reality, the upstream range and upstream
#'    gap width defines a multiplier, and all upstream coordinates are
#'    compressed through zero. Similarly, all downstream coordinates
#'    are compressed to 10 billion, which is roughly 3 times the size
#'    of the human genome.
#' @param nBreaks the default number of x-axis coordinate breaks used
#'    in ggplot labeling.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @seealso `grl2df()`, `test_junc_wide_gr`
#'
#' @examples
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' data(test_exon_wide_gr);
#' # To plot a simple GRanges object
#' widedf <- grl2df(test_exon_wide_gr);
#' ggWide <- ggplot(widedf, aes(x=x, y=y, group=id, fill=feature_type)) +
#'    geom_polygon() +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    xlab("chr1") +
#'    ggtitle("exons (introns as-is)")
#' print(ggWide);
#'
#' # Now compress the introns keeping axis labels
#' ref2c <- make_ref2compressed(test_exon_wide_gr,
#'    nBreaks=10);
#' ggWide2 <- ggWide +
#'    scale_x_continuous(trans=ref2c$trans_grc) +
#'    xlab("chr1 (compressed introns)") +
#'    ggtitle("exons (compressed introns)")
#' print(ggWide2);
#'
#' @export
make_ref2compressed <- function
(gr,
 gapWidth=200,
 keepValues=FALSE,
 upstream=50000,
 upstreamGapWidth=gapWidth*3,
 downstream=50000,
 downstreamGapWidth=gapWidth*3,
 nBreaks=7,
 verbose=FALSE,
...)
{
   ## Purpose is to use a GRanges object to train a coordinate
   ## compression function, which assigns the gaps to a fixed
   ## width, but allows GRanges features to keep their original
   ## width.
   if (length(gapWidth) == 0) {
      gapWidth <- max(c(10,
         round(median(width(reduce(gr))) / 2)));
      if (verbose) {
         printDebug("make_ref2compressed(): ",
            "determined gapWidth:",
            format(gapWidth,
               scientific=FALSE,
               big.mark=","));
      }
   }
   if (length(upstreamGapWidth) == 0) {
      upstreamGapWidth <- gapWidth * 3;
   }
   if (length(downstreamGapWidth) == 0) {
      downstreamGapWidth <- gapWidth * 3;
   }

   ## Define the exon-exon distance, to see if they are adjacent or not
   #disGR <- disjoin(gr);
   disGR <- reduce(gr);
   if (is.null(names(disGR))) {
      names(disGR) <- makeNames(rep("disGR", length(disGR)));
   }
   if (length(disGR) > 1) {
      exonGaps <- sapply(head(seq_along(disGR), -1), function(i1){
         i2 <- i1 + 1;
         distance(disGR[i1], disGR[i2]) > 0;
      });
      names(exonGaps) <- head(names(disGR), -1);
   } else {
      exonGaps <- NULL;
   }

   ## Now reset coordinates using fixed gap width
   #newWidths <- head(intercalate(width(disGR),
   #   (gapWidth)*exonGaps), -1);
   newWidths <- suppressWarnings(
      head(intercalate(width(disGR),
         (gapWidth)*exonGaps),
         length(disGR)*2-1)
      );

   ##
   if (verbose) {
      printDebug("make_ref2compressed(): ",
         "newCoords");
   }
   newCoords <- cumsum(newWidths);
   newCoordsM <- matrix(c(0, newCoords), ncol=2, byrow=TRUE);
   newCoordsM[,1] <- newCoordsM[,1] + 1;

   lookupCoordDF <- jamba::mixedSortDF(unique(data.frame(
      refCoord=c(start(disGR),end(disGR)),
      coord=c(newCoordsM[,1],newCoordsM[,2])
   )));
   if (upstream > 0) {
      ## Add one upstream point
      refCoord1 <- head(lookupCoordDF$refCoord,1);
      coord1 <- head(lookupCoordDF$coord,1);
      upExtended1 <- (refCoord1 - upstream);
      upComp1 <- (coord1 - upstreamGapWidth);
      lookupCoordDF <- rbind(
         data.frame(refCoord=upExtended1,
            coord=upComp1),
         lookupCoordDF);

      ## Extend the upstream gap to 1
      minRefCoord <- 1;
      if (upExtended1 > minRefCoord) {
         upExtended1 <- minRefCoord;
         newRefDiff <- (refCoord1 - minRefCoord);
         upComp1 <- coord1 - (upstreamGapWidth * newRefDiff) / upstream;
         lookupCoordDF <- rbind(
            data.frame(refCoord=upExtended1,
               coord=upComp1),
            lookupCoordDF);
      }
   }
   if (downstream > 0) {
      ## Add one downstream point
      refCoord2 <- tail(lookupCoordDF$refCoord,1);
      coord2 <- tail(lookupCoordDF$coord,1);
      downExtended2 <- (refCoord2 + downstream);
      downComp2 <- (coord2 + downstreamGapWidth);
      lookupCoordDF <- rbind(
         lookupCoordDF,
         data.frame(refCoord=downExtended2,
            coord=downComp2)
      );
      ## Extend the upstream gap to 10 Gb
      maxRefCoord <- 3e10;
      if (downExtended2 < maxRefCoord) {
         downExtended2 <- maxRefCoord;
         newRefDiff <- (maxRefCoord - refCoord2);
         downComp2 <- coord2 + (downstreamGapWidth * newRefDiff) / downstream;
         lookupCoordDF <- rbind(
            lookupCoordDF,
            data.frame(refCoord=downExtended2,
               coord=downComp2)
         );
      }
   }
   ## TODO: expand the range of coordinates so approxfun() will not fail
   if (verbose) {
      printDebug("make_ref2compressed(): ",
         "approxfun");
   }
   ref2compressed <- approxfun(x=lookupCoordDF[,1],
      y=lookupCoordDF[,2],
      method="linear");
   compressed2ref <- approxfun(x=lookupCoordDF[,2],
      y=lookupCoordDF[,1],
      method="linear");

   ## Add a convenience function for GRanges objects
   if (verbose) {
      printDebug("make_ref2compressed(): ",
         "ref2compressedGR");
   }
   ref2compressedGR <- function(gr) {
      if (length(names(gr)) > 0) {
         GRnames <- names(gr);
      }
      if (all(c("refStart","refEnd") %in% colnames(values(GR)))) {
         ## Re-compress the original reference coordinates
         ranges(gr) <- IRanges(
            start=values(gr)[,"refStart"],
            end=values(gr)[,"refEnd"]
         );
      } else {
         values(gr)[,"refStart"] <- start(gr);
         values(gr)[,"refEnd"] <- end(gr);
      }
      ranges(gr) <- IRanges(
         start=ref2compressed(start(gr)),
         end=ref2compressed(end(gr))
      );
      if (length(GRnames) > 0) {
         names(gr) <- GRnames;
      }
      return(gr);
   }
   ## Add a convenience function for GRangesList objects
   ref2compressedGRL <- function(grl) {
      if (length(names(grl)) == 0) {
         names(grl) <- makeNames(rep("grl", length(grl)));
      }
      grlc1 <- ref2compressedGR(grl@unlistData);
      grlc <- GenomicRanges::split(grlc1,
         factor(
            rep(names(grl), elementNROWS(grl)),
            levels=names(grl)));
      grlc;
   }
   if (verbose) {
      printDebug("make_ref2compressed(): ",
         "retVals");
   }

   ## Custom breaks function
   breaks_gr <- function(x, limits=NULL, n=nBreaks, xFixed=lookupCoordDF[,1], verbose=FALSE, ...) {
      xvals <- unique(sort(xFixed));
      xvals <- xvals[xvals >= min(x) & xvals <= max(x)];
      if (verbose) {
         printDebug("breaks_gr(): ",
            "x:", x);
         printDebug("breaks_gr(): ",
            "limits:", limits);
         printDebug("breaks_gr(): ",
            "n:", n);
      }
      if (n > length(xvals)) {
         return(xvals);
      }
      idx1 <- round(seq.int(from=1, to=length(xvals), length.out=n));
      return(xvals[idx1]);
   }
   minor_breaks <- function(b, limits=NULL, n=2, xFixed=lookupCoordDF[,1], verbose=FALSE, ...) {
      if (verbose) {
         printDebug("minor_breaks(): ",
            "x:", x);
         printDebug("minor_breaks(): ",
            "limits:", limits);
         printDebug("minor_breaks(): ",
            "n:", n);
      }
      nUse <- length(b) * n;
      ref2compressed(breaks_gr(compressed2ref(b),
         limits=compressed2ref(limits),
         n=nUse,
         xFixed=xFixed,
         verbose=verbose));
   }
   labels <- function(x, n=10) {
      scales::comma_format()(breaks_gr(x, n=n))
   }
   retVals <- list();

   ## Make the custom trans function
   trans_grc <- scales::trans_new(name="compressed_gr",
      transform=ref2compressed,
      inverse=compressed2ref,
      breaks=breaks_gr,
      minor_breaks=minor_breaks,
      format=scales::comma_format(),
      domain=range(lookupCoordDF[,1]));

   ## Make the actual scale_x_gr_compressed() function
   scale_x_grc <- function(..., trans=trans_grc){
      scale_x_continuous(name="grc",
         breaks=waiver(),#trans_grc$breaks,
         minor_breaks=waiver(),#trans_grc$minor_breaks,
         labels=function(...){scales::comma(...)},#trans_grc$labels,
         trans=trans,
         ...)
   }

   ## Functions required by scales::trans_new()
   retVals$transform <- ref2compressed;
   retVals$inverse <- compressed2ref;

   ## Convenience functions for GenomicRanges objects
   retVals$gr <- ref2compressedGR;
   retVals$grl <- ref2compressedGRL;

   #retVals$breaks <- breaks_gr;
   #retVals$minor_breaks <- minor_breaks;
   #retVals$format <- scales::comma_format;
   #retVals$labels <- labels;

   ## ggplot2 functions to compress the x-axis
   retVals$scale_x_grc <- scale_x_grc;
   retVals$trans_grc <- trans_grc;

   ## TODO: custom breaks() function that prioritizes choosing labels
   ## at borders of the compressed ranges, then sensible labels between them


   attr(retVals, "lookupCoordDF") <- lookupCoordDF;
   attr(retVals, "gapWidth") <- gapWidth;
   ## GR might be used to allow adding a feature then refreshing the function
   attr(retVals, "gr") <- gr;

   ## TODO: Convert reference coordinates to feature-based coordinates,
   ## e.g. to coordinates relative to the length of a transcript

   ## Note the format to convert any GRanges coordinates
   ## newCoordsGR <- ref2compressedGR(GR,
   ##    ref2compressed=ref2compressed);

   return(retVals);
}

#' Simplify XY coordinates to minimal line segments
#'
#' Simplify XY coordinates to minimal line segments
#'
#' This function takes a numeric matrix of x,y coordinates
#' and returns the minimal matrix of x,y coordinates that represents
#' the same line segments. It is intended in cases where there
#' is a long repeated line segment that could be represented by
#' far fewer points.
#'
#' @param xy numeric matrix with two columns representing x,y coordinates.
#' @param minN integer value to define the minimal number of repeated
#'    values before the compression is used. By definition, three consecutive
#'    points must have the same slope in order for compression to be
#'    effective, otherwise the original coordinates will be returned.
#' @param restrictDegrees numeric vector of degrees to restrict the
#'    simplification. For exmample `restrictDegrees=c(0,180)` will only
#'    simplify horizontal lines.
#' @param ... additional arguments are ignored.
#'
#' @family jam spatial functions
#'
#' @param xy numeric matrix of two columns x,y
#' @param minN integer minimum number of consecutive points to
#'    cause compression to occur. It requires at least 3 points
#'    inherent to the algorithm.
#' @param restrictDegrees numeric vector of allowed angles in degrees
#'    (range 0 to 360) of angles allowed to be compressed, when supplied
#'    all other angles will not be compressed.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' xy <- cbind(
#'    x=c(1,1:15,1),
#'    y=c(0,0,0,0,1,2,3,4,5,5,5,5,6,0,0));
#' par("mfrow"=c(1,1));
#' plot(xy, cex=2);
#' points(simplifyXY(xy), pch=20, col="orange", cex=3);
#' points(simplifyXY(xy, restrictDegrees=c(0,180)), pch=17, col="purple");
#' legend("topleft", pch=c(1, 20, 17),# box.col="black", border="black",
#'    pt.cex=c(2, 3, 1),
#'    col=c("black", "orange", "purple"),
#'    legend=c("all points", "simplified points", "only simplified horizontal"));
#'
#' @export
simplifyXY <- function
(xy,
 minN=3,
 restrictDegrees=NULL,
 ...)
{
   ## Purpose is to take a matrix of x,y coordinates and return
   ## the minimal x,y coordinates that describes the same line segments.
   ## It represents only the first and last point for each consecutive
   ## line segments sharing the same angle.
   xDiff <- xy[,1] - c(head(xy[,1], 1), head(xy[,1], -1));
   yDiff <- xy[,2] - c(head(xy[,2], 1), head(xy[,2], -1));

   ## First compress simple repeated xy coordinates
   if (any(tail(xDiff, -1) == 0 & tail(yDiff, -1) == 0)) {
      repXY <- which(tail(xDiff, -1) == 0 & tail(yDiff, -1) == 0) + 1;
      xy <- xy[-repXY,,drop=FALSE];
      xDiff <- xDiff[-repXY];
      yDiff <- yDiff[-repXY];
   }

   ## Next determine the angle from point to point,
   ## two non-repeated points sharing the same angle must be
   ## on the same line, so we only need the first and last
   ## points to define the segment.
   xyAngle <- atan2(x=xDiff,
      y=yDiff);
   if (length(restrictDegrees) > 0) {
      xyAngle <- rad2deg(xyAngle);
   }
   yRle <- Rle(xyAngle);

   if (any(runLength(yRle) >= minN)) {
      rl <- runLength(yRle);
      rv <- runValue(yRle);
      rlDF <- data.frame(start=cumsum(c(1, head(rl, -1))),
         end=cumsum(rl),
         rl=rl,
         rv=rv);
      if (length(restrictDegrees) > 0) {
         ## && any(rlDF$rv %in% restrictDegrees)
         #
         yWhich <- (rlDF$rv %in% restrictDegrees & rlDF$rl > 1);
         yKeep <- unique(c(1,
            cumsum(
               rep(ifelse(yWhich, rlDF$rl, 1),
                  ifelse(yWhich, 1, rlDF$rl)))
            ));
      } else {
         yKeep <- unique(c(1, rlDF$end));
      }
      xyNew <- xy[yKeep,,drop=FALSE];
      if (1 == 2) {
         rlDF1 <- data.frame(idx=c(
            cumsum(c(1, head(rl, -1))),
            cumsum(rl)),
            rl=rep(rl, 2),
            rv=rep(rv, 2)
         )
         rlDF2 <- unique(mixedSortDF(
            data.frame(idx=c(rlDF[,1], rlDF[,2]),
               rl=rep(rl, 2),
               rv=rep(rv, 2)),
            byCols=1));
         xyNew <- xy[rlDF2$idx,,drop=FALSE];
      }
      xyNew;
   } else {
      xy;
   }
}

#' Compress genome coordinates of a matrix of polygons
#'
#' Compress genome coordinates of a matrix of polygons
#'
#' This function takes a two-column numeric matrix of polygons
#' where the x coordinate is the genomic position, and y coordinate
#' is the coverage. It uses `ref2compressed$transform` to convert
#' coordinates to compressed coordinates, as output from
#' `make_ref2compressed()`.
#'
#' For regions that have been compresssed, it then compresses the
#' y-coordinate information to roughly one y value per integer in
#' compressed coordinate space, using the runmax across the window of
#' coverages compressed to this value.
#'
#' @family jam spatial functions
#'
#' @return data.frame with the same colnames, with reduced rows
#'    for polygons where the coordinate compression defined in
#'    `ref2c` is above `minRatio` ratio. The compressed coverage roughly
#'    represents the max coverage value for each point.
#'
#' @param polyM data.frame including `x`,`y` coordinates of polygons,
#'    `cov` and `gr` to group each polygon.
#' @param ref2c list output from `make_ref2compressed()`.
#' @param minRatio the minimum ratio of coordinate compression required
#'    before polygon resolution is downsampled.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
compressPolygonM <- function
(polyM,
 ref2c,
 minRatio=3,
 verbose=FALSE,
 ...)
{
   ## Purpose is to compress coordinates in coverage polygons
   ##
   ## General strategy is to determine the relative scaling of
   ## each polygon, and downsample polygons proportional to the
   ## amount of compression on the x-axis.
   if (any(is.na(polyM[,1]))) {
      data_style <- "base";
      polyMlengths <- diff(unique(c(0, which(is.na(polyM[,1])), nrow(polyM))));
      polyMratios <- shrinkMatrix(polyDF$ratio,
         groupBy=polyDF$n,
         shrinkFunc=function(x){median(x, na.rm=TRUE)})$x;
   } else {
      data_style <- "fortify";
      idrows <- pasteByRow(polyM[,c("cov","gr")]);
      idrows <- factor(idrows, levels=unique(idrows));
      polyML <- split(polyM, idrows);
      polyMratios <- unlist(lapply(polyML, function(i){
         xr1 <- range(i[,1]);
         xc1 <- ref2c$transform(xr1);
         xv1 <- rev(sort(c(diff(xr1)+1, diff(xc1)+1)));
         xv1[1] / xv1[2];
      }));
      polyMlengths <- sdim(polyML)[,1];
   }
   if (verbose) {
      printDebug("compressPolygonM(): ",
         "polyMratios:", head(polyMratios, 20));
   }

   ## data.frame describing the compression and polygon
   polyDF <- as.data.frame(polyM);
   polyDF$newX <- ref2c$transform(polyDF$x);
   polyDF$ratio <- c(NA, diff(polyDF$newX) / diff(polyDF$x));
   polyMrepN <- rep(seq_along(polyMlengths), polyMlengths);
   polyDF$n <- polyMrepN;

   polyDF$medianRatio <- rep(polyMratios,
      polyMlengths);

   polyDFL <- split(polyDF, polyMrepN);
   whichComp <- which(polyMratios >= minRatio);
   whichNorm <- which(polyMratios < minRatio);
   if (length(whichComp) == 0) {
      return(polyM);
   }
   ## Compress those polygons whose coordinates get compressed
   polyDFLnew <- lapply(nameVector(whichComp), function(k){
      iDF <- polyDFL[[k]];
      baseline <- iDF[1,"y"];
      iDF <- iDF[!is.na(iDF[,1]),,drop=FALSE];
      iDFu <- iDF[match(unique(iDF$x), iDF$x),,drop=FALSE];

      iRange <- range(iDF$x, na.rm=TRUE);
      iRangeNew <- range(iDF$newX, na.rm=TRUE);
      iMedRatio <- head(iDF$medianRatio, 1);
      iN <- ceiling(diff(iRange)/iMedRatio);
      iMultiple <- ceiling(iMedRatio);
      iSeqNew <- seq(from=iRange[1],
         to=iRange[2],
         by=iMedRatio);
      iSeqSub <- seq(from=iRange[1],
         to=iRange[2],
         length.out=iN);
      #iSeqNew <- seq(from=iRangeNew[1],
      #   to=iRangeNew[2],
      #   by=iMedRatio);
      #iSeqSub <- seq(from=iRangeNew[1],
      #   to=iRangeNew[2],
      #   length.out=iN);
      iSeqSub[iSeqSub > max(iSeqNew)] <- max(iSeqNew);

      ## use approx() to fill in holes
      #iYnew <- approx(x=iDFu$newX,
      iYnew <- approx(x=iDFu$x,
         y=iDFu$y,
         xout=iSeqNew)$y;
      ## Take running max value, then use approx
      iYrunmax <- caTools::runmax(iYnew,
         k=iMultiple,
         endrule="keep");
      if (1 == 2) {
         printDebug("head(iDFu):");print(head(iDFu));
         printDebug("head(iYnew):", head(iYnew));
         printDebug("head(iMultiple):", head(iMultiple));
         printDebug("head(iN):", head(iN));
         printDebug("head(iSeqNew):", head(iSeqNew));
         printDebug("head(iYrunmax):", head(iYrunmax));
         printDebug("head(iSeqSub):", head(iSeqSub));
      }
      iYsub <- approx(x=iSeqNew,
         y=iYrunmax,
         xout=iSeqSub)$y;
      if (head(iYsub, 1) != baseline) {
         iYsub <- c(baseline, iYsub);
         iSeqSub <- c(head(iSeqSub, 1), iSeqSub);
      }
      if (tail(iYsub, 1) != baseline) {
         iYsub <- c(iYsub, baseline);
         iSeqSub <- c(iSeqSub, tail(iSeqSub, 1));
      }
      iM <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         x=iSeqSub,
         y=iYsub,
         gr=head(iDF$gr,1),
         cov=head(iDF$cov, 1));
   });
   newPolyDFL <- list();
   newPolyDFL[names(polyDFL)[whichNorm]] <- lapply(nameVector(whichNorm), function(k){
      polyDFL[[k]][,c("x","y","cov","gr"),drop=FALSE];
   });
   newPolyDFL[names(polyDFLnew)] <- polyDFLnew;
   newPolyDFL <- newPolyDFL[mixedSort(names(newPolyDFL))];
   newPolyM <- rbindList(newPolyDFL);
   #plot(newPolyM, pch=".");
   #polygon(newPolyM, col=rainbowJam(46));
   return(newPolyM);
}

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
#' @param gr GRanges where `colnames(values(gr))` is present in `covNames`,
#'    and contains data with class `NumericList`.
#' @param covNames character vector contained in `colnames(values(gr))`.
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
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @seealso `test_cov_wide_gr()` for examples
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

   if (!igrepHas("granges", class(gr))) {
      stop("Input gr should be a GRanges object.");
   }
   if (length(covNames) == 0) {
      covNames <- provigrep(c("pos|[+]$", "neg|[-]$"),
         colnames(values(gr)));
   }
   if (length(covNames) == 0) {
      stop("covNames must be colnames(values(gr)).");
   }
   retVals <- list();

   ## Define compressed coordinate space
   #ref2c <- make_ref2compressed(
   #   subset(gr, feature_type %in% "exon"),
   #   gapWidth=gapWidth);

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   baselineV <- nameVector(
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
      printDebug("exoncov2polygon(): ",
         "covNames:",
         covNames);
   }
   covPolyL <- lapply(nameVector(covNames), function(iName){
      polyL <- lapply(nameVectorN(gr), function(iGRname){
         yVals1a <- unlist(values(gr[iGRname])[[iName]]);
         iBase <- baselineV[iGRname];
         yVals1 <- yVals1a + iBase;
         ## Note: we define xVals1 width using yVals1 width,
         ## because this GRanges might have compressed coordinates
         ## and therefore the width(gr) is not an accurate measure
         ## of the actual width of coverage data
         xVals1 <- seq(from=start(gr[iGRname]),
            to=end(gr[iGRname]),
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
         if (length(xVals) != length(yVals)) {
            printDebug("length(xVals):", length(xVals));
            printDebug("length(yVals):", length(yVals));
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
         tail(rbindList(lapply(iL, function(iM){
            rbind(cbind(x=NA, y=NA),
               iM)
         })), -1);
      });
      return(covPolyML);
   }
   if ("fortify" %in% coord_style) {
      covPolyML <- lapply(nameVectorN(covPolyL), function(iN){
         iL <- covPolyL[[iN]];
         rbindList(lapply(nameVectorN(iL), function(iN2){
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
      if (length(ref2c) > 0) {
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2c=ref2c,
               coord_style=coord_style);
         });
         retVals$compCovPolyML <- compCovPolyML;
         #return(compCovPolyML);
         covPolyML <- compCovPolyML;
      }

      covPolyDF <- rbindList(covPolyML);
      covPolyDF$cov <- factor(covPolyDF$cov,
         levels=unique(c(covNames, covPolyDF$cov)));
      covPolyDF$gr <- factor(covPolyDF$gr,
         levels=unique(c(names(gr),
            covPolyDF$gr)));
      ## Optionally add sample_id
      if (length(sample_id) > 0) {
         cov2sample <- nameVector(sample_id, covNames);
         covPolyDF$sample_id <- cov2sample[covPolyDF$cov];
      }
      return(covPolyDF);
   }

      ## Compress polygon
      if (length(ref2c) > 0) {
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
         if (exists(compCovPolyML)) {
            par("mfrow"=c(2,1));
         }
         polyCol <- colorjam::rainbowJam(length(covPolyL[[1]]));
         plot(rbindList(covPolyML),
            pch=".",
            xaxt="n",
            col="transparent");
         polygon(rbindList(covPolyML),
            border=polyCol,
            col=polyCol);
         if (exists(compCovPolyML)) {
            plot(rbindList(compCovPolyML),
               pch=".",
               xaxt="n",
               col="transparent");
            polygon(compCovPolyML[[1]],
               border=polyCol,
               col=polyCol);
            polygon(compCovPolyML[[2]],
               border=polyCol,
               col=polyCol);
            par("mfrow"=c(1,1));
         }
      }
   }
   return(retVals);
}

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
#' @return DataFrame object, whose colnames are defined using
#'    `makeNames(basename(bwUrls))`. Each column is type `NumericList`,
#'    which is a list of numeric coverage values.
#'
#' @family jam GRanges functions
#' @family RNA-seq functions
#'
#' @param gr GRanges object
#' @param bwUrls character vector of full file paths or web URLs
#'    to bigWig files, suitable for use by `rtracklayer::import()`.
#' @param addGaps logical indicating whether gaps between GRanges
#'    should be added to the query. Gaps are determined using
#'    `getGRgaps()`.
#' @param feature_type_colname,gap_feature_type,default_feature_type
#'    When `addGaps=TRUE` a
#'    new column named using `feature_type_colname` is added to `values(gr)`,
#'    whose value for gap regions is `gap_feature_type`. When
#'    `feature_type_colname` is already present in `gr` it is not modified,
#'    otherwise the column is created with value `default_feature_type`.
#'    By default, this function adds a column `"feature_type"` with
#'    value `"gap"`.
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
 verbose=FALSE,
 ...)
{
   ## Purpose is to get coverage from bigWig files,
   ## returning GRanges with columns containing NumericList
   ## get coverage across the regions of interest.
   ##
   ## TODO: consider splitting bwUrls into bwUrlsPos, bwUrlsNeg
   ## in order to allow strand-specificity
   if (!igrepHas("GRanges", class(gr))) {
      stop("gr must be a GRanges object.");
   }
   if (length(names(bwUrls)) == 0) {
      names(bwUrls) <- makeNames(
         gsub("[.](bw|bigWig)$",
            "",
            ignore.case=TRUE,
            basename(bwUrls)));
   }
   default_feature_type <- head(c(default_feature_type, "exon"), 1);
   gap_feature_type <- head(c(gap_feature_type, "gap"), 1);
   if (addGaps) {
      if (verbose) {
         printDebug("getGRcoverageFromBw(): ",
            "addGRgaps()");
      }
      newValues <- list(feature_type=gap_feature_type);
      names(newValues)[1] <- feature_type_colname;
      if (!feature_type_colname %in% colnames(values(gr))) {
         values(gr)[[feature_type_colname]] <- default_feature_type;
      }
      gr <- addGRgaps(gr,
         newValues=newValues,
         ...);
   }
   ## Iterate bwUrls and get coverage from each
   covL <- lapply(nameVectorN(bwUrls), function(iBw){
      bwUrl <- bwUrls[[iBw]];
      if (verbose) {
         printDebug("getGRcoverageFromBw(): ",
            "Importing bwUrl:",
            bwUrl);
      }
      cov1 <- tryCatch({
         rtracklayer::import(bwUrl,
            selection=rtracklayer::BigWigSelection(gr),
            as="NumericList");
      }, error=function(e){
         ## Note: errors occur most commonly when the file is not available
         warnText <- paste0("getGRcoverageFromBw(): ",
            "BigWig file not accessible:'",
            iBw,
            "', returning NULL.");
         warning(warnText);
         NULL;
      })
   });
   values(gr)[,names(covL)] <-S4Vectors::DataFrame(covL);
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
#' @param gr GRanges object containing coverage data in columns
#'    containing NumericList class data.
#' @param covNames character vector of `colnames(values(gr))`
#'    that contain coverage in NumericList format.
#' @param covName character vector with length equal to
#'    `length(covNames)` representing the `sample_id` for each
#'    `covNames` entry.
#' @param strands character vector, or NULL, indicating the strand
#'    for which the coverage data was obtained. When `NULL` the strand
#'    is inferred by the presence of any negative values.
#' @param scaleFactors numeric vector length equal to `length(covNames)`
#'    of values to multiply by each coverage result, in order to
#'    normalized coverages to each other.
#' @param verbose logical indicating whether to print verbose output.
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
   if (any(!covNames %in% colnames(values(gr)))) {
      stop("combineGRcoverage() requires covNames be present in colnames(values(gr)).");
   }
   if (length(scaleFactors) == 0) {
      scaleFactors <- 1;
   }
   if (length(covName) == 0) {
      covName <- "cov";
   }
   covName <- rep(covName, length.out=length(covNames));
   names(covName) <- covNames;
   if (length(scaleFactors) != length(covNames)) {
      scaleFactors <- rep(scaleFactors, length.out=length(covNames));
   }
   names(scaleFactors) <- covNames;
   if (length(strands) == 0) {
      strands <- factor(sapply(seq_along(covNames), function(i){
         iCov <- covNames[[i]];
         ifelse(any(any(values(gr)[[iCov]] * scaleFactors[i] < 0) &
               all(values(gr)[[iCov]] * scaleFactors[i] <= 0)),
            "-",
            "+")
      }), levels=c("+", "-"));
   }
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
         printDebug("combineGRcoverage(): ",
            "iCovN:", iCovN);
         printDebug("combineGRcoverage(): ",
            "iCovV:", iCovV);
      }
      iCovX <- Reduce("+",
         lapply(iCovV, function(i){
            values(gr)[[i]] * scaleFactors[i]
         }));
      if (verbose) {
         printDebug("combineGRcoverage(): ",
            "iCovN:",
            iCovN);
      }
      values(gr)[[iCovN]] <- iCovX;
   }
   gr <- gr[,!colnames(values(gr)) %in% covNames];
   attr(gr, "covNames") <- names(covNamesL);
   attr(gr, "covName") <- cPasteUnique(covNameL);
   return(gr);
}

#' Prepare Sashimi plot data
#'
#' Prepare Sashimi plot data
#'
#' This function is the workhorse function used to produce
#' Sashimi plots, and is intended to be a convenient wrapper
#' function for several other individual functions.
#'
#' At a minimum, a Sashimi plot requires three things:
#'
#' 1. Exons, usually from a gene of interest.
#' 2. RNA-seq coverage data.
#' 3. Splice junction data.
#'
#' There is some required pre-processing before running
#' `prepareSashimi()`:
#'
#' * Prepare flattened exons by gene using `flattenExonsByGene()`
#' and corresponding data, including `exonsByGene`, `cdsByGene`,
#' and `tx2geneDF`. Verify the gene exon model data using
#' `gene2gg()`.
#' * Find file paths, or web URLs, for a set of bigWig coverage
#' files, representing RNA-seq coverage for each strand, for
#' the samples of interest. Test the coverage data using
#' `getGRcoverageFromBw()` for a small set of GRanges data.
#' * Find file paths, or web URLs, for a set of BED6 or BED12
#' format files, note that it cannot currently use bigBed format
#' due to limitations in the `rtracklayer` package.
#' Test the splice junction data using `rtracklayer::import()`
#' for a small range of GRanges features, then send the data
#' to `spliceGR2junctionDF()` to prepare a data.frame summary.
#'
#' The basic input for coverage and junction data is a data.frame,
#' which defines each file path or url, the type of data
#' `"bw"` or `"junction"`, and the biological sample `"sample_id"`.
#' Any file path compatible with `rtracklayer::import()` will
#' work, including web URLs and local files. When using a web URL
#' you may need to use `"https://"` format to force the use
#' of secure web requests, but this requirement varies by country.
#'
#' @return list containing `ggSashimi` a ggplot2 graphical object
#'    containing a full Sashimi plot; `ggCov` the RNA-seq coverage
#'    subset of the Sashimi plot; `ggJunc` the splice junction
#'    subsset of the Sashimi plot; `ref2c` the output of
#'    `make_ref2compressed()` used for ggplot2 coordinate
#'    visualization; `covDF`, `juncDF` data.frame objects
#'    with the raw data used to create ggplot2 objects;
#'    `covGR`, `juncGR` the GRanges objects used to create
#'    the data.frames; `gr` the GRanges object representing the
#'    exons for the gene of interest; `juncLabelDF` the data.frame
#'    containing exon label coordinates used to add labels to
#'    the splice junction arcs.
#'
#' @param flatExonsByGene GRangesList named by gene, whose GRanges
#'    elements are flattened, disjoint, non-overlapping genomic ranges
#'    per gene.
#' @param filesDF data.frame with columns `url`, `sample_id`, `type`,
#'    where: `url` is any valid file path or URL compatible with
#'    `base::read.table()`; `sample_id` is an identified representing
#'    a biological sample, used to group common files together;
#'    `type` is one of `"bw"` for bigWig coverage, `"junction"` for
#'    BED12 format splice junctions.
#' @param gene character string of the gene to prepare, which must be
#'    present in `names(flatExonsByGene)`.
#' @param gapWidth numeric value of the fixed width to use for
#'    gaps (introns) between exon features. If `NULL` then
#'    `getGRgaps()` will use the default based upon the median exon
#'    width.
#' @param addGaps logical indicating whether to include gap regions
#'    in the coverage plot, for example including introns or intergenic
#'    regions. When `compressGR=TRUE` then gaps regions are
#'    down-sampled using running maximum signal with roughly the same
#'    x-axis resolution as uncompressed regions.
#' @param baseline numeric vector named by `names(flatExonsByGene)`
#'    where baseline is used to adjust the y-axis baseline position
#'    above or below zero.
#' @param compressGR logical indicating whether to compress GRanges
#'    coordinates in the output data, where gaps/introns are set
#'    to a fixed width. When `ref2c` is not supplied, and
#'    `compressGR=TRUE`, then `ref2c` is created using
#'    `make_ref2compressed()`.
#' @param ref2c list object output from `make_ref2compressed()` used
#'    to compress axis coordinates, to compress polygon coverage
#'    data in compressed regions, and to adjust splice junction arcs
#'    using compressed coordinates.
#' @param gap_feature_type the default feature_type value to use for
#'    gaps when `addGaps=TRUE`.
#' @param doStackJunctions logical indicating whether to stack
#'    junction arcs at each end, this argument is passed to
#'    `grl2df()` which calls `stackJunctions()`.
#' @param covGR GRanges object containing coverage data in columns
#'    stored as NumericList class, where `colnames(values(covGR))`
#'    are present in `filesDF$url` when `files$type %in% "coverage_gr"`.
#' @param juncGR GRanges object containing splice junctions, where
#'    `"score"` is used for the abundance of splice junction reads,
#'    and `"sample_id"` is used to define the biological `sample_id`.
#' @param include_strand character value, one of `"both"`, `"+"`,
#'    `"-"` indicating the strandedness of coverage and junctions
#'    to display. The default `"both"` shows coverage on both strands,
#'    otherwise coverage is filtered either by filename (presence of
#'    `"pos"`, `"+"`, or `"plus"` indicates positive strand), or
#'    by detecting strandedness by positive/negative coverage scores.
#'    Detecting by filename is intended to avoid retrieving coverage
#'    in the R-shiny app, to help efficiency.
#' @param verbose logical indicating whether to print verbose output.
#' @param do_shiny_progress logical indicating whether to send
#'    progress updates to a running shiny app, using the
#'    `shiny::withProgress()` and `shiny::incProgress()` methods.
#'    This function only calls `shiny::incProgress()` and
#'    assumes the `shiny::withProgress()` has already been
#'    initialized.
#' @param ... additional arguments are passed to `make_ref2compressed()`,
#'    `getGRcoverageFromBw()`, `exoncov2polygon()`.
#'
#' @family RNA-seq functions
#' @family jam plot functions
#' @family splicejam core functions
#'
#' @examples
#' # The active example below uses sample data
#' suppressPackageStartupMessages(library(GenomicRanges));
#'
#' data(test_exon_gr);
#' data(test_junc_gr);
#' data(test_cov_gr);
#' filesDF <- data.frame(url="sample_A",
#'    type="coverage_gr",
#'    sample_id="sample_A");
#' sh1 <- prepareSashimi(GRangesList(TestGene1=test_exon_gr),
#'    filesDF=filesDF,
#'    gene="TestGene1",
#'    covGR=test_cov_gr,
#'    juncGR=test_junc_gr);
#' plotSashimi(sh1);
#'
#' @export
prepareSashimi <- function
(flatExonsByGene=NULL,
 filesDF=NULL,
 gene,
 sample_id=NULL,
 minJunctionScore=10,
 gapWidth=200,
 addGaps=TRUE,
 baseline=0,
 compressGR=TRUE,
 ref2c=NULL,
 gap_feature_type="intron",
 default_feature_type="exon",
 feature_type_colname="feature_type",
 exon_label_type=c("none", "repel", "mark"),
 junc_label_type=c("repel", "mark", "none"),
 return_data=c("ggCov", "ggJunc", "ggSashimi", "covDF", "juncDF", "ref2c", "all"),
 include_strand=c("both", "+", "-"),
 junc_color=alpha2col("goldenrod3", 0.7),
 junc_fill=alpha2col("goldenrod1", 0.4),
 doStackJunctions=TRUE,
 coord_method=c("coord", "scale", "none"),
 scoreFactor=1,
 scoreArcFactor=0.2,
 scoreArcMinimum=100,
 covGR=NULL,
 juncGR=NULL,
 do_shiny_progress=FALSE,
 verbose=FALSE,
 ...)
{
   ## Purpose it to wrapper several functions used to prepare various
   ## types of data for Sashimi plots
   ##
   ## TODO: Allow filtering junction scores based upon the exon coverage
   ## so the threshold is proportional to the coverage.
   junc_label_type <- match.arg(junc_label_type);
   exon_label_type <- match.arg(exon_label_type);
   include_strand <- match.arg(include_strand);
   coord_method <- match.arg(coord_method);
   retVals <- list();
   #if (length(return_data) > 1 || "all" %in% return_data) {
   #}

   ## Validate filesDF
   if (length(filesDF) == 0 ||
         !jamba::igrepHas("tibble|tbl|data.*frame", class(filesDF))) {
      stop("filesDF must be a data.frame or equivalent");
   }
   if (!all(c("url","sample_id","type") %in% colnames(filesDF))) {
      stop("filesDF must contain colnames 'url', 'sample_id', and 'type'.");
   }
   ## validate other input
   if (!igrepHas("GRangesList", class(flatExonsByGene))) {
      stop("flatExonsByGene must be GRangesList")
   }
   if (length(sample_id) == 0) {
      sample_id <- unique(filesDF$sample_id);
   }

   ############################################
   ## Extract exons
   gr <- flatExonsByGene[gene]@unlistData;
   if (any(c("all") %in% return_data)) {
      retVals$gr <- gr;
   }

   ## Compress GRanges coordinates
   if (length(ref2c) == 0) {
      if (compressGR) {
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "running make_ref2compressed.");
         }
         ref2c <- make_ref2compressed(gr=gr,
            gapWidth=gapWidth,
            ...);
      } else {
         ref2c <- NULL;
         coord_method <- "none";
      }
   }
   if (any(c("all", "ref2c") %in% return_data)) {
      retVals$ref2c <- ref2c;
   }

   ############################################
   ## Get coverage data
   bwFilesDF <- filesDF[filesDF$type %in% "bw" &
      filesDF$sample_id %in% sample_id,,drop=FALSE];
   ## Subset by filename if include_strand is not "both"
   if (!"both" %in% include_strand) {
      if (igrepHas("pos|[+]|plus", bwFilesDF$url)) {
         if ("+" %in% include_strand) {
            bwFilesDF <- bwFilesDF[jamba::igrep("pos|[+]|plus", bwFilesDF$url),,drop=FALSE];
         } else {
            bwFilesDF <- bwFilesDF[jamba::unigrep("pos|[+]|plus", bwFilesDF$url),,drop=FALSE];
         }
      }
   }
   bwUrls <- nameVector(bwFilesDF[,c("url","url")]);
   bwSamples <- nameVector(bwFilesDF[,c("sample_id","url")]);
   bwUrlsL <- split(bwUrls, unname(bwSamples));
   if ("scale_factor" %in% colnames(bwFilesDF)) {
      bwScaleFactors <- rmNA(naValue=1,
         bwFilesDF$scale_factor);
   } else {
      bwScaleFactors <- rep(1, length(bwUrls));
   }
   names(bwScaleFactors) <- names(bwUrls);
   if (verbose) {
      if (any(bwScaleFactors != 1)) {
         printDebug("prepareSashimi(): ",
            "bwScaleFactors:",
            bwScaleFactors);
      }
      printDebug("prepareSashimi(): ",
         "bwUrls:");
      if (length(bwUrls) > 0) {
         print(data.frame(bwUrls));
      } else {
         printDebug("No coverage files", fgText="red");
      }
      printDebug("GRanges:");
      print(gr);
   }
   if (length(covGR) > 0) {
      ## Input is pre-processed coverage data
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "Preparing coverage from covGR");
      }
      if (do_shiny_progress) {
         ##
         shiny::incProgress(1/4,
            detail=paste0("Preparing GR coverage data for ", gene));
      }
      covGRuse <- covGR[names(covGR) %in% names(gr)];
      if (length(covGRuse) == 0) {
         warning("Supplied coverage covGR did not have names matching gr");
      } else {
         covDF <- filesDF[filesDF$type %in% "coverage_gr" &
               filesDF$url %in% colnames(values(covGR)) &
               filesDF$sample_id %in% sample_id,,drop=FALSE];
         if (nrow(covDF) == 0) {
            stop("Supplied coverage covGR does not have colnames in filesDF$url with filesDF$type == 'coverage_gr'");
         }
         covUrls <- nameVector(covDF[,c("url","url")]);
         covSamples <- nameVector(covDF[,c("sample_id","url")]);
         covGRuse <- covGRuse[,covUrls];
         if ("scale_factor" %in% colnames(covDF)) {
            covScaleFactors <- rmNA(naValue=1,
               covDF$scale_factor);
         } else {
            covScaleFactors <- rep(1, length(covUrls));
         }
         ## Combine coverage per strand
         covGR2 <- combineGRcoverage(covGRuse,
            covName=covSamples,
            scaleFactors=covScaleFactors,
            covNames=names(covUrls));
         ## Obtain the new set of covNames
         covNames <- attr(covGR2, "covNames");
         covName <- attr(covGR2, "covName");
         ## Create polygon data.frame
         covDF <- exoncov2polygon(covGR2,
            ref2c=ref2c,
            covNames=covNames,
            sample_id=covName,
            coord_style="fortify",
            ...);
         if (any(c("all", "covDF") %in% return_data)) {
            retVals$covDF <- covDF;
         }
      }
   }
   if (length(bwUrls) > 0) {
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "Preparing coverage from bigWig filesDF");
      }
      if (do_shiny_progress) {
         ##
         shiny::incProgress(1/4,
            detail=paste0("Preparing bw coverage data for ", gene));
      }
      covGR <- getGRcoverageFromBw(gr=gr,
         bwUrls=bwUrls,
         addGaps=addGaps,
         gap_feature_type=gap_feature_type,
         default_feature_type=default_feature_type,
         feature_type_colname=feature_type_colname,
         verbose=verbose,
         ...);
      ## Combine coverage per strand
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "Combining coverage by sample_id");
      }
      covGR2 <- combineGRcoverage(covGR,
         covName=bwSamples,
         scaleFactors=bwScaleFactors,
         covNames=names(bwUrls));
      #retVals$covGR <- covGR2;
      ## Obtain the new set of covNames
      covNames <- attr(covGR2, "covNames");
      covName <- attr(covGR2, "covName");
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covGR <- covGR;
         retVals$covGR2 <- covGR2;
      }

      ## Create polygon data.frame
      covDF <- exoncov2polygon(covGR2,
         ref2c=ref2c,
         covNames=covNames,
         sample_id=covName,
         coord_style="fortify",
         ...);
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covDF <- covDF;
      }
      ## Add feature_type data
      if (feature_type_colname %in% colnames(values(covGR2)) &&
            !feature_type_colname %in% colnames(covDF)) {
         covDF[,feature_type_colname] <- values(covGR2)[match(
            as.character(covDF$gr),
            names(covGR2)),feature_type_colname];
      }

      ########################################
      ## Optional exon labels
      covDFsub <- (as.character(covDF$gr) %in% names(gr));
      covDFlab <- covDF[covDFsub,,drop=FALSE];

      exonLabelDF1 <- shrinkMatrix(covDFlab[,c("x","y")],
         groupBy=pasteByRowOrdered(covDFlab[,c("gr", "sample_id")], sep=":!:"),
         shrinkFunc=function(x){mean(range(x))});
      exonLabelDF1[,c("gr","sample_id")] <- rbindList(
         strsplit(as.character(exonLabelDF1$groupBy), ":!:"));
      exonLabelDF1$gr <- factor(exonLabelDF1$gr, levels=unique(exonLabelDF1$gr));
      exonLabelDF <- renameColumn(exonLabelDF1,
         from="groupBy",
         to="gr_sample");
      retVals$exonLabelDF <- exonLabelDF;

   }

   ############################################
   ## Load Junctions
   ##
   ## Pre-existing junctions supplied as GRanges
   if (length(juncGR) > 0) {
      if (!"sample_id" %in% colnames(values(juncGR))) {
         stop("juncGR must have 'sample_id' in colnames(values(juncGR)).");
      }
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "Preparing junctions from juncGR");
      }
      if (do_shiny_progress) {
         ##
         shiny::incProgress(2/4,
            detail=paste0("Preparing GR junction data for ", gene));
      }
      juncGRuse <- juncGR[values(juncGR)[["sample_id"]] %in% sample_id];
      if (length(juncGRuse) == 0) {
         warning("Supplied junctions juncGR did not have names matching sample_id");
      } else {
         ## Create junction summary data.frame
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "running spliceGR2junctionDF for juncGR()");
         }
         juncDF1 <- spliceGR2junctionDF(spliceGRgene=juncGR,
            exonsGR=gr,
            sampleColname="sample_id");
         ## Subset junctions by minimum score
         if (length(minJunctionScore) > 0 && minJunctionScore > 0) {
            juncDF1 <- subset(juncDF1, abs(score) >= minJunctionScore);
         }
      }
   } else {
      juncDF1 <- NULL;
   }
   ##
   ## Junctions available from filesDF type %in% "junction"
   ##
   juncFilesDF <- filesDF[filesDF$type %in% "junction" &
      filesDF$sample_id %in% sample_id,c("url", "sample_id", "scale_factor"),drop=FALSE];
   juncUrls <- nameVector(juncFilesDF[,c("url", "sample_id")]);
   juncSamples <- nameVector(juncFilesDF[,c("sample_id","sample_id")]);
   juncUrlsL <- split(juncUrls, juncSamples);
   if (!"scale_factor" %in% colnames(juncFilesDF)) {
      juncFilesDF$scale_factor <- 1;
   }
   juncScaleFactors <- nameVector(juncFilesDF[,c("scale_factor","sample_id")]);
   if (verbose) {
      printDebug("prepareSashimi(): ",
         "Applying scale_factor values to junction scores,",
         "juncScaleFactors:",
         format(juncScaleFactors, digits=2));
   }
   if (verbose) {
      printDebug("prepareSashimi(): ",
         "juncUrls:");
      if (length(juncUrls) > 0) {
         print(data.frame(juncUrls));
      } else {
         printDebug("No junction files", fgText="red");
      }
   }
   if (length(juncUrls) > 0) {
      if (do_shiny_progress) {
         ##
         shiny::incProgress(2/4,
            detail=paste0("Preparing BED junction data for ", gene));
      }
      juncBedGR <- GRangesList(lapply(nameVectorN(juncUrls), function(iBedName){
         iBed <- juncUrls[[iBedName]];
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "Importing bed:",
               iBed,
               " with scale_factor:",
               format(digits=1, juncScaleFactors[iBedName]),
               " for sample_id:", juncSamples[iBedName]);
         }
         bed1 <- rtracklayer::import(iBed,
            which=range(gr));
         ## Assign score, apply scale_factor
         ## Consider rounding the score to integer value?
         values(bed1)$score <- as.numeric(as.character(values(bed1)$name)) * juncScaleFactors[iBedName];

         values(bed1)[,c("juncNames")] <- iBedName;
         values(bed1)[,c("sample_id")] <- juncSamples[iBedName];

         ## Subset junctions to require one end within the region of interest
         bed1 <- subset(bed1,
            (
               overlapsAny(flank(bed1, -1, start=TRUE), range(gr)) |
               overlapsAny(flank(bed1, -1, start=FALSE), range(gr))
            )
         );
         bed1;
      }))@unlistData;
      ## Create junction summary data.frame
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "running spliceGR2junctionDF for juncBedGR()");
      }
      juncDF1f <- spliceGR2junctionDF(spliceGRgene=juncBedGR,
         exonsGR=gr,
         sampleColname="sample_id");
      if (length(juncDF1) > 0) {
         if (verbose) {
            printDebug("prepareSashimi(): ",
               "Appending BED-derived junctions to supplied juncGR-derived data.");
         }
         juncDF1 <- rbind(juncDF1, juncDF1f);
      } else {
         juncDF1 <- juncDF1f;
      }
      ## Subset junctions by minimum score
      if (length(minJunctionScore) > 0 && minJunctionScore > 0) {
         juncDF1 <- subset(juncDF1, abs(score) >= minJunctionScore);
      }
   }
   if (exists("juncDF1") && length(juncDF1) > 0 && nrow(juncDF1) > 0) {
      juncGR <- as(renameColumn(juncDF1, from="ref", to="seqnames"), "GRanges");
      names(juncGR) <- makeNames(values(juncGR)[,"nameFromToSample"]);
      ## Convert junctions to polygons usable by geom_diagonal_wide()
      #      juncPolyDF <- grl2df(setNames(GRangesList(subset(juncGR, score > minJunctionScore)), sample_id),
      if (length(baseline) == 0) {
         baseline <- 0;
      }
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "calling grl2df() on juncGR");
         print(head(juncGR));
         print(head(GenomicRanges::split(juncGR, values(juncGR)[["sample_id"]])));
         printDebug("prepareSashimi(): ",
            "calling grl2df() on juncGR");
      }
      juncDF <- grl2df(
         GenomicRanges::split(juncGR, values(juncGR)[["sample_id"]]),
         shape="junction",
         ref2c=ref2c,
         scoreFactor=1,
         scoreArcFactor=scoreArcFactor,
         scoreArcMinimum=scoreArcMinimum,
         baseline=baseline,
         doStackJunctions=doStackJunctions,
         verbose=verbose,
         ...);
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "called grl2df() on juncGR");
      }
      if (!"sample_id" %in% colnames(juncDF)) {
         juncDF <- renameColumn(juncDF,
            from="grl_name",
            to="sample_id");
      }
      if (verbose) {
         printDebug("prepareSashimi(): ",
            "dim(juncDF):", dim(juncDF));
      }
      juncDF <- juncDF[,!colnames(juncDF) %in% c("grl_name"),drop=FALSE];

      if (any(c("all", "juncDF") %in% return_data)) {
         retVals$juncGR <- juncGR;
         retVals$juncDF1 <- juncDF1;
         retVals$juncDF <- juncDF;
      }

      ## define junction label positions
      #juncLabelDF1 <- subset(mutate(juncCoordDF, id_name=makeNames(id)), grepl("_v1_v3$", id_name));
      if (do_shiny_progress) {
         ##
         shiny::incProgress(3/4,
            detail=paste0("Preparing junction labels for ", gene));
      }
      juncLabelDF1 <- subset(plyr::mutate(juncDF, id_name=makeNames(id)),
         grepl("_v1_v[23]$", id_name));
      shrink_colnames <- intersect(c("x","y","score","junction_rank"),
         colnames(juncLabelDF1));
      juncLabelDF <- renameColumn(
         shrinkMatrix(juncLabelDF1[,shrink_colnames,drop=FALSE],
            groupBy=juncLabelDF1[,"nameFromToSample"]),
         from="groupBy",
         to="nameFromToSample");
      juncLabelDF[,c("nameFromTo","sample_id")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromToSample"], ":!:"));
      juncLabelDF[,c("nameFrom", "nameTo")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromTo"], " "));

      if (any(c("all", "juncLabelDF") %in% return_data)) {
         retVals$juncLabelDF <- juncLabelDF;
      }
   } else {
      juncDF1 <- NULL;
      juncDF <- NULL;
      juncLabelDF <- NULL;
   }
   if (do_shiny_progress) {
      ##
      shiny::incProgress(4/4,
         detail=paste0("Sashimi data is ready for ", gene));
   }

   return(retVals);
}

#' maximum overlapping internal junction score
#'
#' maximum overlapping internal junction score
#'
#' This function is intended for internal use, and calculates the
#' maximum score for junction GRanges for the special case where:
#'
#' * a junction range overlaps the start or end of another junction
#' * the start or end of the other junction is internal, that is
#' it does not overlap the same start or end.
#' * overlaps are constrained to the same `"sample_id"` stored
#' in `values(juncGR)$sample_id`.
#'
#' The goal is to determine the highest score contained
#' inside each junction region, so that the junction arc height
#' can be defined higher in order to minimize junction arc
#' overlaps.
#'
#' @return numeric vector named by `names(juncGR)` whose values
#'    are the maximum score of internal overlapping junction ends.
#'
#' @param juncGR GRanges containing splice junctions, with numeric
#'    column `scoreColname` representing the abundance of splice
#'    junction-spanning sequence reads.
#' @param scoreColname,sampleColname colnames in `values(juncGR)`
#'    which define the junction score, and the `"sample_id"`.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
internal_junc_score <- function
(juncGR,
 scoreColname="score",
 sampleColname="sample_id",
 verbose=FALSE,
   ...)
{
   ## Purpose is to sum scores for junctions which overlap
   ## junction ends inside the junction boundary.
   if (verbose) {
      printDebug("internal_junc_score(): ",
         "Defining flank ends of each junction.");
   }
   juncEndsGR <- c(
      flank(juncGR[,c(scoreColname,sampleColname)],
         start=TRUE,
         width=1),
      flank(juncGR[,c(scoreColname,sampleColname)],
         start=FALSE,
         width=1));
   values(juncEndsGR)$side <- rep(c("start", "end"), each=length(juncGR));
   values(juncEndsGR)$id <- rep(seq_along(juncGR), 2);
   juncEndsDF <- as.data.frame(unname(juncEndsGR));
   juncGroup <- pasteByRow(juncEndsDF[,c(sampleColname, "seqnames", "start", "strand", "side")]);
   juncEndsRed <- shrinkMatrix(juncEndsDF[[scoreColname]],
      groupBy=juncGroup,
      shrinkFunc=sum);
   if (verbose) {
      printDebug("internal_junc_score(): ",
         "dim(juncEndsRed):", dim(juncEndsRed));
   }
   juncEndsRefGR <- juncEndsGR[match(juncEndsRed$groupBy, juncGroup)];
   values(juncEndsRefGR)[[scoreColname]] <- juncEndsRed$x;
   if (verbose) {
      printDebug("internal_junc_score(): ",
         "length(juncEndsRefGR):", length(juncEndsRefGR));
   }

   fo1 <- findOverlaps(juncGR, juncEndsRefGR);
   fo1df <- data.frame(q=queryHits(fo1),
      s=subjectHits(fo1));
   if (verbose) {
      printDebug("internal_junc_score(): ",
         "dim(fo1df):", dim(fo1df));
   }
   if (nrow(fo1df) == 0) {
      intScore <- nameVector(rep(0, length(juncGR)),
         names(juncGR));
   } else {
      fo1df$qSample <- values(juncGR[fo1df$q])[[sampleColname]];
      fo1df$sSample <- values(juncEndsRefGR[fo1df$s])[[sampleColname]];
      fo1df$sScore <- values(juncEndsRefGR[fo1df$s])[[scoreColname]];

      fo1dfuse <- subset(fo1df,
         qSample == sSample);

      intScore <- rmNA(naValue=0,
         max(List(split(fo1dfuse$sScore, fo1dfuse$q)))[as.character(seq_along(juncGR))]);
      names(intScore) <- names(juncGR);
   }
   return(intScore);
}
