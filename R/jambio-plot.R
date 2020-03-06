




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
         round(median(GenomicRanges::width(GenomicRanges::reduce(gr))) / 2)));
      if (verbose) {
         jamba::printDebug("make_ref2compressed(): ",
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
   disGR <- GenomicRanges::reduce(gr);
   if (is.null(names(disGR))) {
      names(disGR) <- jamba::makeNames(rep("disGR", length(disGR)));
   }
   if (length(disGR) > 1) {
      exonGaps <- sapply(head(seq_along(disGR), -1), function(i1){
         i2 <- i1 + 1;
         GenomicRanges::distance(disGR[i1], disGR[i2]) > 0;
      });
      names(exonGaps) <- head(names(disGR), -1);
   } else {
      exonGaps <- NULL;
   }

   ## Now reset coordinates using fixed gap width
   #newWidths <- head(intercalate(width(disGR),
   #   (gapWidth)*exonGaps), -1);
   newWidths <- suppressWarnings(
      head(intercalate(GenomicRanges::width(disGR),
         (gapWidth)*exonGaps),
         length(disGR)*2-1)
      );

   ##
   if (verbose) {
      jamba::printDebug("make_ref2compressed(): ",
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
      jamba::printDebug("make_ref2compressed(): ",
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
      jamba::printDebug("make_ref2compressed(): ",
         "ref2compressedGR");
   }
   ref2compressedGR <- function(gr) {
      if (length(names(gr)) > 0) {
         GRnames <- names(gr);
      }
      if (all(c("refStart","refEnd") %in% colnames(GenomicRanges::values(GR)))) {
         ## Re-compress the original reference coordinates
         ranges(gr) <- IRanges::IRanges(
            start=GenomicRanges::values(gr)[,"refStart"],
            end=GenomicRanges::values(gr)[,"refEnd"]
         );
      } else {
         GenomicRanges::values(gr)[,"refStart"] <- GenomicRanges::start(gr);
         GenomicRanges::values(gr)[,"refEnd"] <- GenomicRanges::end(gr);
      }
      ranges(gr) <- IRanges::IRanges(
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
         names(grl) <- jamba::makeNames(rep("grl", length(grl)));
      }
      grlc1 <- ref2compressedGR(grl@unlistData);
      grlc <- GenomicRanges::split(grlc1,
         factor(
            rep(names(grl), S4Vectors::elementNROWS(grl)),
            levels=unique(names(grl))));
      grlc;
   }
   if (verbose) {
      jamba::printDebug("make_ref2compressed(): ",
         "retVals");
   }

   ## Custom breaks function
   breaks_gr <- function
   (x,
    limits=NULL,
    n=nBreaks,
    xFixed=lookupCoordDF[,1],
    verbose=FALSE,
    ...)
   {
      xvals <- unique(sort(xFixed));
      xvals <- xvals[xvals >= min(x) & xvals <= max(x)];
      if (verbose) {
         jamba::printDebug("breaks_gr(): ",
            "x:", x);
         jamba::printDebug("breaks_gr(): ",
            "limits:", limits);
         jamba::printDebug("breaks_gr(): ",
            "n:", n);
      }
      if (n > length(xvals)) {
         return(xvals);
      }
      idx1 <- round(seq.int(from=1, to=length(xvals), length.out=n));
      return(xvals[idx1]);
   }
   minor_breaks <- function
   (b,
    limits=NULL,
    n=2,
    xFixed=lookupCoordDF[,1],
    verbose=FALSE,
    ...)
   {
      if (verbose) {
         jamba::printDebug("minor_breaks(): ",
            "x:", x);
         jamba::printDebug("minor_breaks(): ",
            "limits:", limits);
         jamba::printDebug("minor_breaks(): ",
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
 minRatio=5,
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
      jamba::printDebug("compressPolygonM(): ",
         "data_style:", data_style);
      jamba::printDebug("compressPolygonM(): ",
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
   sdim_polyDFL <- sdim(polyDFL)$rows;
   whichComp <- which(polyMratios >= minRatio & sdim_polyDFL > 5);
   whichNorm <- which(polyMratios < minRatio | sdim_polyDFL <= 5);
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
         jamba::printDebug("head(iDFu):");print(head(iDFu));
         jamba::printDebug("head(iYnew):", head(iYnew));
         jamba::printDebug("head(iMultiple):", head(iMultiple));
         jamba::printDebug("head(iN):", head(iN));
         jamba::printDebug("head(iSeqNew):", head(iSeqNew));
         jamba::printDebug("head(iYrunmax):", head(iYrunmax));
         jamba::printDebug("head(iSeqSub):", head(iSeqSub));
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

   if (!igrepHas("granges", class(gr))) {
      stop("Input gr should be a GRanges object.");
   }
   if (length(covNames) == 0) {
      covNames <- provigrep(c("pos|[+]$", "neg|[-]$"),
         colnames(GenomicRanges::values(gr)));
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
      jamba::printDebug("exoncov2polygon(): ",
         "covNames:",
         covNames);
   }
   covPolyL <- lapply(nameVector(covNames), function(iName){
      polyL <- lapply(nameVectorN(gr), function(iGRname){
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
         if (length(xVals) != length(yVals)) {
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
      if (length(ref2c) > 0 && compress_introns) {
         t1 <- Sys.time();
         compCovPolyML <- lapply(covPolyML, function(iL){
            iML <- compressPolygonM(iL,
               ref2c=ref2c,
               coord_style=coord_style,
               verbose=verbose);
         });
         t2 <- Sys.time();
         if (verbose) {
            jamba::printDebug("exoncov2polygon(): ",
               "Completed polygon compression:",
               format(t2 - t1));
         }
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
         if (exists("compCovPolyML")) {
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
#'    `jamba::makeNames(basename(bwUrls))`. Each column is type `NumericList`,
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
   if (!igrepHas("GRanges", class(gr))) {
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
         warning(warnText);
         NULL;
      });
      cov1;
   }
   if (use_memoise) {
      import_or_null_m <- memoise::memoise(import_or_null,
         cache=memoise::cache_filesystem(memoise_coverage_path));
   }
   ## Iterate bwUrls and get coverage from each
   if (verbose) {
      jamba::printDebug("getGRcoverageFromBw(): ",
         "bwUrls:");
      print(bwUrls);
   }
   covL <- lapply(nameVectorN(bwUrls), function(iBw){
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
#' @param gr GRanges object containing coverage data in columns
#'    containing NumericList class data.
#' @param covNames character vector of `colnames(values(gr))`
#'    representing columns in `gr` that contain coverage data
#'    in NumericList format, for example data prepared
#'    with `getGRcoverageFromBw()`.
#' @param covName character vector with length equal to
#'    `length(covNames)` representing the `sample_id` for each
#'    `covNames` entry.
#' @param strands character vector, or NULL, indicating the strand
#'    for which the coverage data was obtained. When `NULL` the strand
#'    is inferred by the presence of any negative values.
#' @param scaleFactors numeric vector length equal to `length(covNames)`
#'    or expanded to that length. Values are multiplied by each
#'    coverage data result, intended to apply a normalization
#'    to each coverage value. A `-1` value can also be used
#'    to flip the score of negative strand data, in the event the
#'    source coverage data is scored only using positive values.
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
   #   stop("combineGRcoverage() found no covNames present in colnames(values(gr)).");
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
            ifelse(any(any(GenomicRanges::values(gr)[[iCov]] * scaleFactors[i] < 0) &
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
   gr <- gr[,keep_gr_colnames];
   attr(gr, "covNames") <- names(covNamesL);
   attr(gr, "covName") <- jamba::cPasteU(covNameL);
   attr(gr, "some_null") <- some_null;
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
#' @param compress_introns logical indicating whether to compress
#'    the coverage polygon coordinates to approximately the same
#'    number of pixels per inch as the exon polygons. This option
#'    greatly reduces the size of the polygon, since introns are
#'    already about 50 to 100 times wider than exons, and when
#'    `compressGR` is `TRUE`, the introns are visibly compressed
#'    to a fixed width on the x-axis. The data has many more
#'    x-axis coordinates than the data visualization, this argument
#'    is intended to reduce the intron coordinates accordingly.
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
#'    `shiny::withProgress()` and `shiny::setProgress()` methods.
#'    This function only calls `shiny::setProgress()` and
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
 compress_introns=TRUE,
 ref2c=NULL,
 gap_feature_type="intron",
 default_feature_type="exon",
 feature_type_colname="feature_type",
 exon_label_type=c("none", "repel", "mark"),
 junc_label_type=c("repel", "mark", "none"),
 return_data=c("df", "ref2c"),
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
 use_memoise=FALSE,
 memoise_coverage_path="coverage_memoise",
 memoise_junction_path="junctions_memoise",
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

   ## Validate filesDF
   if (length(filesDF) == 0 ||
         !jamba::igrepHas("tibble|tbl|data.*frame", class(filesDF))) {
      stop("filesDF must be a data.frame or equivalent");
   }
   if (!all(c("url","sample_id","type") %in% colnames(filesDF))) {
      stop("filesDF must contain colnames 'url', 'sample_id', and 'type'.");
   }
   if (!"scale_factor" %in% colnames(filesDF)) {
      filesDF[,"scale_factor"] <- rep(1, nrow(filesDF));
   }
   ## validate other input
   if (!igrepHas("GRangesList", class(flatExonsByGene))) {
      stop("flatExonsByGene must be GRangesList")
   }
   if (length(sample_id) == 0) {
      sample_id <- unique(as.character(filesDF$sample_id));
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
            jamba::printDebug("prepareSashimi(): ",
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
   retVals$some_null <- FALSE;

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
   if (!"scale_factor" %in% colnames(bwFilesDF)) {
      bwFilesDF$scale_factor <- rep(1, nrow(bwFilesDF));
   }
   bwScaleFactors <- rmNA(naValue=1,
      bwFilesDF$scale_factor);
   names(bwScaleFactors) <- names(bwUrls);
   if (length(bwScaleFactors) > 0) {
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "bwScaleFactors:",
            paste0(basename(names(bwScaleFactors)),
               ":",
               format(bwScaleFactors, digits=2, trim=TRUE)),
            sep=", ");
      }
   }
   if (verbose) {
      jamba::printDebug("prepareSashimi(): ",
         "bwUrls:");
      if (length(bwUrls) > 0) {
         print(data.frame(bwUrls));
      } else {
         jamba::printDebug("No coverage files", fgText="red");
      }
      #printDebug("GRanges:");
      #print(gr);
   }
   ##
   ## Where possible, re-use covGR coverage supplied as GRanges
   ##
   if (length(covGR) > 0) {
      ## Input is pre-processed coverage data
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Preparing coverage from covGR");
      }
      if (do_shiny_progress) {
         ##
         shiny::setProgress(0/4,
            detail=paste0("Preparing GR coverage data for ", gene));
      }
      covGRuse <- covGR[names(covGR) %in% names(gr)];
      if (length(covGRuse) == 0) {
         warning("Supplied coverage covGR did not have names matching gr");
      } else {
         covfilesDF <- filesDF[filesDF$type %in% "coverage_gr" &
               filesDF$url %in% colnames(GenomicRanges::values(covGR)) &
               filesDF$sample_id %in% sample_id,,drop=FALSE];
         if (nrow(covfilesDF) == 0) {
            stop("Supplied coverage covGR does not have colnames in filesDF$url with filesDF$type == 'coverage_gr'");
         }
         covUrls <- nameVector(covfilesDF[,c("url","url")]);
         covSamples <- nameVector(covfilesDF[,c("sample_id","url")]);
         covGRuse <- covGRuse[,covUrls];
         if ("scale_factor" %in% colnames(covfilesDF)) {
            covScaleFactors <- rmNA(naValue=1,
               covfilesDF$scale_factor);
         } else {
            covScaleFactors <- rep(1, length(covUrls));
         }

         ## Combine coverage per strand
         covGR2 <- combineGRcoverage(covGRuse,
            covName=covSamples,
            scaleFactors=covScaleFactors,
            covNames=names(covUrls),
            verbose=verbose);
         ## Obtain the new set of covNames
         covNames <- attr(covGR2, "covNames");
         covName <- attr(covGR2, "covName");
         ## some_null is TRUE when any underlying coverage file failed to return coverage
         some_null <- attr(covGR2, "some_null");
         if (length(some_null) && some_null) {
            retVals$some_null <- some_null;
         }

         ## Create polygon data.frame
         covDF <- exoncov2polygon(covGR2,
            ref2c=ref2c,
            covNames=covNames,
            sample_id=covName,
            coord_style="fortify",
            compress_introns=compress_introns,
            verbose=verbose,
            ...);
         ## Enforce ordered factor levels for sample_id
         ## which also forces empty factor levels if applicable
         covDF$sample_id <- factor(
            as.character(covDF$sample_id),
            levels=unique(sample_id)
         );
         if (any(c("all", "covDF") %in% return_data)) {
            retVals$covDF <- covDF;
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
         exonLabelDF1$gr <- factor(exonLabelDF1$gr,
            levels=unique(exonLabelDF1$gr));
         exonLabelDF <- renameColumn(exonLabelDF1,
            from="groupBy",
            to="gr_sample");
         if (any(c("all", "covDF") %in% return_data)) {
            retVals$exonLabelDF <- exonLabelDF;
         }
      }
   }
   ##
   ## load coverage from bigWig files defined in filesDF
   ##
   if (length(bwUrls) > 0) {
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Preparing coverage from bigWig filesDF");
      }
      if (do_shiny_progress) {
         ##
         shiny::setProgress(1/4,
            detail=paste0("Preparing bw coverage data for ", gene));
      }
      ## Note that coverage is not scaled at this step
      covGR <- getGRcoverageFromBw(gr=gr,
         bwUrls=bwUrls,
         addGaps=addGaps,
         gap_feature_type=gap_feature_type,
         default_feature_type=default_feature_type,
         feature_type_colname=feature_type_colname,
         use_memoise=use_memoise,
         memoise_coverage_path=memoise_coverage_path,
         do_shiny_progress=do_shiny_progress,
         verbose=verbose,
         ...);
      ## Combine coverage per strand
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Combining coverage by sample_id");
      }

      covGR2 <- combineGRcoverage(covGR,
         covName=bwSamples,
         scaleFactors=bwScaleFactors,
         covNames=names(bwUrls),
         verbose=verbose);
      #retVals$covGR <- covGR2;
      ## Obtain the new set of covNames
      covNames <- attr(covGR2, "covNames");
      covName <- attr(covGR2, "covName");
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covGR <- covGR;
         retVals$covGR2 <- covGR2;
      }
      ## some_null is TRUE when any underlying coverage file failed to return coverage
      some_null <- attr(covGR2, "some_null");
      if (length(some_null) && some_null) {
         retVals$some_null <- some_null;
      }

      ## Create polygon data.frame
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Calling exoncov2polygon()");
      }
      covDF <- exoncov2polygon(covGR2,
         ref2c=ref2c,
         covNames=covNames,
         sample_id=covName,
         coord_style="fortify",
         compress_introns=compress_introns,
         verbose=verbose,
         ...);
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$covDF <- covDF;
      }
      ## Add feature_type data
      if (feature_type_colname %in% colnames(GenomicRanges::values(covGR2)) &&
            !feature_type_colname %in% colnames(covDF)) {
         covDF[,feature_type_colname] <- GenomicRanges::values(covGR2)[match(
            as.character(covDF$gr),
            names(covGR2)),feature_type_colname];
      }
      ## Enforce ordered factor levels for sample_id
      ## which also forces empty factor levels if applicable
      covDF$sample_id <- factor(
         as.character(covDF$sample_id),
         levels=unique(sample_id)
      );
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "head(covDF$sample_id):");
         print(head(covDF$sample_id));
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
      exonLabelDF1$gr <- factor(exonLabelDF1$gr,
         levels=unique(exonLabelDF1$gr));
      exonLabelDF <- renameColumn(exonLabelDF1,
         from="groupBy",
         to="gr_sample");
      if (any(c("all", "covDF") %in% return_data)) {
         retVals$exonLabelDF <- exonLabelDF;
      }

   }

   ############################################
   ## Load Junctions
   ##
   ## Pre-existing junctions supplied as GRanges
   if (length(juncGR) > 0) {
      if (!"sample_id" %in% colnames(GenomicRanges::values(juncGR))) {
         stop("juncGR must have 'sample_id' in colnames(values(juncGR)).");
      }
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
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
            jamba::printDebug("prepareSashimi(): ",
               "running spliceGR2junctionDF for juncGR()");
         }
         juncDF1 <- spliceGR2junctionDF(spliceGRgene=juncGR,
            exonsGR=gr,
            sampleColname="sample_id");
         ## Subset junctions by minimum score
         if (length(minJunctionScore) > 0 &&
               minJunctionScore > 0 &&
               length(juncDF1) > 0) {
            juncDF1 <- subset(juncDF1, abs(score) >= minJunctionScore);
         }
         if (length(juncDF1) == 0 || nrow(juncDF1) == 0) {
            juncDF1 <- NULL;
         } else {
            ## Enforce ordered factor levels for sample_id
            ## which also forces empty factor levels if applicable
            juncDF1$sample_id <- factor(
               as.character(juncDF1$sample_id),
               levels=unique(sample_id)
            );
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
   juncSamples <- nameVector(juncFilesDF[,c("sample_id", "sample_id")]);
   juncUrlsL <- split(juncUrls, juncSamples);
   if (!"scale_factor" %in% colnames(juncFilesDF)) {
      juncFilesDF$scale_factor <- rep(1, nrow(juncFilesDF));
   }
   juncScaleFactors <- rmNA(naValue=1,
      nameVector(juncFilesDF[,c("scale_factor", "sample_id")]));
   if (verbose && length(juncScaleFactors) > 0) {
      jamba::printDebug("prepareSashimi(): ",
         "juncScaleFactors: ",
         paste0(names(juncScaleFactors),
            ":",
            format(juncScaleFactors, digits=2, trim=TRUE)),
         sep=", ");
   }
   if (verbose) {
      jamba::printDebug("prepareSashimi(): ",
         "juncUrls:");
      if (length(juncUrls) > 0) {
         print(data.frame(juncUrls));
      } else {
         jamba::printDebug("No junction files", fgText="red");
      }
   }
   if (length(juncUrls) > 0) {
      if (use_memoise) {
         import_juncs_m <- memoise::memoise(import_juncs_from_bed,
            cache=memoise::cache_filesystem(memoise_junction_path));
      }
      if (do_shiny_progress) {
         shiny::incProgress(2/4,
            detail=paste0("Importing junctions for ", gene));
      }
      juncBedList <- lapply(nameVectorN(juncUrls), function(iBedName){
         iBed <- juncUrls[[iBedName]];
         iBedNum <- match(iBedName, jamba::makeNames(names(juncUrls)));
         iBedPct <- (iBedNum - 1) / length(juncUrls);
         if (do_shiny_progress) {
            if (!is.na(iBedPct)) {
               shiny::setProgress(
                  value=2/4 + iBedPct/4,
                  detail=paste0("Importing junctions (",
                     iBedNum,
                     " of ",
                     length(juncUrls),
                     ") for ", gene));
            }
         }
         if (verbose) {
            jamba::printDebug("prepareSashimi(): ",
               "progress:", format(digits=2, 2/4 + iBedPct/4),
               "Importing bed:",
               iBed,
               " with scale_factor:",
               format(digits=1, juncScaleFactors[iBedName]),
               " for sample_id:", juncSamples[iBedName]);
         }
         ## Question: should scale_factor be applied here within the cache,
         ## or should cache just store the data without scale_factor, so
         ## the scale_factor can be applied independently?
         ## I think we know the answer is to cache the data, apply
         ## scale_factor separately, to allow the scale_factor to be
         ## changed as needed.
         if (use_memoise) {
            if (verbose) {
               import_juncs_m_cached <- memoise::has_cache(import_juncs_m)(
                  iBed,
                  juncNames=iBedName,
                  sample_id=juncSamples[iBedName],
                  scale_factor=1,
                  #scale_factor=juncScaleFactors[iBedName],
                  gr=gr);
               jamba::printDebug("prepareSashimi():",
                  "import_juncs_m_cached: ",
                  import_juncs_m_cached);
            }
            bed1 <- import_juncs_m(
               iBed,
               juncNames=iBedName,
               sample_id=juncSamples[iBedName],
               scale_factor=1,
               #scale_factor=juncScaleFactors[iBedName],
               use_memoise=TRUE,
               memoise_junction_path=memoise_junction_path,
               gr=gr);
         } else {
            bed1 <- import_juncs_from_bed(iBed,
               juncNames=iBedName,
               sample_id=juncSamples[iBedName],
               scale_factor=1,
               #scale_factor=juncScaleFactors[iBedName],
               gr=gr);
         }
         ## Apply scale_factor
         if (length(bed1) > 0) {
            values(bed1)$score <- values(bed1)$score * juncScaleFactors[iBedName];
         }
         bed1;
      });
      juncBedList <- juncBedList[lengths(juncBedList) > 0];
      if (length(juncBedList) == 0) {
         juncDF1 <- NULL;
      } else {
         juncBedGR <- GRangesList(juncBedList)@unlistData;
         ## Create junction summary data.frame
         if (verbose) {
            jamba::printDebug("prepareSashimi(): ",
               "running spliceGR2junctionDF for juncBedGR()");
         }
         if (do_shiny_progress) {
            shiny::setProgress(3/4,
               detail=paste0("Combining junction data for ", gene));
         }
         juncDF1f <- spliceGR2junctionDF(spliceGRgene=juncBedGR,
            exonsGR=gr,
            sampleColname="sample_id");
         if (length(juncDF1) > 0) {
            if (verbose) {
               jamba::printDebug("prepareSashimi(): ",
                  "Appending BED-derived junctions to supplied juncGR-derived data.");
            }
            juncDF1 <- rbind(juncDF1, juncDF1f);
         } else {
            juncDF1 <- juncDF1f;
         }
         if (length(juncDF1) == 0 || nrow(juncDF1) == 0) {
            juncDF1 <- NULL;
         }
         ## Subset junctions by minimum score
         if (length(minJunctionScore) > 0 &&
               minJunctionScore > 0 &&
               length(juncDF1) > 0) {
            juncDF1 <- subset(juncDF1, abs(score) >= minJunctionScore);
         }
         if (length(juncDF1) == 0 || nrow(juncDF1) == 0) {
            juncDF1 <- NULL;
         } else {
            ## Enforce ordered factor levels for sample_id
            ## which also forces empty factor levels if applicable
            juncDF1$sample_id <- factor(
               as.character(juncDF1$sample_id),
               levels=unique(sample_id)
            );
         }
      }
   }
   if (exists("juncDF1") &&
         length(juncDF1) > 0 &&
         length(dim(juncDF1)) > 1) {
      juncGR <- as(
         renameColumn(juncDF1,
            from="ref",
            to="seqnames"),
         "GRanges");
      names(juncGR) <- jamba::makeNames(GenomicRanges::values(juncGR)[,"nameFromToSample"]);

      if (length(baseline) == 0) {
         baseline <- 0;
      }
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "calling grl2df() on juncGR");
         #print(head(juncGR));
         #print(head(GenomicRanges::split(juncGR,
         #   GenomicRanges::values(juncGR)[["sample_id"]])));
      }
      if (do_shiny_progress) {
         shiny::setProgress(3.5/4,
            detail=paste0("Calculating junction stacking for ", gene));
      }
      ## Convert junctions to polygons usable by geom_diagonal_wide()
      juncDF <- grl2df(
         GenomicRanges::split(juncGR,
            GenomicRanges::values(juncGR)[["sample_id"]]),
         shape="junction",
         ref2c=ref2c,
         scoreFactor=1,
         scoreArcFactor=scoreArcFactor,
         scoreArcMinimum=scoreArcMinimum,
         baseline=baseline,
         doStackJunctions=doStackJunctions,
         verbose=verbose,
         ...);
      if (!"sample_id" %in% colnames(juncDF)) {
         juncDF <- renameColumn(juncDF,
            from="grl_name",
            to="sample_id");
      }
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "dim(juncDF):", dim(juncDF));
      }
      juncDF <- juncDF[,!colnames(juncDF) %in% c("grl_name"),drop=FALSE];
      ## Enforce ordered factor levels for sample_id
      ## which also forces empty factor levels if applicable
      juncDF$sample_id <- factor(
         as.character(juncDF$sample_id),
         levels=unique(sample_id)
      );

      if (any(c("all", "juncDF") %in% return_data)) {
         retVals$juncGR <- juncGR;
         retVals$juncDF1 <- juncDF1;
         retVals$juncDF <- juncDF;
      }

      ## define junction label positions
      #juncLabelDF1 <- subset(mutate(juncCoordDF, id_name=makeNames(id)), grepl("_v1_v3$", id_name));
      if (do_shiny_progress) {
         ##
         shiny::setProgress(3.8/4,
            detail=paste0("Preparing junction label coordinates for ", gene));
      }
      juncLabelDF1 <- subset(plyr::mutate(juncDF, id_name=jamba::makeNames(id)),
         grepl("_v1_v[23]$", id_name));
      shrink_colnames <- intersect(c("x","y","score","junction_rank"),
         colnames(juncLabelDF1));
      ## Define junction placement at max position, in stranded fashion
      juncLabelDF_y <- renameColumn(
         shrinkMatrix(juncLabelDF1[,"y",drop=FALSE],
            shrinkFunc=function(x){max(abs(x))*sign(max(x))},
            groupBy=juncLabelDF1[,"nameFromToSample"]),
         from="groupBy",
         to="nameFromToSample");
      juncLabelDF <- renameColumn(
         shrinkMatrix(juncLabelDF1[,shrink_colnames,drop=FALSE],
            groupBy=juncLabelDF1[,"nameFromToSample"]),
         from="groupBy",
         to="nameFromToSample");
      juncLabelDF$y <- juncLabelDF_y$y;
      juncLabelDF[,c("nameFromTo","sample_id")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromToSample"], ":!:"));
      juncLabelDF[,c("nameFrom", "nameTo")] <- rbindList(
         strsplit(juncLabelDF[,"nameFromTo"], " "));
      ## Enforce ordered factor levels for sample_id
      ## which also forces empty factor levels if applicable
      juncLabelDF$sample_id <- factor(
         as.character(juncLabelDF$sample_id),
         levels=unique(sample_id)
      );
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
      shiny::setProgress(4/4,
         detail=paste0("Sashimi data is ready for ", gene));
   }

   ## Merge data.frame entries together
   ## exon coverage
   if (exists("covDF") && length(covDF) > 0) {
      covDF$type <- "coverage";
      covDF$name <- pasteByRow(covDF[,c("gr", "cov", "sample_id")], sep=" ");
      covDF$feature <- covDF$gr;
      covDF$row <- seq_len(nrow(covDF));
      ## define name as factor to maintain the drawing order
      covDF$name <- factor(covDF$name,
         levels=unique(covDF$name));
   } else {
      covDF <- NULL;
   }
   ## unclear how best to define "color_by" column at this step
   ## exon labels
   if (exists("exonLabelDF") && length(exonLabelDF) > 0) {
      exonLabelDF$type <- "exon_label";
      exonLabelDF$name <- pasteByRow(exonLabelDF[,c("gr","sample_id")], sep=" ");
      exonLabelDF$feature <- exonLabelDF$gr;
      exonLabelDF$row <- seq_len(nrow(exonLabelDF));
      exonLabelDF$color_by <- NA;
      ## define name as factor to maintain the drawing order
      exonLabelDF$name <- factor(exonLabelDF$name,
         levels=unique(exonLabelDF$name));
   } else {
      exonLabelDF <- NULL;
   }
   ## junctions
   if (exists("juncDF") && length(juncDF) > 0) {
      juncDF$type <- "junction";
      juncDF$name <- pasteByRow(juncDF[,c("nameFromTo", "sample_id")], sep=" ");
      juncDF$feature <- juncDF$nameFromTo;
      juncDF$row <- seq_len(nrow(juncDF));
      junction_spans <- abs(
         shrinkMatrix(juncDF$x, groupBy=juncDF$name, min, returnClass="matrix") -
            shrinkMatrix(juncDF$x, groupBy=juncDF$name, max, returnClass="matrix"))[,1];
      juncDF$junction_span <- junction_spans[as.character(juncDF$name)];
      ## order the name column using junction_rank
      ## has affect on drawing order, making lower junction_rank
      ## drawn last, since they are typically small and otherwise
      ## easily obscured by the higher rank and larger junctions.
      if (length(unique(juncDF$junction_rank)) > 1) {
         junc_name_levels <- unique(
            mixedSortDF(juncDF,
               byCols=c("-junction_rank", "-junction_span", "row"))$name);
         juncDF$name <- factor(juncDF$name,
            levels=junc_name_levels);
      }
   } else {
      juncDF <- NULL;
   }
   ## junction labels
   if (exists("juncLabelDF") && length(juncLabelDF) > 0) {
      juncLabelDF$type <- "junction_label";
      juncLabelDF$name <- pasteByRow(juncLabelDF[,c("nameFromTo", "sample_id")], sep=" ");
      juncLabelDF$feature <- juncLabelDF$nameFromTo;
      juncLabelDF$row <- seq_len(nrow(juncLabelDF));
      ## define name as factor to maintain the drawing order
      juncLabelDF$name <- factor(juncLabelDF$name,
         levels=unique(juncLabelDF$name));
   } else {
      juncLabelDF <- NULL;
   }
   ## create a list of data.frames
   cjL <- list();
   cjL$coverage <- covDF;
   cjL$junction <- juncDF;
   cjL$junction_label <- juncLabelDF;
   cjL$exon_label <- exonLabelDF;
   ## Merge into one data.frame, then re-order
   if (length(cjL) == 0) {
      cjDF <- NULL;
   } else if (length(cjL) == 1) {
      cjDF <- cjL[[1]];
   } else {
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Merging ",
            length(cjL),
            " data.frames into df.");
      }
      cjDF <- jamba::mergeAllXY(cjL);
      cjDF <- mixedSortDF(cjDF,
         byCols=c("type","row"));
   }
   ## order columns by presence of NA values
   na_ct <- apply(cjDF, 2, function(i){
      sum(is.na(i))
   });
   cjDF <- cjDF[,order(na_ct),drop=FALSE];

   ## Add ref2c as an attribute to cjDF just to help keep it available
   attr(cjDF, "ref2c") <- ref2c;

   if (any(c("all", "df") %in% return_data)) {
      retVals$df <- cjDF;
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
#' @family RNA-seq functions
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
      jamba::printDebug("internal_junc_score(): ",
         "Defining flank ends of each junction.");
   }
   ## If there is no sampleColname column, ignore
   if (length(sampleColname) > 0 &&
         !sampleColname %in% colnames(GenomicRanges::values(juncGR))) {
      stop(paste0("The sampleColname '",
         sampleColname,
         "' is not present in colnames(values(juncGR))."));
   }
   juncEndsGR <- c(
      flank(juncGR[,c(scoreColname,sampleColname)],
         start=TRUE,
         width=1),
      flank(juncGR[,c(scoreColname,sampleColname)],
         start=FALSE,
         width=1));
   GenomicRanges::values(juncEndsGR)$side <- rep(c("start", "end"), each=length(juncGR));
   GenomicRanges::values(juncEndsGR)$id <- rep(seq_along(juncGR), 2);
   juncEndsDF <- as.data.frame(unname(juncEndsGR));
   juncGroupColnames <- intersect(
      c(sampleColname, "seqnames", "start", "strand", "side"),
      colnames(juncEndsDF));
   juncGroup <- pasteByRow(juncEndsDF[, juncGroupColnames, drop=FALSE]);
   juncEndsRed <- shrinkMatrix(juncEndsDF[[scoreColname]],
      groupBy=juncGroup,
      shrinkFunc=sum);
   if (verbose) {
      jamba::printDebug("internal_junc_score(): ",
         "dim(juncEndsRed):", dim(juncEndsRed));
   }
   juncEndsRefGR <- juncEndsGR[match(juncEndsRed$groupBy, juncGroup)];
   GenomicRanges::values(juncEndsRefGR)[[scoreColname]] <- juncEndsRed$x;
   if (verbose) {
      jamba::printDebug("internal_junc_score(): ",
         "length(juncEndsRefGR):", length(juncEndsRefGR));
   }

   fo1 <- GenomicRanges::findOverlaps(juncGR, juncEndsRefGR);
   fo1df <- data.frame(q=S4Vectors::queryHits(fo1),
      s=S4Vectors::subjectHits(fo1));
   if (verbose) {
      jamba::printDebug("internal_junc_score(): ",
         "dim(fo1df):", dim(fo1df));
   }
   if (nrow(fo1df) == 0) {
      intScore <- nameVector(
         rep(0, length(juncGR)),
         names(juncGR));
   } else {
      if (length(sampleColname) == 0) {
         fo1df$qSample <- rep("sample_id", length(fo1df$q));
         fo1df$sSample <- rep("sample_id", length(fo1df$q));
      } else {
         fo1df$qSample <- GenomicRanges::values(juncGR[fo1df$q])[[sampleColname]];
         fo1df$sSample <- GenomicRanges::values(juncEndsRefGR[fo1df$s])[[sampleColname]];
      }
      fo1df$sScore <- GenomicRanges::values(juncEndsRefGR[fo1df$s])[[scoreColname]];

      fo1dfuse <- subset(fo1df,
         qSample == sSample);
      if (verbose) {
         jamba::printDebug("internal_junc_score(): ",
            "dim(fo1dfuse):", dim(fo1dfuse));
      }
      if (nrow(fo1dfuse) == 0) {
         intScore <- nameVector(
            rep(0, length(juncGR)),
            names(juncGR));
      } else {
         intScore <- rmNA(naValue=0,
            max(
               List(
                  split(fo1dfuse$sScore, fo1dfuse$q)
               )
            )[as.character(seq_along(juncGR))]);
         names(intScore) <- names(juncGR);
      }
   }
   return(intScore);
}

#' Import splice junction data from BED file
#'
#' Import splice junction data from BED file
#'
#' This function is intended to be called internally by
#' `prepareSashimi()`, and is provided primarily to enable
#' use of `memoise::memoise()` to cache results.
#'
#' This function uses `rtracklayer::import()` to import data,
#' so any file or URL compatible with that function is acceptable.
#' Also, any compatible file format can be used, for example
#' a gzipped BED file. However, "bigBed" format cannot be used,
#' since the rtracklayer package does not support that format.
#'
#' @family jam data import functions
#'
#' @param iBed path or URL to one BED file containing splice
#'    junction data. The score is expected to be stored in the
#'    name column (column 3), primarily because the score column
#'    is sometimes restricted to maximum value 1000. However if
#'    the name column cannot be converted to numeric without
#'    creating NA values, then the score column will be used.
#' @param juncNames the name of the junction source file
#' @param sample_id character string representing the sample
#'    identifier.
#' @param scale_factor numeric value used to adjust the raw
#'    score, applied by multiplying the scale_factor by each score.
#' @param gr GRanges representing the overall range for which
#'    junction data will be retrieved. Note that any junctions
#'    that span this range, but do not start or end inside this
#'    range, will be removed.
#'
#' @export
import_juncs_from_bed <- function
(iBed,
 juncNames,
 sample_id,
 scale_factor,
 use_memoise=FALSE,
 memoise_junction_path="junctions_memoise",
 gr=NULL,
 verbose=FALSE,
 ...)
{
   # use_memoise=TRUE invokes slightly different cache strategy:
   # - load fill iBed data (using cached version if present)
   # - subset by range(gr)
   if (use_memoise) {
      import_m <- memoise::memoise(rtracklayer::import.bed,
         cache=memoise::cache_filesystem(memoise_junction_path));
      import_has_cache <- memoise::has_cache(import_m)(iBed);
      if (verbose) {
         jamba::printDebug("import_juncs_from_bed():",
            "import_has_cache:", import_has_cache);
      }
      bed1 <- tryCatch({
         import_m(iBed);
      }, error=function(e){
         warnText <- paste0("getGRcoverageFromBw():",
            "import_juncs_from_bed() error:",
            "file not accessible:'",
            iBed,
            "', returning NULL.");
         warning(warnText);
         NULL;
      })
      if (length(bed1) == 0) {
         jamba::printDebug("import_juncs_from_bed(): ",
            "Repairing junction cache.",
            fgText=c("darkorange","seagreen3"));
         if (compareVersion(as.character(packageVersion("memoise")), "1.1.0.900") >= 0) {
            #import_has_cache <- memoise::has_cache(import_m)(iBed);
            memoise::drop_cache(import_m)(iBed);
            bed1 <- import_m(iBed);
            bed1 <- tryCatch({
               import_m(iBed);
            }, error=function(e){
               NULL;
            })
         } else {
            bed1 <- tryCatch({
               if (length(gr) > 0) {
                  rtracklayer::import(iBed,
                     which=range(gr));
               } else {
                  rtracklayer::import(iBed);
               }
            }, error=function(e){
               NULL;
            });
         }
      }
      if (length(bed1) > 0 && length(gr) > 0) {
         bed1 <- subsetByOverlaps(bed1,
            ranges=range(gr));
      } else {
         jamba::printDebug("import_juncs_from_bed(): ",
            "Failed to Repair junction cache.",
            fgText=c("darkorange","red"));
         bed1 <- NULL;
         return(bed1);
      }
   } else {
      bed1 <- tryCatch({
         if (length(gr) > 0) {
            rtracklayer::import(iBed,
               which=range(gr));
         } else {
            rtracklayer::import(iBed);
         }
      }, error=function(e){
         NULL;
      })
   }
   if (length(bed1) == 0) {
      return(bed1);
   }

   ## Assign score and apply scale_factor
   ##
   ## By default if name has numeric values, use them as scores,
   ## since the score column is sometimes restricted to integer
   ## values with a maximum value 1000.
   if (!any(is.na(as.numeric(as.character(GenomicRanges::values(bed1)$name))))) {
      GenomicRanges::values(bed1)$score <- as.numeric(as.character(GenomicRanges::values(bed1)$name)) * scale_factor;
   } else {
      GenomicRanges::values(bed1)$score <- as.numeric(as.character(GenomicRanges::values(bed1)$score)) * scale_factor;
   }
   ## Assign annotation values
   GenomicRanges::values(bed1)[,c("juncNames")] <- juncNames;
   GenomicRanges::values(bed1)[,c("sample_id")] <- sample_id;
   ## Subset junctions to require either start or end to be contained
   ## within the region of interest (filters out phantom mega-junctions)
   bed1 <- subset(bed1,
      (
         IRanges::overlapsAny(GenomicRanges::flank(bed1, -1, start=TRUE), range(gr)) |
         IRanges::overlapsAny(GenomicRanges::flank(bed1, -1, start=FALSE), range(gr))
      )
   );
   bed1;
}

