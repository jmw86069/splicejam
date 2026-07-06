




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
#' @param gr `GRanges` object containing regions not to compress. Regions
#'    which are unstranded gaps are compressed to fixed width.
#' @param gapWidth integer value used for fixed gap width, or when
#'    NULL the gap width is defined as 3 times the median feature width.
#' @param keepValues `logical` indicating whether to keep feature values
#'    in the GRanges data, default FALSE.
#' @param upstream,downstream,upstreamGapWidth,downstreamGapWidth
#'    `integer` number of bases to extend upstream and downstream the
#'    overall range of features, and the number of bases to compress
#'    that extended range, respectively. Any coordinates upstream
#'    by the 'upstream' amount will be compressed into 'upstreamGapWidth'
#'    bases in visual space. This ratio represents a visual compression
#'    multiplier. The compression is applied to upstream coordinate zero,
#'    and to downstream coordinate 10 billion.
#'    Defaults are 50000 for upstream and downstream, and three times
#'    the gap width `3 * gapWidth`.
#' @param nBreaks `integer` number of x-axis coordinate breaks used
#'    in ggplot labeling, default 7.
#' @param ignore.strand `logical` default TRUE, new as of 0.0.89.900,
#'    will create compressed coordinates by ignoring strand-specificity.
#'    Use FALSE for previous behavior, although it would be limited
#'    to scenarios where two features on opposite strands should not
#'    be used to produce one common set of "gaps" between the genome
#'    region features. In most cases, typical splicejam functionality
#'    is unaffected, but implies gaps should ignore strandedness.
#' @param verbose `logical` indicating whether to print verbose output.
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
 ignore.strand=TRUE,
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
   disGR <- GenomicRanges::reduce(gr,
      ignore.strand=ignore.strand);
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
      refCoord=c(GenomicRanges::start(disGR),
         GenomicRanges::end(disGR)),
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
         GenomicRanges::ranges(gr) <- IRanges::IRanges(
            start=GenomicRanges::values(gr)[,"refStart"],
            end=GenomicRanges::values(gr)[,"refEnd"]
         );
      } else {
         GenomicRanges::values(gr)[,"refStart"] <- GenomicRanges::start(gr);
         GenomicRanges::values(gr)[,"refEnd"] <- GenomicRanges::end(gr);
      }
      GenomicRanges::ranges(gr) <- IRanges::IRanges(
         start=ref2compressed(GenomicRanges::start(gr)),
         end=ref2compressed(GenomicRanges::end(gr))
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
   (x=NULL,
    limits=NULL,
    n=nBreaks,
    xFixed=lookupCoordDF[,1],
    verbose=FALSE,
    ...)
   {
      xvals <- jamba::rmNA(unique(sort(xFixed)));
      if (length(x) == 0) {
         x <- xFixed;
      }
      if (length(limits) == 0) {
         limits <- range(x, na.rm=TRUE)
      }
      xvals <- xvals[xvals >= limits[1] & xvals <= limits[2]];
      if (verbose) {
         jamba::printDebug("breaks_gr(): ",
            "x:", x);
         jamba::printDebug("breaks_gr(): ",
            "limits:", limits);
         jamba::printDebug("breaks_gr(): ",
            "n:", n);
      }
      if (n >= length(xvals)) {
         return(xvals);
      }
      idx1 <- round(seq(from=1, to=length(xvals), length.out=n));
      return(xvals[idx1]);
   }
   minor_breaks <- function
   (b=NULL,
    limits=NULL,
    n=5,
    xFixed=lookupCoordDF[,1],
    compressed=FALSE,
    verbose=FALSE,
    ...)
   {
      if (length(n) == 0) {
         n <- 5
      }
      if (verbose) {
         jamba::printDebug("minor_breaks(): ",
            "b:", b);
         jamba::printDebug("minor_breaks(): ",
            "limits:", limits);
         jamba::printDebug("minor_breaks(): ",
            "n:", n);
      }
      if (length(b) == 0) {
         b <- breaks_gr(limits=limits)
      }
      xFixed_range <- range(head(tail(xFixed, -1), -1));
      use_b <- b;
      do_comp <- FALSE;
      if (isTRUE(compressed)) {
         use_b <- compressed2ref(b)
         limits <- compressed2ref(limits)
      }
      nUse <- length(b) * n;
      minor_br <- breaks_gr(x=use_b,
         limits=limits,
         n=nUse,
         xFixed=xFixed,
         verbose=verbose);
      minor_br <- setdiff(minor_br, b);
      if (isTRUE(compressed)) {
         minor_br <- ref2compressed(minor_br);
      }
      minor_br <- setdiff(minor_br, b);
      return(minor_br)
   }
   labels <- function(x, n=10) {
      scales::comma_format()(breaks_gr(x, n=n))
   }
   retVals <- list();

   ## Make the custom trans function
   trans_grc <- scales::new_transform(name="compressed_gr",
      transform=ref2compressed,
      inverse=compressed2ref,
      breaks=function(x, n=7){breaks_gr(x=x, n=n)},
      minor_breaks=function(limits, b, n=5){minor_breaks(b=b, limits=limits, n=n, compressed=FALSE)},
      format=scales::comma_format(),
      domain=range(lookupCoordDF[,1]));

   ## Make the actual scale_x_gr_compressed() function
   scale_x_grc <- function(..., trans=trans_grc){
      ggplot2::scale_x_continuous(name="grc",
         # breaks=ggplot2::waiver(),#trans_grc$breaks,
         # minor_breaks=ggplot2::waiver(),#trans_grc$minor_breaks,
         breaks=trans$breaks,
         minor_breaks=trans$minor_breaks,
         labels=function(...){scales::comma(...)},#trans_grc$labels,
         trans=trans,
         ...)
   }

   ## Functions required by scales::trans_new()
   retVals$transform <- ref2compressed;
   retVals$inverse <- compressed2ref;
   retVals$breaks <- breaks_gr;
   retVals$minor_breaks <- minor_breaks;

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
      xyAngle <- jamba::rad2deg(xyAngle);
   }
   yRle <- S4Vectors::Rle(xyAngle);

   if (any(S4Vectors::runLength(yRle) >= minN)) {
      rl <- S4Vectors::runLength(yRle);
      rv <- S4Vectors::runValue(yRle);
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
         rlDF2 <- unique(jamba::mixedSortDF(
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
   #
   if (inherits(polyM[[1]], c("list", "AsIs"))) {
      # x,y columns are AsIs list of numeric vectors, one row per polygon
      data_style <- "fortify";
      iseq <- seq_len(nrow(polyM));
      polyMratios <- unlist(lapply(iseq, function(irow){
         xr1 <- range(polyM[irow, 1][[1]]);
         xc1 <- ref2c$transform(xr1);
         xv1 <- rev(sort(c(diff(xr1)+1, diff(xc1)+1)));
         xv1[1] / xv1[2];
      }))
      polyMlengths <- lengths(polyM[[1]]);
      #
      # polyDF <- polyM;
      # polyDF$newX <- (lapply(polyDF$x, function(ix){
      #    ref2c$transform(ix)
      # }))
      # polyDF$ratio <- c(NA, diff(polyDF$newX) / diff(polyDF$x));
      # which to adjust
      whichComp <- which(polyMratios >= minRatio & polyMlengths > 5);
      newpolyMs <- jamba::rbindList(lapply(whichComp, function(icomp){
         xvals <- unname(polyM[icomp, 1][[1]])
         yvals <- unname(polyM[icomp, 2][[1]])
         
         xr1 <- range(polyM[icomp, 1][[1]]);
         xc1 <- ref2c$transform(xr1);
         xuse <- seq_along(xvals)[c(-1, -1 * length(xvals) + c(1, 0))]
         xseq <- round(seq(from=xc1[1]+0.0,
            to=xc1[2]-0.0,
            length.out=length(xvals) - 3))
         # new x - use mean x-value
         xseq1 <- sapply(split(xvals[xuse], xseq), mean, na.rm=TRUE);
         xseq1[1] <- head(xvals[xuse], 1);
         xseq1[length(xseq1)] <- tail(xvals[xuse], 1)
         # new y - use max y-value
         yseq1 <- sapply(split(yvals[xuse], xseq), max, na.rm=TRUE);
         newy <- c(0, yseq1, 0, 0)
         newx <- c(xseq1[1], xseq1, tail(xseq1, 1), xseq1[1])
         newpolyM <- polyM[icomp, , drop=FALSE];
         newpolyM[1, 1][[1]] <- list(newx);
         newpolyM[1, 2][[1]] <- list(newy);
         newpolyM
      }))
      polyM[whichComp, ] <- newpolyMs;
      return(polyM);
      #
   } else if (any(is.na(polyM[,1]))) {
      data_style <- "base";
      polyMlengths <- diff(unique(c(0, which(is.na(polyM[,1])), nrow(polyM))));
      polyMratios <- shrinkMatrix(polyDF$ratio,
         groupBy=polyDF$n,
         shrinkFunc=function(x){median(x, na.rm=TRUE)})$x;
   } else {
      data_style <- "fortify";
      idrows <- jamba::pasteByRow(polyM[, c("cov", "gr"), drop=FALSE]);
      idrows <- factor(idrows, levels=unique(idrows));
      polyML <- split(polyM, idrows);
      polyMratios <- unlist(lapply(polyML, function(i){
         xr1 <- range(i[,1]);
         xc1 <- ref2c$transform(xr1);
         xv1 <- rev(sort(c(diff(xr1)+1, diff(xc1)+1)));
         xv1[1] / xv1[2];
      }));
      polyMlengths <- jamba::sdim(polyML)[,1];
   }
   if (verbose > 1) {
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
   sdim_polyDFL <- jamba::sdim(polyDFL)$rows;
   whichComp <- which(polyMratios >= minRatio & sdim_polyDFL > 5);
   whichNorm <- which(polyMratios < minRatio | sdim_polyDFL <= 5);
   if (length(whichComp) == 0) {
      return(polyM);
   }
   ## Compress those polygons whose coordinates get compressed
   polyDFLnew <- lapply(jamba::nameVector(whichComp), function(k){
      iDF <- polyDFL[[k]];
      baseline <- iDF[1,"y"];
      iDF <- iDF[!is.na(iDF[,1]), , drop=FALSE];
      iDFu <- iDF[match(unique(iDF$x), iDF$x), , drop=FALSE];

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
      iSeqSub[iSeqSub > max(iSeqNew)] <- max(iSeqNew);

      ## use approx() to fill in holes
      iYnew <- approx(x=iDFu$x,
         y=iDFu$y,
         xout=iSeqNew)$y;
      ## Take running max value, then use approx
      iYrunmax <- caTools::runmax(iYnew,
         k=iMultiple,
         endrule="keep");
      if (verbose > 2) {
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
   newPolyDFL[names(polyDFL)[whichNorm]] <- lapply(jamba::nameVector(whichNorm),
      function(k){
      polyDFL[[k]][, c("x", "y", "cov", "gr"), drop=FALSE];
   });
   newPolyDFL[names(polyDFLnew)] <- polyDFLnew;
   newPolyDFL <- newPolyDFL[jamba::mixedSort(names(newPolyDFL))];
   newPolyM <- jamba::rbindList(newPolyDFL);
   return(newPolyM);
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
#' @return `list` containing `ggSashimi` a ggplot2 graphical object
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
#' @param sjenv `environment` (or `list`) containing argument values
#'    indicated below. This input is intended to be compatible with
#'    output from `sashimiDataConstants()`. It currently recognizes
#'    the following, although function arguments take priority when
#'    provided.
#'    * 'flatExonsByGene'
#'    * 'ref2c' (optional)
#'    * 'filesDF'
#' @param flatExonsByGene `GRangesList` named by gene, specifically
#'    named by 'gene_name' as in the 'tx2geneDF' `data.frame`.
#'    `GRanges` elements are disjoint, non-overlapping genomic ranges
#'    per gene and should represent only exon regions. All non-exon
#'    regions are assumed to be introns or gaps.
#' @param filesDF `data.frame` with columns:
#'    * `url`: any valid path compatible with `data.table::fread()`
#'    * `sample_id`: `character` string representing a biological sample,
#'    and often may be used for sample grouping, to combine coverage
#'    and/or junction read counts for all samples in the group.
#'    * `type`: `character` string with 'bw' for bigwig coverage,
#'    'junction' for BED12 splice junctions, or 'coverageGR' when
#'    coverage data are supplied as `GRanges` objects named by 'sample_id',
#'    using argument 'covGR'.
#' @param gene `character` string of the gene to prepare, which must be
#'    present in `names(flatExonsByGene)`.
#' @param sample_id `character` default NULL, used to subset entries
#'    in 'filesDF' to include only specific matching sample_id entries.
#'    When NULL, it uses all entries.
#' @param minJunctionScore `numeric` default 10, to require junctions
#'    to have at least this total score per sample_id to be displayed.
#'    This filter it useful to hide spurious junction counts which may
#'    only have very few supporting read counts, and is typically
#'    adjusted relative to the total coverage of the gene.
#'    For example, 'GAPDH' gene often has extremely high coverage,
#'    and may be associated with numerous spurious junctions with
#'    10 to 50 reads, which are not useful to display when coverage
#'    is in the 10,000 order of magnitude.
#' @param gapWidth `numeric` value of the fixed width to use for
#'    gaps (introns) between exon features. If `NULL` then
#'    `getGRgaps()` will use the default based upon the median exon
#'    width.
#' @param addGaps `logical` default TRUE, whether to include gap regions
#'    in the coverage plot, for example including introns or intergenic
#'    regions. When `compressGR=TRUE` then gaps regions are
#'    down-sampled using running maximum signal with roughly the same
#'    x-axis resolution as uncompressed regions.
#' @param baseline `numeric` vector named by `names(flatExonsByGene)`
#'    where baseline is used to adjust the y-axis baseline position
#'    above or below zero. Default zero 0 for all practical purposes.
#' @param compressGR `logical` default TRUE, whether to compress GRanges
#'    coordinates in the output data, where gaps/introns are set
#'    to a fixed width. When `ref2c` is not supplied, and
#'    `compressGR=TRUE`, then `ref2c` is created using
#'    `make_ref2compressed()`.
#' @param compress_introns `logical`default TRUE, whether to compress
#'    the coverage polygon coordinates to approximately the same
#'    number of pixels per inch as the exon polygons. This option
#'    greatly reduces the size of the polygon, since introns are
#'    already about 50 to 100 times wider than exons, and when
#'    `compressGR` is `TRUE`, the introns are visibly compressed
#'    to a fixed width on the x-axis. The data has many more
#'    x-axis coordinates than the data visualization, this argument
#'    is intended to reduce the intron coordinates accordingly.
#' @param ref2c `list` object output from `make_ref2compressed()`
#'    used to compress axis coordinates, to compress polygon coverage
#'    data in compressed regions, and to adjust splice junction arcs
#'    using compressed coordinates.
#' @param gap_feature_type `character` string, default 'intron',
#'    with the feature_type representing gaps when `addGaps=TRUE`.
#' @param default_feature_type `character` string, default 'exon'
#'    with the default feature_type representing non-gap regions.
#' @param feature_type_colname `character` string with the column
#'    name of `values()` to represent the feature type (exon/intron).
#' @param exon_label_type `character` string indicating the type
#'    of label for exons:
#'    * 'none' (default): no exon label is displayed
#'    * 'repel': use ggrepel to display the exon label
#'    * 'mark': use `ggforce::geom_mark_rect()` to display exon label
#' @param junc_label_type `character` string indicating the type
#'    of label for junction counts:
#'    * 'repel' (default): use ggrepel to display the exon label
#'    * 'mark': use `ggforce::geom_mark_rect()` to display exon label
#'    * 'none': no exon label is displayed
#' @param return_data `character` string, default 'df', with format
#'    of data to return: 'df' is `data.frame` and 'ref2c' returns
#'    only the 'ref2c' axis compression as `list`. Typically only 'df'
#'    is used.
#' @param include_strand `character` string, default 'both' with
#'    strandedness of data to include. Coverage strandedness is
#'    determined by filename, assuming some substring which matches
#'    'plus|pos|+' represents positive strand. By convention in
#'    splicejam, scores are negative for negative strand coverage,
#'    and negative strand junction counts, however not required.
#'    * `"both"`: include coverage and junctions from both strands,
#'    * `"+"`: include coverage and junctions from only the '+' strand,
#'    * `"-"`: include coverage and junctions from only the '-' strand.
#' @param junc_color,junc_fill `character` valid R colors used for
#'    junction ribbon arcs. The alpha transparency is maintained,
#'    and when junctions are stacked, the actual color is slightly
#'    adjusted to vary from light to dark in order of stacked
#'    junction read counts. Default uses 'goldenrod3' border, and
#'    'goldenrod1' for the fill color.
#' @param doStackJunctions `logical` default TRUE, whether to stack
#'    junction arcs at each end, this argument is passed to
#'    `grl2df()` which calls `stackJunctions()`.
#'    The purpose of stacking junctions is to place each ribbon
#'    junction arc starting above the previous ribbon arc, for each
#'    junction start, and junction end. In this way, the total
#'    height of junction ribbon widths should represent the total
#'    junction read count/score at each junction boundary.
#' @param coord_method `character` string indicating how the x-axis
#'    coordinates are represented for use in ggplot2 plot functions.
#'    * 'coord' (default): stores coordinates unchanged, without
#'    regard to intron/gap compression. It allows the axis to be
#'    displayed using ggplot2 techniques, optionally using 'ref2c'
#'    to compress introns/gaps visually, keeping numeric values as-is.
#'    * 'scale': stores numeric coordinates after compression when
#'    'compressGR=TRUE', which means these values are not genome
#'    coordinates, instead are display coordinates. The ggplot2
#'    axis labels are expected to be adjusted using 'ref2c' for
#'    display purposes.
#'    * 'none': genome coordinates are stored as-is with no anticipated
#'    modifications by 'ref2c'. This approach is effectively
#'    equivalent to 'coord'.
#' @param scoreFactor `numeric` default 1, scalar multiplied by each
#'    score, expected to be provided as named `numeric` vector
#'    named by 'sample_id'.
#' @param scoreArcFactor `numeric` default 0.2, used to adjust the
#'    magnitude of splice ribbon arc curvature. Use higher values to
#'    increase the height of the midpoint of the ribbon arc relative
#'    to each end of the ribbon. This value is a scalar, roughly applied
#'    as '1 + scoreArcFactor' such that 0.2 has the effect of 
#'    'score * 1.2', a 20% addition.
#' @param scoreArcMinimum `numeric` default 100, the minimum height of
#'    the midpoint of the junction ribbon arc above the flat mean value
#'    between the start and end of the ribbon. It is used as minimum,
#'    when the `scoreArcFactor` provide too low a value. Useful when
#'    there may be intermediate coverage, or intervening exon, to have
#'    the arc curve above the coverage.
#' @param covGR `GRanges` default NULL, optional coverage data in columns
#'    stored as `NumericList`, where `colnames(GenomicRanges::values(covGR))`
#'    are present in 'filesDF' for rows with
#'    `filesDF$type %in% "coverage_gr"`. The name should match 'url' or
#'    'sample_id' in 'filesDF'.
#' @param juncGR `GRanges` default NULL, optional splice junctions,
#'    `'score'` is used for the abundance of splice junction reads,
#'    and `'sample_id'` is used to define the biological sample
#'    as defined in 'filesDF'.
#' @param use_memoise `logical` default FALSE, passed to
#'    `getGRcoverageFromBw()` and `import_juncs_from_bed()` as needed,
#'    indicating whether to use memoise to cache intermediate data
#'    results. Default is TRUE for `launchSashimiApp()` and FALSE
#'    otherwise. In general, TRUE would be a valid default option.
#' @param memoise_coverage_path,memoise_junction_path `character`
#'    string with default folder path to store memoise cache files,
#'    used with `use_memoise=TRUE`.
#' @param do_shiny_progress `logical` default FALSE, whether to send
#'    progress updates to a running shiny app using the
#'    `shiny::withProgress()` and `shiny::setProgress()` methods.
#'    This function calls `shiny::setProgress()` and
#'    assumes that `shiny::withProgress()` has already been
#'    initialized. It also does not close the progress, instead
#'    pushes that responsibility to the calling function.
#' @param verbose `logical` whether to print verbose output.
#' @param ... additional arguments are passed to `make_ref2compressed()`,
#'    `getGRcoverageFromBw()`, `exoncov2polygon()`.
#'
#' @family jam RNA-seq functions
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
(sjenv=NULL,
 flatExonsByGene=NULL,
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
 exon_label_type=c("none",
    "repel",
    "mark"),
 junc_label_type=c("repel",
    "mark",
    "none"),
 return_data=c("df",
    "ref2c"),
 include_strand=c("both",
    "+",
    "-"),
 junc_color=jamba::alpha2col("goldenrod3", 0.7),
 junc_fill=jamba::alpha2col("goldenrod1", 0.4),
 doStackJunctions=TRUE,
 coord_method=c("coord",
    "scale",
    "none"),
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
   asSeconds <- function(x, digits=2, ...){
      format(difftime(as.POSIXct(x), 0), digits=digits, ...)
   }

   # recognize environment as input
   if (inherits(sjenv, "environment")) {
      # To consider:
      # - flatExonsByTx, tx2geneDF, detectedTx -> new flatExonsByGene
      if ('flatExonsByGene' %in% ls(sjenv) &&
         length(flatExonsByGene) == 0) {
         if (verbose) {
            jamba::printDebug("prepareSashimi(): ",
               "Using ", "sjenv$flatExonsByGene")
         }
         flatExonsByGene <- sjenv$flatExonsByGene;
      }
      if ('ref2c' %in% ls(sjenv) &&
         length(ref2c) == 0) {
         if (verbose) {
            jamba::printDebug("prepareSashimi(): ",
               "Using ", "sjenv$ref2c")
         }
         ref2c <- sjenv$ref2c;
      }
      if ('filesDF' %in% ls(sjenv) &&
         length(filesDF) == 0) {
         if (verbose) {
            jamba::printDebug("prepareSashimi(): ",
               "Using ", "sjenv$filesDF")
         }
         filesDF <- sjenv$filesDF;
      }
   }
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
   if (!jamba::igrepHas("GRangesList", class(flatExonsByGene))) {
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
      if (any(grepl("pos|[+]|plus", ignore.case=TRUE, bwFilesDF$url))) {
         if ("+" %in% include_strand) {
            bwFilesDF <- subset(bwFilesDF,
               grepl("pos|[+]|plus",
                  ignore.case=TRUE,
                  basename(url)));
         } else {
            bwFilesDF <- subset(bwFilesDF,
               !grepl("pos|[+]|plus",
                  ignore.case=TRUE,
                  basename(url)));
         }
      }
   }
   bwUrls <- jamba::nameVector(
      bwFilesDF[, c("url","url"), drop=FALSE]);
   bwSamples <- jamba::nameVector(
      bwFilesDF[, c("sample_id", "url"), drop=FALSE]);
   bwUrlsL <- split(bwUrls, unname(bwSamples));
   if (!"scale_factor" %in% colnames(bwFilesDF)) {
      bwFilesDF$scale_factor <- rep(1, nrow(bwFilesDF));
   }
   bwScaleFactors <- jamba::rmNA(naValue=1,
      bwFilesDF$scale_factor);
   names(bwScaleFactors) <- names(bwUrls);
   if (length(bwScaleFactors) > 0) {
      if (verbose > 1) {
         jamba::printDebug("prepareSashimi(): ",
            "bwScaleFactors:",
            paste0(basename(names(bwScaleFactors)),
               ":",
               format(bwScaleFactors, digits=2, trim=TRUE)),
            sep=", ");
      }
   }
   if (verbose > 1) {
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
            stop(paste0("Supplied coverage covGR does not have colnames ",
               "in filesDF$url with filesDF$type == 'coverage_gr'"));
         }
         covUrls <- jamba::nameVector(
            covfilesDF[, c("url", "url"), drop=FALSE]);
         covSamples <- jamba::nameVector(
            covfilesDF[, c("sample_id", "url"), drop=FALSE]);
         covGRuse <- covGRuse[, covUrls];
         if ("scale_factor" %in% colnames(covfilesDF)) {
            covScaleFactors <- jamba::rmNA(naValue=1,
               covfilesDF$scale_factor);
         } else {
            covScaleFactors <- rep(1, length.out=length(covUrls));
         }

         ## Combine coverage per strand
         if (do_shiny_progress) {
            ##
            shiny::setProgress(1/4,
               detail=paste0("Combining bw coverage by sample_id"));
         }
         covGR2 <- combineGRcoverage(covGRuse,
            covName=covSamples,
            scaleFactors=covScaleFactors,
            covNames=names(covUrls),
            verbose=verbose > 1);
         ## Obtain the new set of covNames
         covNames <- attr(covGR2, "covNames");
         covName <- attr(covGR2, "covName");
         ## some_null is TRUE when any underlying coverage file failed to return coverage
         some_null <- attr(covGR2, "some_null");
         if (length(some_null) && some_null) {
            retVals$some_null <- some_null;
         }

         ## 0.0.82.900 - add some verbose output
         if (verbose > 1) {
            jamba::printDebug("prepareSashimi(): ",
               "covName:");
            print(covName);
            jamba::printDebug("prepareSashimi(): ",
               "covNames:");
            print(covNames);
            jamba::printDebug("prepareSashimi(): ",
               "some_null:", some_null);
         }

         ## Create polygon data.frame
         if (verbose) jamba::printDebug("prepareSashimi(): ", "Calling exoncov2polygon() #1");
         st7 <- system.time({
            covDF <- exoncov2polygon(covGR2,
               ref2c=ref2c,
               covNames=covNames,
               sample_id=covName,
               coord_style="fortify",
               compress_introns=compress_introns,
               verbose=verbose > 1,
               ...);
         })
         if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st7["elapsed"]));
         ## Enforce ordered factor levels for sample_id
         ## which also forces empty factor levels if applicable
         if (length(covDF) > 0) {
            covDF$sample_id <- factor(
               as.character(covDF$sample_id),
               levels=unique(sample_id)
            );
         }
         if (any(c("all", "covDF") %in% return_data)) {
            retVals$covDF <- covDF;
         }
         ########################################
         ## Optional exon labels
         exonLabelDF <- NULL;
         if (length(covDF) > 0 && nrow(covDF) > 0) {
            covDFsub <- (as.character(covDF$gr) %in% names(gr));
            covDFlab <- subset(covDF, covDFsub);

            exonLabelDF1 <- shrinkMatrix(covDFlab[, c("x", "y"), drop=FALSE],
               groupBy=jamba::pasteByRowOrdered(
                  covDFlab[, c("gr", "sample_id"), drop=FALSE], sep=":!:"),
               shrinkFunc=function(x){mean(range(x))});
            exonLabelDF1[,c("gr","sample_id")] <- jamba::rbindList(
               strsplit(as.character(exonLabelDF1$groupBy), ":!:"));
            exonLabelDF1$gr <- factor(exonLabelDF1$gr,
               levels=unique(exonLabelDF1$gr));
            exonLabelDF <- jamba::renameColumn(exonLabelDF1,
               from="groupBy",
               to="gr_sample");
         }
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
            "Preparing bw coverage data for ", gene);
      }
      if (do_shiny_progress) {
         ##
         shiny::setProgress(1/4,
            detail=paste0("Preparing bw coverage data for ", gene));
      }
      ## Note that coverage is not scaled at this step
      st15 <- system.time({
         covGR <- getGRcoverageFromBw(gr=gr,
            bwUrls=bwUrls,
            addGaps=addGaps,
            gap_feature_type=gap_feature_type,
            default_feature_type=default_feature_type,
            feature_type_colname=feature_type_colname,
            use_memoise=use_memoise,
            memoise_coverage_path=memoise_coverage_path,
            do_shiny_progress=do_shiny_progress,
            verbose=verbose > 1,
            ...);
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st15["elapsed"]));
      ## Combine coverage per strand
      if (verbose > 1) {
         jamba::printDebug("prepareSashimi(): ",
            "Combining bw coverage by sample_id");
      }
      if (do_shiny_progress) {
         ##
         shiny::setProgress(1/4,
            detail=paste0("Combining bw coverage by sample_id"));
      }
      st16 <- system.time({
         covGR2 <- combineGRcoverage(covGR,
            covName=bwSamples,
            scaleFactors=bwScaleFactors,
            covNames=names(bwUrls),
            verbose=verbose > 1);
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st16["elapsed"]));
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
            "Converting coverage to polygons.");
      }
      if (do_shiny_progress) {
         ##
         shiny::setProgress(2.5/4,
            detail=paste0("Converting coverage to polygons."));
      }
      st8 <- system.time({
         covDF <- exoncov2polygon(covGR2,
            ref2c=ref2c,
            covNames=covNames,
            sample_id=covName,
            coord_style="fortify",
            compress_introns=compress_introns,
            verbose=verbose > 1,
            ...);
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st8["elapsed"]));
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
      if (verbose > 1) {
         jamba::printDebug("prepareSashimi(): ",
            "head(covDF$sample_id):");
         print(head(covDF$sample_id));
      }

      ########################################
      ## Optional exon labels
      covDFsub <- (as.character(covDF$gr) %in% names(gr));
      covDFlab <- covDF[covDFsub, , drop=FALSE];

      if (verbose) jamba::printDebug("prepareSashimi(): ", "Creating exon labels.");#
      # st9a <- system.time({
      #    DT <- data.table(
      #       data.frame(check.names=FALSE,
      #          stringsAsFactors=FALSE,
      #          covDFlab[, c("x", "y", "gr", "sample_id"), drop=FALSE]),
      #       key=c("gr", "sample_id"))
      #    shrinkFunc <- function(x){mean(range(x))}
      #    byDT <- DT[,lapply(.SD, shrinkFunc),
      #       by=c("gr", "sample_id")]
      # })
      # if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st9a["elapsed"]));
      st9 <- system.time({
         exonLabelDF1 <- shrinkMatrix(covDFlab[, c("x", "y"), drop=FALSE],
            groupBy=jamba::pasteByRowOrdered(
               covDFlab[,c("gr", "sample_id"), drop=FALSE], sep=":!:"),
            shrinkFunc=function(x){mean(range(x))});
         exonLabelDF1[,c("gr","sample_id")] <- jamba::rbindList(
            strsplit(as.character(exonLabelDF1$groupBy), ":!:"));
         exonLabelDF1$gr <- factor(exonLabelDF1$gr,
            levels=unique(exonLabelDF1$gr));
         exonLabelDF <- jamba::renameColumn(exonLabelDF1,
            from="groupBy",
            to="gr_sample");
         if (any(c("all", "covDF") %in% return_data)) {
            retVals$exonLabelDF <- exonLabelDF;
         }
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st9["elapsed"]));

   }

   ############################################
   ## Load Junctions
   ##
   ## Pre-existing junctions supplied as GRanges
   if (length(juncGR) > 0) {
      if (!"sample_id" %in% colnames(GenomicRanges::values(juncGR))) {
         stop("juncGR must have 'sample_id' in colnames(GenomicRanges::values(juncGR)).");
      }
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Preparing junctions from juncGR");
      }
      if (do_shiny_progress) {
         ##
         shiny::incProgress(3/4,
            detail=paste0("Preparing GR junction data for ", gene));
      }
      juncGRuse <- juncGR[GenomicRanges::values(juncGR)[["sample_id"]] %in% sample_id];
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
      filesDF$sample_id %in% sample_id,
      c("url", "sample_id", "scale_factor"), drop=FALSE];
   juncUrls <- jamba::nameVector(
      juncFilesDF[, c("url", "sample_id"), drop=FALSE]);
   juncSamples <- jamba::nameVector(
      juncFilesDF[, c("sample_id", "sample_id"), drop=FALSE]);
   juncUrlsL <- split(juncUrls, juncSamples);
   if (!"scale_factor" %in% colnames(juncFilesDF)) {
      juncFilesDF$scale_factor <- rep(1, nrow(juncFilesDF));
   }
   juncScaleFactors <- jamba::rmNA(naValue=1,
      jamba::nameVector(
         juncFilesDF[, c("scale_factor", "sample_id"), drop=FALSE]));
   if (verbose > 1 && length(juncScaleFactors) > 0) {
      jamba::printDebug("prepareSashimi(): ",
         "juncScaleFactors: ",
         paste0(names(juncScaleFactors),
            ":",
            format(juncScaleFactors, digits=2, trim=TRUE)),
         sep=", ");
   }
   if (verbose > 1) {
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
         shiny::incProgress(3/4,
            detail=paste0("Importing junctions for ", gene));
      }
      if (verbose) jamba::printDebug("prepareSashimi(): ", "Importing junctions for ", gene);
      st10 <- system.time({
         juncBedList <- lapply(jamba::nameVectorN(juncUrls), function(iBedName){
            iBed <- juncUrls[[iBedName]];
            iBedNum <- match(iBedName, jamba::makeNames(names(juncUrls)));
            iBedPct <- (iBedNum - 1) / length(juncUrls);
            if (do_shiny_progress) {
               if (!is.na(iBedPct)) {
                  shiny::setProgress(
                     value=3/4 + iBedPct/6,
                     detail=paste0("Importing junctions (",
                        iBedNum,
                        " of ",
                        length(juncUrls),
                        ") for ", gene));
               }
            }
            if (verbose > 1) {
               jamba::printDebug("prepareSashimi(): ",
                  "progress:", format(digits=2, 3/4 + iBedPct/6),
                  " Importing bed:",
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
               if (verbose > 1) {
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
                  use_memoise=TRUE,
                  memoise_junction_path=memoise_junction_path,
                  gr=gr);
            } else {
               bed1 <- import_juncs_from_bed(iBed,
                  juncNames=iBedName,
                  sample_id=juncSamples[iBedName],
                  scale_factor=1,
                  gr=gr);
            }
            ## Apply scale_factor
            if (length(bed1) > 0) {
               GenomicRanges::values(bed1)$score <- (
                  GenomicRanges::values(bed1)$score * juncScaleFactors[iBedName]);
            }
            bed1;
         });
      });
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st10["elapsed"]));
      juncBedList <- juncBedList[lengths(juncBedList) > 0];
      if (length(juncBedList) == 0) {
         juncDF1 <- NULL;
      } else {
         juncBedGR <- GenomicRanges::GRangesList(juncBedList)@unlistData;
         # 23jan2022 bug with duplicate unlistData names in rare edge cases
         # causes problem when coercing with as.data.frame()
         names(juncBedGR) <- jamba::makeNames(names(juncBedGR));
         ## Create junction summary data.frame
         if (verbose > 1) {
            jamba::printDebug("prepareSashimi(): ",
               "running spliceGR2junctionDF for juncBedGR()");
         }
         if (verbose) jamba::printDebug("prepareSashimi(): ", "Combining junction data for ", gene);
         if (do_shiny_progress) {
            shiny::setProgress(3.7/4,
               detail=paste0("Combining junction data for ", gene));
         }
         st11 <- system.time({
            juncDF1f <- spliceGR2junctionDF(spliceGRgene=juncBedGR,
               exonsGR=gr,
               sampleColname="sample_id");
         })
         if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st11["elapsed"]));
         if (length(juncDF1) > 0) {
            if (verbose > 1) {
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
         jamba::renameColumn(juncDF1,
            from="ref",
            to="seqnames"),
         "GRanges");
      names(juncGR) <- jamba::makeNames(
         GenomicRanges::values(juncGR)[, "nameFromToSample"]);

      if (length(baseline) == 0) {
         baseline <- 0;
      }
      if (verbose) {
         jamba::printDebug("prepareSashimi(): ",
            "Calculating junction stacking for ", gene);
      }
      if (do_shiny_progress) {
         shiny::setProgress(3.8/4,
            detail=paste0("Calculating junction stacking for ", gene));
      }
      ## Convert junctions to polygons usable by geom_diagonal_wide()
      st12 <- system.time({
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
            verbose=verbose > 1,
            ...);
      });
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st12["elapsed"]));
      if (!"sample_id" %in% colnames(juncDF)) {
         juncDF <- jamba::renameColumn(juncDF,
            from="grl_name",
            to="sample_id");
      }
      if (verbose > 1) {
         jamba::printDebug("prepareSashimi(): ",
            "dim(juncDF):", dim(juncDF));
      }
      juncDF <- juncDF[, !colnames(juncDF) %in% c("grl_name"), drop=FALSE];
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
      if (verbose) jamba::printDebug("prepareSashimi(): ", "Preparing junction label coordinates.");
      #juncLabelDF1 <- subset(mutate(juncCoordDF, id_name=jamba::makeNames(id)), grepl("_v1_v3$", id_name));
      if (do_shiny_progress) {
         ##
         shiny::setProgress(3.9/4,
            detail=paste0("Preparing junction label coordinates for ", gene));
      }
      st13 <- system.time({
         juncLabelDF1 <- subset(
            plyr::mutate(juncDF, id_name=jamba::makeNames(id)),
            grepl("_v1_v[23]$", id_name));
         shrink_colnames <- intersect(c("x", "y", "score", "junction_rank"),
            colnames(juncLabelDF1));
         ## Define junction placement at max position, in stranded fashion
         juncLabelDF_y <- jamba::renameColumn(
            shrinkMatrix(juncLabelDF1[, "y", drop=FALSE],
               shrinkFunc=function(x){max(abs(x))*sign(max(x))},
               groupBy=juncLabelDF1[, "nameFromToSample"]),
            from="groupBy",
            to="nameFromToSample");
         juncLabelDF <- jamba::renameColumn(
            shrinkMatrix(juncLabelDF1[, shrink_colnames, drop=FALSE],
               groupBy=juncLabelDF1[, "nameFromToSample"]),
            from="groupBy",
            to="nameFromToSample");
         juncLabelDF$y <- juncLabelDF_y$y;
         juncLabelDF[, c("nameFromTo", "sample_id")] <- jamba::rbindList(
            strsplit(juncLabelDF[, "nameFromToSample"], ":!:"));
         juncLabelDF[,c("nameFrom", "nameTo")] <- jamba::rbindList(
            strsplit(juncLabelDF[, "nameFromTo"], " "));
         ## Enforce ordered factor levels for sample_id
         ## which also forces empty factor levels if applicable
         juncLabelDF$sample_id <- factor(
            as.character(juncLabelDF$sample_id),
            levels=unique(sample_id)
         );
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st13["elapsed"]));
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
      covDF$name <- jamba::pasteByRow(
         covDF[, c("gr", "cov", "sample_id"), drop=FALSE], sep=" ");
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
      exonLabelDF$name <- jamba::pasteByRow(
         exonLabelDF[, c("gr", "sample_id"), drop=FALSE], sep=" ");
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
      juncDF$name <- jamba::pasteByRow(
         juncDF[, c("nameFromTo", "sample_id"), drop=FALSE], sep=" ");
      juncDF$feature <- juncDF$nameFromTo;
      juncDF$row <- seq_len(nrow(juncDF));
      junction_spans <- abs(
         shrinkMatrix(juncDF$x,
            groupBy=juncDF$name, min, returnClass="matrix") -
         shrinkMatrix(juncDF$x,
            groupBy=juncDF$name, max, returnClass="matrix"))[,1];
      juncDF$junction_span <- junction_spans[as.character(juncDF$name)];
      ## order the name column using junction_rank
      ## has affect on drawing order, making lower junction_rank
      ## drawn last, since they are typically small and otherwise
      ## easily obscured by the higher rank and larger junctions.
      if (length(unique(juncDF$junction_rank)) > 1) {
         junc_name_levels <- unique(
            jamba::mixedSortDF(juncDF,
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
      juncLabelDF$name <- jamba::pasteByRow(
         juncLabelDF[, c("nameFromTo", "sample_id"), drop=FALSE], sep=" ");
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
   cjL <- cjL[jamba::sdim(cjL)$rows > 0];
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
      st14 <- system.time({
         cjL <- lapply(cjL, function(ijL){
            if (!inherits(ijL$x, c("AsIs", "list"))) {
               ijL$x <- I(as.list(ijL$x));
            }
            if (!inherits(ijL$y, c("AsIs", "list"))) {
               ijL$y <- I(as.list(ijL$y));
            }
            if ("junction_rank" %in% colnames(ijL) &&
               inherits(ijL$junction_rank, c("factor"))) {
               ijL$junction_rank <- as.numeric(ijL$junction_rank);
            }
            ijL
         })
         
         cjDF <- cjL[[1]];
         
         for (ij in 2:length(cjL)) {
            # print(ij);# debug
            # print(jamba::sdim(cjL[[ij]]));# debug
            cjDF <- suppressMessages({
               dplyr::full_join(cjDF, cjL[[ij]])
            })
         }
         # cjDF <- jamba::mergeAllXY(cjL, by=NULL);
         cjDF <- jamba::mixedSortDF(cjDF,
            byCols=c("type","row"));
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st14["elapsed"]));
   }
   ## order columns by presence of NA values
   na_ct <- apply(cjDF, 2, function(i){
      sum(is.na(i))
   });
   cjDF <- cjDF[, order(na_ct), drop=FALSE];

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
#' @family jam RNA-seq functions
#'
#' @return numeric vector named by `names(juncGR)` whose values
#'    are the maximum score of internal overlapping junction ends.
#'
#' @param juncGR GRanges containing splice junctions, with numeric
#'    column `scoreColname` representing the abundance of splice
#'    junction-spanning sequence reads.
#' @param scoreColname,sampleColname colnames in `values(juncGR)`
#'    which define the junction score, and the `"sample_id"`.
#' @param minScore optional `numeric` value indicating the minimum
#'    score to use, which can be useful for example if there is
#'    no matching `scoreColname`. Note that values are matched using
#'    absolute value, but the original sign is maintained.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
internal_junc_score <- function
(juncGR,
 scoreColname="score",
 sampleColname="sample_id",
 minScore=0,
 verbose=FALSE,
   ...)
{
   ## Purpose is to sum scores for junctions which overlap
   ## junction ends inside the junction boundary.
   if (verbose) {
      jamba::printDebug("internal_junc_score(): ",
         "Defining flank ends of each junction.");
      jamba::printDebug("internal_junc_score(): ",
         "sampleColname:",
         sampleColname);
      jamba::printDebug("internal_junc_score(): ",
         "scoreColname:",
         scoreColname);
   }
   if (!scoreColname %in% colnames(GenomicRanges::values(juncGR))) {
      if (length(minScore) > 0) {
         if (verbose) {
            jamba::printDebug("internal_junc_score(): ",
               "Defined temporary score column '",
               scoreColname,
               "' with minScore=",
               minScore);
         }
         GenomicRanges::values(juncGR)[[scoreColname]] <- rep(minScore,
            length.out=length(juncGR));
      } else {
         stop(paste0("The score column ",
            scoreColname,
            "' was not found, and no minScore was provided."));
      }
   }
   if (length(minScore) > 0) {
      if (any(abs(GenomicRanges::values(juncGR)[[scoreColname]]) < abs(minScore))) {
         GenomicRanges::values(juncGR)[[scoreColname]] <- jamba::noiseFloor(
            abs(GenomicRanges::values(juncGR)[[scoreColname]]),
            minimum=abs(minScore)) *
            sign(GenomicRanges::values(juncGR)[[scoreColname]] + 1e-99);
         if (verbose) {
            jamba::printDebug("internal_junc_score(): ",
               "Applied minScore:",
               minScore);
         }
      }
   }
   ## If there is no sampleColname column, ignore
   if (length(sampleColname) > 0 &&
         !sampleColname %in% colnames(GenomicRanges::values(juncGR))) {
      stop(paste0("The sampleColname '",
         sampleColname,
         "' is not present in colnames(GenomicRanges::values(juncGR))."));
   }
   juncEndsGR <- c(
      GenomicRanges::flank(juncGR[,c(scoreColname,sampleColname)],
         start=TRUE,
         width=1),
      GenomicRanges::flank(juncGR[,c(scoreColname,sampleColname)],
         start=FALSE,
         width=1));
   GenomicRanges::values(juncEndsGR)$side <- rep(c("start", "end"),
      each=length(juncGR));
   GenomicRanges::values(juncEndsGR)$id <- rep(seq_along(juncGR), 2);
   juncEndsDF <- as.data.frame(unname(juncEndsGR));
   juncGroupColnames <- intersect(
      c(sampleColname, "seqnames", "start", "strand", "side"),
      colnames(juncEndsDF));
   juncGroup <- jamba::pasteByRow(juncEndsDF[, juncGroupColnames, drop=FALSE]);
   juncEndsRed <- shrinkMatrix(juncEndsDF[[scoreColname]],
      groupBy=juncGroup,
      shrinkFunc=sum);
   if (verbose) {
      jamba::printDebug("internal_junc_score(): ",
         "dim(juncEndsRed):", dim(juncEndsRed));
   }
   juncEndsRefGR <- juncEndsGR[match(juncEndsRed$groupBy, juncGroup)];
   GenomicRanges::values(juncEndsRefGR)[[scoreColname]] <- juncEndsRed$x;
   if (verbose > 1) {
      jamba::printDebug("internal_junc_score(): ",
         "length(juncEndsRefGR):", length(juncEndsRefGR));
   }

   fo1 <- GenomicRanges::findOverlaps(juncGR, juncEndsRefGR);
   fo1df <- data.frame(q=S4Vectors::queryHits(fo1),
      s=S4Vectors::subjectHits(fo1));
   if (verbose > 1) {
      jamba::printDebug("internal_junc_score(): ",
         "dim(fo1df):", dim(fo1df));
   }
   if (nrow(fo1df) == 0) {
      intScore <- jamba::nameVector(
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
      if (verbose > 1) {
         jamba::printDebug("internal_junc_score(): ",
            "dim(fo1dfuse):", dim(fo1dfuse));
      }
      if (nrow(fo1dfuse) == 0) {
         intScore <- jamba::nameVector(
            rep(0, length(juncGR)),
            names(juncGR));
      } else {
         intScore <- jamba::rmNA(naValue=0,
            max(
               S4Vectors::List(
                  split(fo1dfuse$sScore, fo1dfuse$q)
               )
            )[as.character(seq_along(juncGR))]);
         names(intScore) <- names(juncGR);
      }
   }
   return(intScore);
}

#' Import splice junction data from BED or SJ.out.tab file
#'
#' Import splice junction data from BED or SJ.out.tab file
#'
#' This function is intended to be called internally by
#' `prepareSashimi()`, and is provided primarily to enable
#' use of `memoise::memoise()` to cache results.
#'
#' This function was refactored in version 0.0.69.900 to
#' handle either BED format, or `"SJ.out.tab"` junction
#' format as produced by STAR alignment. The method uses
#' `data.table::fread()` then if there are 9 columns,
#' it assumes the format is `"SJ.out.tab"`. Otherwise it
#' coerces the `data.frame` with `as(bed, "GRanges")`.
#'
#' The BED or `"SJ.out.tab"` file can be gzipped, provided
#' `data.table()` is able to recognieze and import the
#' compression format.
#'
#' Note that bigBed format still cannot be used since the
#' rtracklayer package does not support that format.
#'
#' Also note that junctions with `score=0` are dropped at
#' this step, to prevent propagation of junctions with
#' zero counts.
#'
#'
#' @family jam data import functions
#'
#' @param iBed `character` path or URL to one BED file containing splice
#'    junction data. The score is typially expected to be stored in the
#'    name column (column 3), primarily because the score column
#'    is restricted to maximum value 1000 when used for UCSC tracks.
#'
#'    This process also recognizes names in the form "JUNC000000_1234"
#'    where the read depth/score is interpreted as "1234". When all
#'    entries have this name convention, the values are converted to
#'    numeric and used to populate the "score" column.
#'
#'    If not all values in the "name" column can be converted to numeric
#'    without introducing NA values, the "score" column is used as-is.
#'
#'    If you see scores with maximum value "1000" the "name" field is
#'    probably not being used properly as a numeric score.
#' @param juncNames `character` the name of the junction source file
#' @param sample_id character string representing the sample
#'    identifier.
#' @param scale_factor `numeric` value used to adjust the raw
#'    score, applied by multiplying the scale_factor by each score.
#' @param use_memoise `logical` default FALSE, whether to cache data
#'    using memoise.
#' @param memoise_junction_path `character` default 'junctions_memoise'
#'    with default subdirectory for memoise cache files, used only
#'    when `use_memoise=TRUE`.
#' @param gr `GRanges` representing the overall range for which
#'    junction data will be retrieved. Note that any junctions
#'    that span this range, but do not start or end inside this
#'    range, will be removed.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @import data.table
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
   #
   # New workflow:
   # - import into data.frame
   # - determine BED or SJ.out.tab format
   # - call appropriate method

   # If given GRangesList, use the GRanges portion.
   # - design decision -- the other option is error,
   #   since this is unexpected input.
   if (length(gr) > 0 && inherits(gr, "GRangesList")) {
      gr <- gr@unlistData;
   }
   # helper function to convert data.frame to junctions GRanges
   df_to_junc_gr <- function(bed_df){
      bed_gr <- GenomicRanges::GRanges(
         seqnames="chr1",
         ranges=IRanges::IRanges(
            start=1, end=2),
         strand="*",
         score=1)[0];
      if (nrow(bed_df) == 0) {
         return(bed_gr);
      }
      if (ncol(bed_df) == 9) {
         # SJ.out format
         bed_df <- subset(bed_df, !bed_df[,7] %in% 0);
         bed_gr <- GenomicRanges::GRanges(
            seqnames=c(bed_df[,1]),
            ranges=IRanges::IRanges(
               start=bed_df[,2] - 1,
               end=bed_df[,3] + 0),
            strand=ifelse(bed_df[,4] == 1, "+",
               ifelse(bed_df[,4] == 2, "-", "*")),
            score=bed_df[,7]
         )
      } else {
         k <- seq_len(min(c(ncol(bed_df), 12)))
         colnames(bed_df)[k] <- c("seqnames",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts")[k];
         # check strand/score columns are switched
         if (!all(bed_df$strand %in% c("+", "-", "*")) &&
               all(bed_df$score %in% c("+", "-", "*")) &&
               all(is.numeric(bed_df$strand))) {
            bed_df[, c("score", "strand")] <- bed_df[, c("strand", "score"),
               drop=FALSE];
         }
         # spot-check the strand column
         is_strand <- apply(head(bed_df, 200), 2, function(i){
            all(i %in% c("*", "+", "-"))})
         if (!"strand" %in% names(which(is_strand))) {
            stop("Strand is not in the correct BED column.");
            # is_strand_k <- head(which(is_strand), 1) + c(0, 1);
            # bed_df <- jamba::renameColumn(bed_df,
            #    from=colnames(bed_df)[is_strand_k],
            #    to=colnames(bed_df)[rev(is_strand_k)])
         }
         # convert to GRanges
         bed_gr <- as(bed_df, "GRanges")
         # remove rows with score==0
         if ("score" %in% colnames(bed_df)) {
            bed_gr <- subset(bed_gr, !score %in% 0);
         }
      }
      return(bed_gr)
   }

   ######################################################################
   ## import_junc similar to import_or_null used for bigwig coverage
   # For files/URLs ending 'bb' or 'bigbed' it uses cpp11bigwig
   # for the range requested in gr.
   # Otherwise data.table::fread() for the whole file.
   bed_is_bigbed <- grepl("bb|bigbed$",ignore.case=TRUE, iBed);
   import_junc_or_null <- function(bed, gr=NULL) {
      # accept bigbed format
      if (grepl("(bb|bigbed)$",ignore.case=TRUE, bed)) {
         # use cpp11bigwig
         # use unique() to permit multiple gr to return duplicates
         # the other option is to use range(gr)
         junc_df <- unique(data.frame(check.names=FALSE,
            cpp11bigwig::read_bigbed(bbfile=bed,
               chrom=range(gr))))
               # chrom=GenomicRanges::range(gr)))
      } else {
         junc_df <- data.table::fread(file=bed,
            data.table=FALSE,
            showProgress=FALSE)
      }
      if (ncol(junc_df) < 4) {
         warn_msg <- paste(sep="\n",
            "Junction file only contains 3 columns. It is expected to have",
            "12 columns, in bed12 format.",
            "For bigbed files, the autoSql schema must be encoded",
            "into the file, for example like this:",
            "bedToBigBed -as=bed12.as junc.bed chromsizes.txt junc.bb",
            "",
            "An example bed12.as file is shown below:",
            'table bed12',
            '"Browser extensible data, with extended fields for detail page"',
            '    (',
            '    string chrom;      "Reference sequence chromosome or scaffold"',
            '    uint   chromStart; "Start position in chromosome"',
            '    uint   chromEnd;   "End position in chromosome"',
            '    string name;       "Short Name of item"',
            '    uint   score;      "Score from 0-1000"',
            '    char[1] strand;    "+ or -"',
            '    uint thickStart;   "Start of where display should be thick (start codon)"',
            '    uint thickEnd;     "End of where display should be thick (stop codon)"',
            '    uint reserved;     "Used as itemRgb as of 2004-11-22"',
            '    int blockCount;    "Number of blocks"',
            '    int[blockCount] blockSizes; "Comma separated list of block sizes"',
            '    int[blockCount] chromStarts; "Start positions relative to chromStart"'
         )
         warning(warn_msg);
         junc_df <- NULL;
      }
      return(junc_df);
   }

   bed_df <- NULL;
   if (use_memoise) {
      # function(x){data.table::fread(x, data.table=FALSE, showProgress=FALSE)},
      import_m <- memoise::memoise(
         import_junc_or_null,
         cache=memoise::cache_filesystem(memoise_junction_path));
      # for bigbed we include GR in the cache
      # otherwise we cache the whole file anyway
      if (isTRUE(bed_is_bigbed)) {
         import_has_cache <- memoise::has_cache(import_m)(iBed, gr=gr);
      } else {
         import_has_cache <- memoise::has_cache(import_m)(iBed, gr=NULL);
      }
      if (verbose) {
         jamba::printDebug("import_juncs_from_bed():",
            "import_has_cache:", import_has_cache);
      }
      bed_df <- tryCatch({
         if (isTRUE(bed_is_bigbed)) {
            import_m(iBed, gr=gr);
         } else {
            import_m(iBed, gr=NULL)
         }
      }, error=function(e){
         warnText <- paste0(
            "import_juncs_from_bed() error:",
            "file not accessible: '",
            iBed,
            "', returning NULL.");
         print(warnText);
         print(e);
         warning(warnText);
         NULL;
      });
      if (length(bed_df) == 0) {
         if (verbose) {
            jamba::printDebug("import_juncs_from_bed(): ",
               "Repairing junction cache.",
               fgText=c("darkorange","seagreen3"));
         }
         if (isTRUE(bed_is_bigbed)) {
            memoise::drop_cache(import_m)(iBed, gr=gr);
         } else {
            memoise::drop_cache(import_m)(iBed, gr=NULL)
         }
         bed_df <- tryCatch({
            if (isTRUE(bed_is_bigbed)) {
               import_m(iBed, gr=gr);
            } else {
               import_m(iBed, gr=NULL)
            }   
         }, error=function(e){
            warnText <- paste0(
               "import_juncs_from_bed() error:",
               "file not accessible during repair junction step: '",
               iBed,
               "', returning NULL.");
            print(warnText);
            print(e);
            warning(warnText);
            NULL;
         })
         if (length(bed_df) == 0) {
            jamba::printDebug("import_juncs_from_bed(): ",
               "Failed to Repair junction cache.",
               fgText=c("darkorange","red"));
         }
      }
   } else {
      bed_df <- tryCatch({
         if (verbose) {
            jamba::printDebug("import_juncs_from_bed(): ",
               "Calling import_junc_or_null()",
               fgText=c("darkorange","seagreen3"));
            # print(class(iBed));# debug
            # print(class(gr));# debug
         }
         suppressMessages({
            import_junc_or_null(bed=iBed, gr=gr)
         })
      }, error=function(e){
         jamba::printDebug("Error in import_junc_or_null():");
         print(e);
         NULL;
      });
   }
   # convert to junction GRanges
   if (length(bed_df) > 0) {
      bed1 <- df_to_junc_gr(bed_df);
   } else {
      bed1 <- NULL;
   }
   if (length(bed1) == 0) {
      return(NULL)
   }

   if (FALSE && length(bed1) > 0 && length(gr) > 0) {
      bed1 <- IRanges::subsetByOverlaps(bed1,
         ranges=range(gr));
   }

   ## Assign score and apply scale_factor
   ##
   ## By default if name has numeric values, use them as scores,
   ## since the score column is sometimes restricted to integer
   ## values with a maximum value 1000.
   ##
   ## 10mar2026: try to repair entries using JUNC000001_123 format
   ## which is interpreted as 123.
   if ("name" %in% colnames(GenomicRanges::values(bed1)) &&
         all(grepl("^JUNC[0-9.]+_[0-9.]+$",
            as.character(GenomicRanges::values(bed1)$name)))) {
      if (isTRUE(verbose)) {
         jamba::printDebug("import_juncs_from_bed(): ",
            "Converting JUNC000001_123 to 123 for score");
      }
      GenomicRanges::values(bed1)$score <- as.numeric(
         gsub("^JUNC[0-9.]+_", "",
            as.character(GenomicRanges::values(bed1)$name)));
   }
   if ("name" %in% colnames(GenomicRanges::values(bed1)) &&
         !any(is.na(as.numeric(as.character(GenomicRanges::values(bed1)$name))))) {
      if (isTRUE(verbose)) {
         jamba::printDebug("import_juncs_from_bed(): ",
            "Converting name directly to score");
      }
      GenomicRanges::values(bed1)$score <- as.numeric(
         as.character(GenomicRanges::values(bed1)$name)) * scale_factor;
   } else {
      if (isTRUE(verbose)) {
         jamba::printDebug("import_juncs_from_bed(): ",
            "Using score as-is.");
      }
      GenomicRanges::values(bed1)$score <- as.numeric(
         as.character(GenomicRanges::values(bed1)$score)) * scale_factor;
   }
   
   ## Assign annotation values
   GenomicRanges::values(bed1)[,c("juncNames")] <- juncNames;
   GenomicRanges::values(bed1)[,c("sample_id")] <- sample_id;
   ## Subset junctions to require either start or end to be contained
   ## within the region of interest (filters out phantom mega-junctions)
   if (length(gr) > 0) {
      bed1 <- subset(bed1,
         (
            IRanges::overlapsAny(GenomicRanges::flank(bed1, -1, start=TRUE),
               range(gr)) |
            IRanges::overlapsAny(GenomicRanges::flank(bed1, -1, start=FALSE),
               range(gr))
         )
      );
   }
   bed1;
}

