
#' GRangesList to data.frame for ggplot2
#'
#' GRangesList to data.frame for ggplot2
#'
#' @export
grl2df <- function
(grl,
 keepGRvalues=TRUE,
 keepGRLvalues=FALSE,
 addGaps=TRUE,
 width=0.6,
 widthV=c(`exon`=0.6, `intron`=0.02, `gap`=0.02),
 width_colname="feature_type",
 shape=c("rectangle", "junction"),
 baseline=NULL,
 scoreColname="score",
 scoreMinimum=200,
 scoreFactor=1,
 scoreArcFactor=0.5,
 doStackJunctions=TRUE,
 ref2c=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert GRangesList to tall data.frame
   shape <- match.arg(shape);
   offset <- width / 2;
   if ("GRanges" %in% class(grl)) {
      grl <- GRangesList(list(gr=grl));
   }
   if (addGaps && !"junction" %in% shape) {
      grl <- addGRLgaps(grl,
         ...);
   }
   ## TODO: handle gaps, perhaps by calling this function with
   ## result of getGRLgaps(), but setting addGaps=FALSE, and using
   ## width=0.05
   if ("rectangle" %in% shape) {
      xCoords <- as.vector(rbind(
         start(grl@unlistData),
         end(grl@unlistData),
         end(grl@unlistData),
         start(grl@unlistData)
      ));
      width <- rep(width, length.out=length(grl@unlistData));
      if (width_colname %in% colnames(values(grl@unlistData))) {
         for (wName in names(widthV)) {
            wFound <- (values(grl@unlistData)[[width_colname]] %in% wName);
            if (any(wFound)) {
               printDebug("wName:", wName);
               width[values(grl@unlistData)[[width_colname]] %in% wName] <- widthV[wName];
            }
         }
      }
      offset <- rep(width, each=4) / 2;
      yBaseline <- rep(
         rep(seq_along(grl) - 1,
            elementNROWS(grl)),
         each=4);
      yCoords <- rep(c(-1, -1, 1, 1) * offset,
         length.out=length(xCoords)) +
         yBaseline;
      id <- rep(seq_along(grl@unlistData),
         each=4);
      df <- data.frame(x=xCoords,
         y=yCoords,
         id=id);
      if (length(names(grl)) > 0) {
         grlNames <- factor(
            rep(
               rep(names(grl),
                  elementNROWS(grl)),
               each=4),
            levels=names(grl));
         df$grlNames <- grlNames;
      }
      if (length(names(grl@unlistData)) > 0) {
         grNames <- factor(
            rep(
               names(grl@unlistData),
               each=4),
            levels=names(grl@unlistData));
         df$grNames <- grNames;
      }
      if (keepGRvalues) {
         for (iCol in colnames(values(grl@unlistData))) {
            df[,iCol] <- rep(values(grl@unlistData)[,iCol], each=4);
         }
      }
      if (keepGRLvalues) {
         for (iCol in colnames(values(grl))) {
            df[,iCol] <- rep(
               rep(values(grl)[,iCol],
                  elementNROWS(grl)),
               each=4);
         }
      }
   } else if ("junction" %in% shape) {
      #############################################################
      ## optionally stack junction ends so they do not overlap
      if (!all(scoreFactor == 1)) {
         if (verbose) {
            printDebug("grl2df(): ",
               "scoreFactor:",
               scoreFactor);
         }
         values(grl@unlistData)[[scoreColname]] <- (scoreFactor *
            values(grl@unlistData)[[scoreColname]])
      }
      if (doStackJunctions) {
         if (verbose) {
            printDebug("grl2df(): ",
               "stackJunctions()");
         }
         grlNew <- GRangesList(
            lapply(grl, function(iGR){
               stackJunctions(gr=iGR,
                  scoreColname=scoreColname,
                  baseline=baseline,
                  scoreFactor=1);
            })
         );
         grl <- grlNew;
      }
      ## create two polygons
      xStart <- start(grl@unlistData);
      xEnd <- end(grl@unlistData);
      if (length(ref2c) > 0) {
         xMid <- ref2c$inverse((ref2c$transform(xStart) + ref2c$transform(xEnd)) / 2);
      } else {
         xMid <- (xStart + xMid) / 2;
      }
      xCoords <- as.vector(rbind(
         xStart,
         xMid,
         xMid,
         xStart,
         xMid,
         xEnd,
         xEnd,
         xMid
      ));
      ## y-axis midpoints
      if ("yStart" %in% colnames(values(grl@unlistData))) {
         yStart <- values(grl@unlistData)$yStart;
      } else {
         yStart <- rep(0, length(grl@unlistData));
      }
      if ("yEnd" %in% colnames(values(grl@unlistData))) {
         yEnd <- values(grl@unlistData)$yEnd;
      } else {
         yEnd <- rep(0, length(grl@unlistData));
      }
      yScore <- values(grl@unlistData)[[scoreColname]] * 1;
      yMaxStartEnd <- pmax(yStart, yEnd);
      yHeight <- yScore * scoreArcFactor + scoreMinimum + yMaxStartEnd;
      yCoords <- as.vector(rbind(
         yStart,
         yStart + yHeight,
         yStart + yHeight + yScore,
         yStart + yScore,
         yStart + yHeight,
         yEnd,
         yEnd + yScore,
         yStart + yHeight + yScore
      ));
      ## make unique ids
      id <- rep(
         makeNames(
            rep(seq_along(grl@unlistData), each=2)),
         each=4);
      ## data.frame of coordinates
      df <- data.frame(x=xCoords,
         y=yCoords,
         id=id);
      if (length(names(grl)) > 0) {
         grlNames <- factor(
            rep(
               rep(names(grl),
                  elementNROWS(grl)),
               each=8),
            levels=names(grl));
         df$grlNames <- grlNames;
      }
      if (length(names(grl@unlistData)) > 0) {
         grNames <- factor(
            rep(
               names(grl@unlistData),
               each=8),
            levels=names(grl@unlistData));
         df$grNames <- grNames;
      }
      if (keepGRvalues) {
         for (iCol in colnames(values(grl@unlistData))) {
            df[,iCol] <- rep(values(grl@unlistData)[[iCol]],
               each=8);
         }
      }
      if (keepGRLvalues) {
         for (iCol in colnames(values(grl))) {
            df[,iCol] <- rep(
               rep(values(grl)[,iCol],
                  elementNROWS(grl)),
               each=8);
         }
      }
   }
   return(df);
}


#' Stack the y-axis position of junctions
#'
#' Stack the y-axis position of junctions
#'
#' @export
stackJunctions <- function
(gr,
 scoreColname="score",
 scoreFactor=1,
 baseline=NULL,
 ...)
{
   ## Purpose is to stack junctions by score (width) so they do
   ## not overlap at the start or end of each junction.
   if (!scoreColname %in% colnames(values(gr))) {
      stop("The scoreColname must be present in colnames(values(gr))");
   }

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   exonsFrom <- gsub("[.][-]*[0-9]+$", "",
      values(gr)[,"nameFrom"]);
   exonsTo <- gsub("[.][-]*[0-9]+$", "",
      values(gr)[,"nameTo"]);
   allExons <- mixedSort(unique(
      c(exonsFrom, exonsTo)));
   baselineV <- nameVector(
      rep(0,
         length.out=length(allExons)),
      allExons);
   if (length(baseline) > 0) {
      if (length(names(baseline)) > 0) {
         baselineV[names(baseline)] <- baseline;
      } else {
         baselineV[] <- baseline;
      }
   }

   ## Start position
   if (length(names(gr)) > 0) {
      grNames <- names(gr);
   }
   juncGR1 <- gr[do.call(order, list(start(gr), width(gr)))];
   exonsFrom1 <- gsub("[.][-]*[0-9]+$", "",
      values(juncGR1)[,"nameFrom"]);
   values(juncGR1)[,"yStart"] <- shrinkMatrix(
      values(juncGR1)$score * scoreFactor,
      groupBy=start(juncGR1),
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))})$x +
      baselineV[exonsFrom1];
   ## End position
   juncGR1 <- juncGR1[do.call(order, list(end(juncGR1), width(juncGR1)))];
   exonsTo1 <- gsub("[.][-]*[0-9]+$", "",
      values(juncGR1)[,"nameTo"]);
   values(juncGR1)[,"yEnd"] <- shrinkMatrix(
      values(juncGR1)$score,
      groupBy=end(juncGR1),
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))})$x +
      baselineV[exonsTo1];
   if (exists(grNames)) {
      juncGR1 <- juncGR1[grNames];
   } else {
      juncGR1 <- sort(juncGR1);
   }
   return(juncGR1);
}
