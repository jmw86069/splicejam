
#' GRangesList to data.frame for ggplot2
#'
#' GRangesList to data.frame for ggplot2
#'
#' This function central to other plotting functions for
#' GRanges and GRangesList objects. It currently has two
#' modes:
#'
#' * `shape="rectangle"` is intended for
#' exons/peaks/regions, for example plotting exons in a gene structure,
#' or plotting ChIP-seq peaks.
#' * `shape="junction"` is intended for splice junctions, and
#' returns two polygons that join junction ends with a the middle
#' point raised above both baselines. The polygon height is determined
#' by the score, resulting in visual reinforcement of the number
#' of splice junction reads, usually compared to the sequence read
#' coverage at adjacent exons.
#'
#' An interesting argument is `baseline` which can be a named vector
#' of baseline y-axis values for each GRanges entry in the `grl`
#' GRangesList object. For example, it can be used to shift exons
#' up or down on the y-axis to make alternative exons more visibly
#' distinct. When used for Sashimi plots, it should also be
#' supplied to `prepareSashimi()` or `exoncov2polygon()` so
#' the coverages and splice junctions have consistent y-axis baselines.
#'
#' When chromosome coordinates are compressed (to reduce the visible
#' width of introns) it affects the midpoint of splice junction arcs,
#' therefore `ref2c` should be supplied so the arcs are defined
#' using compresssed coordinates.
#'
#' @return data.frame with `x,y` coordinates, and `id` which is used
#'    to group polygon coordinates when used with `ggplot2::geom_polygon()`
#'    or `ggforce::geom_shape()`. When `shape="rectangle"` the colnames
#'    include `grlNames` which are names of the input GRangesList
#'    `names(grl)`; `grNames` which are names of the GRanges entries; and
#'    other columns from the input GRanges entries. When `shape="junction"`
#'    the data includes two polygons per junction, intended to be used
#'    with `ggforce::geom_diagonal_wide()` for each side in order to
#'    produce a ribbon arc. The data also includes `sample_id` which is
#'    helpful for keeping data distinct when derived from multiple
#'    samples.
#'
#' @param grl GRangesList, or GRanges which will be converted to a
#'    GRangesList of length=1.
#' @param keepGRvalues,keepGRLvalues logical indicating whether the
#'    output data.frame should include column values from GRangesList and
#'    GRanges, if available.
#' @param addGaps logical indicating whether to add gap GRanges between
#'    same-strand GRanges features within each GRangesList element.
#'    When `TRUE` the gaps will be drawn between each GRanges rectangle.
#' @param width numeric value of default width for features when
#'    `shape="rectangle"`.
#' @param widthV numeric vector whose names are column values, using values
#'    from the first available colname from `width_colname`. Some common
#'    defaults are provided. Values are suggested to be between 0 and 1,
#'    since each GRangesList element is separated by 1 y-axis unit.
#' @param width_colname when `widthV` is used to determine width based
#'    upon a column value, the first matching colname of `width_colname`
#'    is used.
#' @param shape character string indicating whether input data should
#'    be processed as rectangular segments or splice junction arcs.
#' @param scoreColname,scoreMinimum,scoreFactor,scoreArcFactor numeric
#'    values used to determine junction ribbon height, the minimum
#'    height of the arc above the starting y-axis values based upon
#'    the score, the scaling factor for score values, and the
#'    relative height of the arc above the starting y-axis values
#'    multiplied by the score.
#' @param doStackJunctions logical indicating whether to stack junctions
#'    at the start and end of junctions sharing the same coordinate,
#'    in order of shortest to longest junction width.
#' @param ref2c optional output from `make_ref2compressed()` used to
#'    compress axis coordinates during junction arc calculations.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are passsed to relevant downstream
#'    functions.
#'
#' @family jam GRanges functions
#' @family jam plot functions
#'
#' @examples
#' gr <- GRanges(seqnames=rep(c("chr1"), 7),
#'    ranges=IRanges(start=c(50, 100, 1300, 2500, 23750, 24900, 25000),
#'       end=c(100, 150, 1450, 2600, 23800, 25000, 25200)),
#'    strand=rep("+", 7),
#'    feature_type=rep(c("noncds", "cds", "noncds"), c(1,5,1)));
#' names(gr) <- makeNames(rep("exon", 7));
#' grldf <- grl2df(gr, addGaps=TRUE);
#'
#' gg1 <- ggplot2::ggplot(grldf, ggplot2::aes(x=x, y=y, group=id)) +
#'    ggforce::geom_shape(ggplot2::aes(fill=feature_type)) +
#'    colorjam::theme_jam()
#' print(gg1);
#'
#' ## For fun, compress the introns and plot again.
#' ## This method uses x-axis breaks at the exon boundaries.
#' ref2c <- make_ref2compressed(gr);
#' gg2 <- gg1 +
#'    ggplot2::scale_x_continuous(trans=ref2c$trans_grc) +
#'    colorjam::theme_jam()
#' print(gg2);
#'
#' ## data can also be plotted using coord_trans()
#' ## the main difference is that x-axis breaks are defined before the
#' ## transformation, which can result in non-optimal placement
#' gg3 <- gg1 +
#'    coord_trans(x=ref2c$trans_grc) + colorjam::theme_jam();
#' print(gg3);
#'
#' @export
grl2df <- function
(grl,
 keepGRvalues=TRUE,
 keepGRLvalues=FALSE,
 addGaps=TRUE,
 width=0.6,
 widthV=c(exon=0.6, cds=0.6, noncds=0.3, intron=0.01, gap=0.01),
 width_colname=c("subclass", "feature_type"),
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
   ## TODO: draw arrows on gap regions, optionally at the end of
   ## multi-rectangle features, to indicate strandedness.
   if ("rectangle" %in% shape) {
      xCoords <- as.vector(rbind(
         start(grl@unlistData),
         end(grl@unlistData),
         end(grl@unlistData),
         start(grl@unlistData)
      ));
      width <- rep(width, length.out=length(grl@unlistData));
      width_colname_use <- head(
         intersect(width_colname,
            colnames(values(grl@unlistData))),
         1);
      if (length(width_colname_use) > 0) {
         for (wName in names(widthV)) {
            wFound <- (values(grl@unlistData)[[width_colname_use]] %in% wName);
            if (any(wFound)) {
               printDebug("wName:", wName);
               width[values(grl@unlistData)[[width_colname_use]] %in% wName] <- widthV[wName];
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

#' Gene GRangesList to ggplot2 grob
#'
#' Gene GRangesList to ggplot2 grob
#'
#' This function is intended to help plot gene and transcript exon
#' models, and is a lightweight wrapper around `grl2df()`.
#'
#' It takes `flatExonsByGene` which is the output from
#' `flattenExonsByGene()`, and essentially plots the end result
#' for review.
#'
#' Alternatively, when `return_type="df"`, the output is
#' the `data.frame` used to produce the ggplot, which allows
#' for more customization.
#'
#' @family jam plot functions
#'
#' @param gene character string of the gene to plot, compared
#'    with `names(flatExonsByGene)` and `values(flatExonsByTx)$gene_name`.
#' @param tx character vector of the transcripts to plot, useful
#'    when specifying specific transcripts. Values are matched with
#'    `names(flatExonsByTx)`.
#' @param flatExonsByGene,flatExonsByTx GRangesList objects, named
#'    by `"gene_name"` or `"transcript_id"` respectively, containing
#'    disjoint (non-overlapping) exons within each GRangesList
#'    element. The data is expected to be in the form provided by
#'    `flattenExonsByGene()`.
#' @param geneColor character color used as the base color for
#'    exons, where the color is varied for each feature type or
#'    subclass.
#' @param labelExons logical indicating whether to print text
#'    labels beneath each exon, using the values in colname
#'    `"gene_nameExon"`. Typically the gene and transcripts are
#'    named using consistent names, in which case one exon label
#'    is placed at the bottom of the lowest transcript for each
#'    unique exon label.
#' @param exonLabelAngle numeric angle in degrees (0 to 360)
#'    indicating how to rotate exon labels, where `90` is
#'    vertical, and `0` is horizontal.
#' @param newValues argument passed to `addGRLgaps()` to fill
#'    column values for newly created gap entries. It is useful
#'    to have `feature_type="gap"` so gaps have a different value
#'    than exons. It is also useful to have `subclass="gap"`
#'    when there are `"cds"` and `"noncds"` entries in the
#'    provided `flatExonsByGene` data.
#' @param gene_order character value indicating whether the
#'    flattened gene model should be plotted `"first"` above the
#'    transcript exon models, or `"last"` and below the
#'    transcript exon models.
#' @param return_type character value indicating whether to return
#'    the ggplot graphic object `"grob"`, or the data.frame
#'    `"df"` used to create the ggplot object.
#' @param ref2c list output from `make_ref2compressed()` which
#'    contains among other things, the `trans_grc` data of
#'    class `"trans"` used in `ggplot2::coord_trans()` or
#'    `ggplot2::scale_x_continuous()`.
#' @param hjust,vjust numeric value to position exon labels
#'    using `ggplot::geom_text()`.
#' @param compressGaps logical indicating whether to compress gaps
#'    between exons. When `ref2c` is supplied, this argument is
#'    ignored and the supplied `ref2c` is used directly.
#' @param ... additional arguments are passed to relevant functions
#'    as needed, including `make_ref2compressed()`.
#'
#' @examples
#' ## Assume we start with flattened gene exons
#' if (1 == 2) {
#'    ## Do not run automated examples until sample data is available
#'
#' ggGria1 <- gene2gg("Gria1",
#'    flatExonsByGene=flatExonsByGeneCds);
#'
#' ## if transcript exons are available
#' ggGria1 <- gene2gg("Gria1",
#'    flatExonsByGene=flatExonsByGene,
#'    flatExonsByTx=flatExonsByTx);
#'
#' @export
gene2gg <- function
(gene=NULL,
 tx=NULL,
 flatExonsByGene=NULL,
 flatExonsByTx=NULL,
 geneColor="dodgerblue",
 labelExons=TRUE,
 exonLabelAngle=70,
 newValues=list(feature_type="gap", subclass="gap", gene_nameExon="gap"),
 gene_order=c("first","last"),
 return_type=c("grob", "df"),
 ref2c=NULL,
 hjust=1.2,
 vjust=0,
 compressGaps=TRUE,
 ...)
{
   ## Purpose is a lightweight wrapper around grl2df() specifically intended
   ## for gene and transcript exon structure
   gene_order <- match.arg(gene_order);
   return_type <- match.arg(return_type);
   if (length(flatExonsByGene) > 0 && length(gene) > 0) {
      if (!any(gene %in% names(flatExonsByGene))) {
         stop("gene was not found in names(flatExonsByGene)");
      }
      grl1a <- flatExonsByGene[gene];
   } else {
      grl1a <- NULL;
   }
   if (length(flatExonsByTx) > 0) {
      grl1 <- NULL;
      if (length(gene) > 0) {
         grl1 <- subset(flatExonsByTxCds2, gene_name %in% gene);
      }
      if (length(tx) > 0) {
         grl1 <- c(grl1,
            flatExonsByTxCds2[names(flatExonsByTxCds2) %in% tx]);
      }
   } else {
      grl1 <- NULL;
   }
   if ("first" %in% gene_order) {
      grl1a1 <- GRangesList(c(grl1, grl1a));
   } else {
      grl1a1 <- GRangesList(c(grl1a, grl1));
   }
   printDebug("class(grl1a1):", class(grl1a1));
   if (length(grl1a1) == 0) {
      stop("no exon models found for the gene and tx arguments given.");
   }
   grl1a1df <- grl2df(grl1a1,
      newValues=newValues,
      ...);

   ## ggplot2
   colorColname <- head(
      intersect(c("subclass", "feature_type"),
         colnames(grl1a1df)),
      1);
   subclassV <- provigrep(c("noncds", "utr", "cds", "exon", "intron", "gap", "."),
      unique(grl1a1df[[colorColname]]));
   colorSubV <- color2gradient(
      nameVector(rep(geneColor, length.out=length(subclassV)),
         subclassV));
   #showColors(colorSubSubclass);
   #colorSubSubclass <- c(cds="navy", noncds="dodgerblue", gap="grey30");
   ## Make a data.frame to label each exon
   if (labelExons) {
      exonLabelDF <- renameColumn(
         from="groupBy",
         to="gene_nameExon",
         shrinkMatrix(grl1a1df[,c("x","y")],
            groupBy=grl1a1df[,"gene_nameExon"]));
      exonLabelDF$y <- shrinkMatrix(grl1a1df[,c("x","y")],
         groupBy=grl1a1df[,"gene_nameExon"],
         shrinkFunc=min)$y;
      exonLabelDF <- subset(exonLabelDF, !grepl("^gap$|,", gene_nameExon));
   }
   if (length(ref2c) == 0) {
      if (compressGaps) {
         ref2c <- make_ref2compressed(grl1a1@unlistData,
            ...);
      } else {
         ref2c <- NULL;
      }
   }
   ## Put it together
   grl1a1gg <- ggplot(grl1a1df,
         aes(x=x, y=y, group=id)) +
      geom_shape(show.legend=FALSE,
         aes(fill=subclass, color=subclass)) +
      theme_jam() +
      ylab("") +
      scale_color_manual(values=makeColorDarker(colorSubV, darkFactor=1.3)) +
      scale_fill_manual(values=alpha2col(colorSubV, alpha=1)) +
      scale_y_continuous(breaks=seq_along(grl1a1)-1,
         limits=c(-0.7, length(grl1a1)-0.5),
         labels=names(grl1a1));
   if (length(ref2c) > 0) {
      grl1a1gg <- grl1a1gg + scale_x_continuous(trans=ref2c$trans_grc);
   }
   if (labelExons) {
      if (1 == 1 || length(ref2c) == 0) {
         grl1a1gg <- grl1a1gg +
            ggrepel::geom_text_repel(
               inherit.aes=FALSE,
               data=exonLabelDF,
               aes(x=x, y=min(y), label=gene_nameExon),
               angle=exonLabelAngle,
               hjust=vjust,
               vjust=hjust,
               segment.color="grey35",
               fill="white",
               size=3,
               direction="x");
      } else {
         grl1a1gg <- grl1a1gg +
            geom_text(
               inherit.aes=FALSE,
               data=exonLabelDF,
               aes(x=x, y=min(y), label=gene_nameExon),
               angle=exonLabelAngle,
               hjust=hjust,
               vjust=vjust,
               fill="white",
               size=3,
               direction="y");
      }
   }

   if (length(gene) > 0) {
      grl1a1gg <- grl1a1gg +
         ggtitle(paste0(cPaste(gene), " exon model"));
   }
   if ("df" %in% return_type) {
      return(grl1a1df);
   }
   grl1a1gg;
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
