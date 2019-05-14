
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
#'    include `grl_name` which are names of the input GRangesList
#'    `names(grl)`; `gr_name` which are names of the GRanges entries; and
#'    other columns from the input GRanges entries. When `shape="junction"`
#'    the data includes two polygons per junction, intended to be used
#'    with `geom_diagonal_wide_arc()` for each side in order to
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
#' @param scoreColname,scoreArcMinimum,scoreFactor,scoreArcFactor numeric
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
#' @family splicejam core functions
#'
#' @examples
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(jamba));
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
 widthV=c(exon=0.6, cds=0.6, noncds=0.3, intron=0.01, gap=0.01, `NA`=0.01),
 width_colname=c("subclass", "feature_type"),
 shape=c("rectangle", "junction"),
 baseline=NULL,
 scoreColname="score",
 scoreArcMinimum=100,
 scoreFactor=1,
 scoreArcFactor=0.5,
 doStackJunctions=TRUE,
 strandedScore=TRUE,
 ref2c=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert GRangesList to tall data.frame
   shape <- match.arg(shape);
   offset <- width / 2;
   if ("GRanges" %in% class(grl)) {
      if (verbose) {
         printDebug("grl2df(): ",
            "converting input to GRangesList");
      }
      grl <- GRangesList(list(gr=grl));
   }
   if (addGaps && !"junction" %in% shape) {
      if (verbose) {
         printDebug("grl2df(): ",
            "calling addGRLgaps()");
      }
      grl <- addGRLgaps(grl,
         verbose=verbose,
         ...);
   }
   ## TODO: handle gaps, perhaps by calling this function with
   ## result of getGRLgaps(), but setting addGaps=FALSE, and using
   ## width=0.05
   ## TODO: draw arrows on gap regions, optionally at the end of
   ## multi-rectangle features, to indicate strandedness.
   if ("rectangle" %in% shape) {
      if (verbose) {
         printDebug("grl2df(): ",
            "rectangle processing");
      }
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
         if (verbose) {
            printDebug("grl2df(): ",
               "width_colname_use:",
               width_colname_use);
         }
         wFound <- match(rmNA(naValue="NA",
            values(grl@unlistData)[[width_colname_use]]),
            names(widthV));
         width[!is.na(wFound)] <- widthV[wFound[!is.na(wFound)]];
      }
      offset <- rep(width, each=4) / 2;
      if (length(baseline) > 0) {
         if (length(baseline) == length(grl)) {
            yBaseline <- rep(
               rep(baseline,
                  elementNROWS(grl)),
               each=4);
         } else {
            yBaseline <- rep(
               rep(baseline,
                  length.out=length(grl@unlistData)),
               each=4);
         }
      } else {
         yBase <- rep(seq_along(grl) - 1,
            elementNROWS(grl));
         ySeqnames <- as.character(unlist(seqnames(grl)));
         abs_rank <- function(x){
            xu <- nameVector(unique(x));
            xr <- rank(xu, ties.method="min");
            unname(xr[as.character(x)] - 1);
         }
         yBase <- shrinkMatrix(yBase,
            groupBy=ySeqnames,
            shrinkFunc=abs_rank)$x;
         yName <- shrinkMatrix(seq_along(yBase),
            groupBy=ySeqnames,
            shrinkFunc=c)$x;
         yBaseline <- rep(
            yBase[yName],
            each=4);
      }
      yCoords <- rep(c(-1, -1, 1, 1) * offset,
         length.out=length(xCoords)) +
         yBaseline;
      id <- rep(seq_along(grl@unlistData),
         each=4);
      df <- data.frame(x=xCoords,
         y=yCoords,
         id=id,
         seqnames=rep(seqnames(grl@unlistData), each=4));
      if (length(names(grl)) > 0) {
         grlNames <- factor(
            rep(
               rep(names(grl),
                  elementNROWS(grl)),
               each=4),
            levels=names(grl));
         df$grl_name <- grlNames;
      }
      if (length(names(grl@unlistData)) > 0) {
         if (any(duplicated(names(grl@unlistData)))) {
            names(grl@unlistData) <- makeNames(rmNA(naValue="",
               names(grl@unlistData)));
         }
         grNames <- factor(
            rep(
               names(grl@unlistData),
               each=4),
            levels=names(grl@unlistData));
         df$gr_name <- grNames;
      }
      if (keepGRvalues) {
         for (iCol in colnames(values(grl@unlistData))) {
            df[,iCol] <- rep(values(grl@unlistData)[,iCol], each=4);
         }
      }
      if (keepGRLvalues) {
         for (iCol in colnames(values(grl))) {
            if (!iCol %in% colnames(df)) {
               df[,iCol] <- rep(
                  rep(values(grl)[,iCol],
                     elementNROWS(grl)),
                  each=4);
            }
         }
      }
   } else if ("junction" %in% shape) {
      #############################################################
      ## optionally stack junction ends so they do not overlap
      if (verbose) {
         printDebug("grl2df(): ",
            "junction processing");
      }
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
         if (length(baseline) == 0) {
            baseline <- 0;
         }
         grlNew <- GRangesList(
            lapply(grl, function(iGR){
               stackJunctions(gr=iGR,
                  scoreColname=scoreColname,
                  baseline=baseline,
                  strandedScore=strandedScore,
                  scoreFactor=1,
                  verbose=verbose,
                  ...);
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
         xMid <- (xStart + xEnd) / 2;
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
      ))+0.5;
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

      ## Calculate the height of each junction arc
      yScore <- values(grl@unlistData)[[scoreColname]] * 1;
      yMaxStartEnd <- pmax(abs(yStart), abs(yEnd));
      if (verbose) {
         printDebug("grl2df(): ",
            "calling internal_junc_score()");
         print(head(grl@unlistData));
      }
      internalMaxScore <- internal_junc_score(grl@unlistData,
         scoreColname=scoreColname,
         sampleColname="sample_id",
         verbose=verbose,
         ...);
      if (verbose) {
         printDebug("grl2df(): ",
            "internalMaxScore:", internalMaxScore);
      }

      #yHeight <- pmax(yMaxStartEnd, internalMaxScore) - internalMaxScore
      arcBaseHeight <- pmax(abs(yScore), scoreArcMinimum) * scoreArcFactor;
      internalArcBaseHeight <- abs(internalMaxScore) * (1.1+scoreArcFactor);
      #yBaseHeight <- noiseFloor(internalArcBaseHeight - yMaxStartEnd,
      #   minimum=scoreArcMinimum*scoreArcFactor);
      yBaseHeight <- noiseFloor(internalArcBaseHeight - yStart,
         minimum=scoreArcMinimum*scoreArcFactor);
      yHeight <- pmax(arcBaseHeight, yBaseHeight);

      yHeightOld1 <- (pmax(abs(yScore), scoreArcMinimum) + yMaxStartEnd) * scoreArcFactor;
      yHeightOld <- abs(yScore) * scoreArcFactor + scoreArcMinimum + yMaxStartEnd;
      yL <- list(yStart=yStart,
         yEnd=yEnd,
         yScore=yScore,
         yMaxStartEnd=yMaxStartEnd,
         yHeightOld=yHeightOld,
         yHeightOld1=yHeightOld1,
         yHeight=yHeight,
         internalMaxScore=internalMaxScore);
      #print(data.frame(yL));
      if (strandedScore) {
         yHeight <- yHeight * ifelse(
            as.character(GenomicRanges::strand(grl@unlistData)) %in% "-",
            -1, 1);
      }
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
         df$grl_name <- grlNames;
      }
      if (length(names(grl@unlistData)) > 0) {
         grNames <- factor(
            rep(
               names(grl@unlistData),
               each=8),
            levels=names(grl@unlistData));
         df$gr_name <- grNames;
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
#' `flattenExonsBy()`, and essentially plots the end result
#' for review.
#'
#' Alternatively, when `return_type="df"`, the output is
#' the `data.frame` used to produce the ggplot, which allows
#' for more customization.
#'
#' @family jam plot functions
#' @family jam ggplot2 functions
#' @family splicejam core functions
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
#'    `flattenExonsBy()`.
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
#' @param exonLabelSize numeric value, compatible with
#'    argument `size` in `ggrepel::geom_text_repel()`, used
#'    to size exon labels when `labelExons=TRUE`.
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
#' @param tx2geneDF data.frame or NULL, optionally used to help
#'    identify matching transcripts for the requested `gene` value,
#'    used when `"gene_name"` is not present in `values(flatExonsByTx)`.
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
#' }
#'
#' @export
gene2gg <- function
(gene=NULL,
 tx=NULL,
 flatExonsByGene=NULL,
 flatExonsByTx=NULL,
 geneColor="dodgerblue",
 labelExons=TRUE,
 exonLabelAngle=90,
 exonLabelSize=3,
 newValues=list(feature_type="gap", subclass="gap", gene_nameExon="gap"),
 gene_order=c("first","last"),
 return_type=c("grob", "df"),
 ref2c=NULL,
 hjust=1.2,
 vjust=0,
 compressGaps=TRUE,
 tx2geneDF=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is a lightweight wrapper around grl2df() specifically intended
   ## for gene and transcript exon structure
   if (!suppressPackageStartupMessages(require(ggplot2))) {
      stop("gene2gg() requires ggplot2.");
   }
   if (!suppressPackageStartupMessages(require(ggforce))) {
      stop("gene2gg() requires ggforce.");
   }
   gene_order <- match.arg(gene_order);
   return_type <- match.arg(return_type);
   if (length(flatExonsByGene) > 0 && length(gene) > 0) {
      if (!any(gene %in% names(flatExonsByGene))) {
         stop("gene was not found in names(flatExonsByGene)");
      }
      if (verbose) {
         printDebug("gene2gg(): ",
            "flatExonsByGene[gene]");
      }
      grl1a <- flatExonsByGene[names(flatExonsByGene) %in% gene];
   } else {
      grl1a <- NULL;
   }
   if (length(flatExonsByTx) > 0 && igrepHas("GRanges", class(flatExonsByTx))) {
      grl1 <- NULL;
      if (length(gene) > 0) {
         if ("gene_name" %in% colnames(values(flatExonsByTx))) {
            if (verbose) {
               printDebug("gene2gg(): ",
                  "values(flatExonsByTx)$gene_name %in% gene");
            }
            grl1 <- subset(flatExonsByTx, gene_name %in% gene);
         } else if (length(tx2geneDF) > 0 &&
            "gene_name" %in% colnames(tx2geneDF)) {
            if (verbose) {
               printDebug("gene2gg(): ",
                  "subset(tx2geneDF, gene_name %in% gene)$transcript_id");
            }
            tx <- unique(c(tx,
               subset(tx2geneDF, gene_name %in% gene)$transcript_id));
         }
      }
      if (length(tx) > 0) {
         grl1 <- GRangesList(c(grl1,
            flatExonsByTx[names(flatExonsByTx) %in% tx]))
         values(grl1)$gene_name <- tx2geneDF[match(names(grl1),
            tx2geneDF$transcript_id),"gene_name"];
         values(grl1@unlistData)$gene_name <- rep(values(grl1)$gene_name,
            elementNROWS(grl1));
         values(grl1)$transcript_id <- names(grl1);
         values(grl1@unlistData)$transcript_id <- rep(values(grl1)$transcript_id,
            elementNROWS(grl1));
      }
   } else {
      grl1 <- NULL;
   }
   if ("first" %in% gene_order) {
      grl1a1 <- GRangesList(c(grl1, grl1a));
   } else {
      grl1a1 <- GRangesList(c(grl1a, grl1));
   }
   if (verbose) {
      printDebug("class(grl1a1):", class(grl1a1));
   }
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
         if (verbose) {
            printDebug("gene2gg(): ",
               "grl1a1:");
            print(grl1a1);
         }
         ref2c <- make_ref2compressed(grl1a1@unlistData,
            ...);
      } else {
         ref2c <- NULL;
      }
   }
   ## Calculate a reasonable y-axis minimum to allow for exon labels,
   ## based upon the number of transcripts being displayed
   ymin <- (-0.1 +
         -1 * (exonLabelSize/10) +
         -1 * ((labelExons*1) * length(grl1a1))/6);
   ## Put it together
   grl1a1gg <- ggplot2::ggplot(grl1a1df,
         aes(x=x,
            y=y,
            group=id)) +
      ggforce::geom_shape(show.legend=FALSE,
         aes(fill=subclass,
            text=paste0(gr_name,
               "<br>", grl_name,
               "<br>",subclass,
               "<br>coord:", scales::comma(x)
            ),
            color=subclass)) +
      colorjam::theme_jam() +
      ylab("") +
      scale_color_manual(values=makeColorDarker(colorSubV,
         darkFactor=1.3)) +
      scale_fill_manual(values=alpha2col(colorSubV,
         alpha=1)) +
      scale_y_continuous(breaks=seq_along(grl1a1)-1,
         limits=c(ymin, length(grl1a1)-0.5),
         labels=names(grl1a1));
   if (length(ref2c) > 0) {
      grl1a1gg <- grl1a1gg +
         scale_x_continuous(trans=ref2c$trans_grc);
   }
   if (labelExons) {
      if (1 == 1 || length(ref2c) == 0) {
         grl1a1gg <- grl1a1gg +
            ggrepel::geom_text_repel(
               inherit.aes=FALSE,
               data=exonLabelDF,
               aes(x=x,
                  y=min(y),
                  text=NULL,
                  label=gene_nameExon),
               angle=exonLabelAngle,
               hjust=vjust,
               vjust=hjust,
               segment.color="grey35",
               #fill="white",
               size=exonLabelSize,
               direction="x");
      } else {
         grl1a1gg <- grl1a1gg +
            geom_text(
               inherit.aes=FALSE,
               data=exonLabelDF,
               aes(x=x,
                  y=min(y),
                  text=NULL,
                  label=gene_nameExon),
               angle=exonLabelAngle,
               hjust=hjust,
               vjust=vjust,
               #fill="white",
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
#' This function is intended to help visualize splice junctions
#' specifically when plotted using `geom_diagonal_wide_arc()`,
#' where the height of the junction arc is defined by the `score`.
#' When two junctions have the same start position, their y-positions
#' are stacked, such that the shorter junction width is placed before
#' longer junction widths. The intention is to reduce visible overlaps.
#'
#' The input data is expected to have annotations similar to
#' those provided by `closestExonToJunctions()`, specifically
#' the columns `"nameFrom"` and `"nameTo"`, see examples below.
#'
#' @family jam plot functions
#' @family jam GRanges functions
#'
#' @return GRanges with colnames `c("yStart", "yEnd")` added
#'    to `values(gr)`, indicating the baseline y-axis position
#'    for the start and end of the junction arc. The score
#'    `values(gr)[[scoreColname]]` will reflect the adjustments
#'    by `scoreFactor`, and if `strandedScore=TRUE` then all
#'    strand `"-"` scores will be negative, all other scores
#'    will be positive.
#'
#' @param gr GRanges object representing splice junctions.
#' @param scoreColname character string matching one of `colnames(values(gr))`
#'    that contains a numeric value representing the abundance of
#'    each splice junction observed.
#' @param scoreFactor numeric value multiplied by the value in `scoreColname`
#'    to allow scaling the junctions across samples. Note that
#'    `scoreFactor` can be a vector, which would be applied to the
#'    vector of scores.
#' @param matchFrom,matchTo optional colnames to use when grouping
#'    junctions at the start and end positions. By default `"nameFrom"`
#'    and `"nameTo"` are used, as output from `closestExonToJunctions()`,
#'    which has the benefit of grouping junctions within the
#'    `spliceBuffer` distance from exon boundaries. If those values
#'    are not present `colnames(values(gr))` then the new
#'    default `c("seqnames", "start", "strand")` is used for
#'    `matchFrom`, and `c("seqnames", "end", "strand")` is used
#'    for `matchTo`. That said, if `matchFrom` or `matchTo` are supplied,
#'    those colnames are used from `as.data.frame(gr)`. Multiple colnames
#'    are allowed.
#' @param strandedScore logical indicating whether to enforce negative
#'    scores for junctions on the `"-"` strand. Note that when `strandedScore`
#'    is true, all `"-"` strand scores will be negative, and all other
#'    scores with be positive.
#' @param baseline numeric vector of length 0, 1 or `length(gr)`, with values
#'    added to the y-axis value for junctions.
#'    If `baseline` has names matching `names(gr)` they will be used for
#'    each `gr` entry; if `baseline` is not named, values are recycled
#'    to `length(gr)`. The purpose is to allow exons to be shifted up
#'    or down on the y-axis, along with associated junctions and
#'    coverage data (see `exoncov2polygon()` for another example.)
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' library(GenomicRanges);
#' library(ggplot2);
#' library(ggforce);
#' library(colorjam);
#' library(jamba);
#' grExons <- GRanges(seqnames=rep("chr1", 4),
#'    ranges=IRanges(
#'       start=c(100, 300, 500,  900),
#'         end=c(200, 400, 750, 1000)),
#'    strand=rep("+", 4));
#' names(grExons) <- makeNames(rep("exon", length(grExons)),
#'    suffix="");
#'
#' grJunc <- GRanges(seqnames=rep("chr1", 5),
#'    ranges=IRanges(start=c(200, 200, 400, 400, 750),
#'       end=c(300, 500, 500, 900, 900)),
#'    strand=rep("+", 5),
#'    score=c(200, 50, 160, 40, 210));
#' names(grJunc) <- makeNames(rep("junc", length(grJunc)));
#'
#' # quick plot showing exons and junctions using rectangles
#' grl <- c(
#'    GRangesList(exons=grExons),
#'    split(grJunc, names(grJunc))
#'    );
#' ggplot(grl2df(grl), aes(x=x, y=y, group=id, fill=feature_type)) +
#'    ggforce::geom_shape() +
#'    scale_y_continuous(breaks=seq_along(grl)-1, labels=names(grl)) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    ggtitle("Schematic of exons and junctions GRanges");
#'
#' # add annotation for closest known exon
#' grJunc <- closestExonToJunctions(grJunc, grExons)$spliceGRgene;
#'
#' # The un-stacked junctions
#' grlJunc2df1 <- grl2df(grJunc,
#'    shape="junction",
#'    doStackJunctions=FALSE);
#' ggplot(grlJunc2df1, aes(x=x, y=y, group=gr_name, fill=gr_name)) +
#'    geom_diagonal_wide_arc(alpha=0.7) +
#'    colorjam::scale_fill_jam() +
#'    colorjam::theme_jam() +
#'    ggtitle("Junctions not stacked at boundaries")
#'
#' # The stacked junctions
#' grJunc2 <- stackJunctions(grJunc);
#' grlJunc2df2 <- grl2df(grJunc2,
#'    scoreArcMinimum=20,
#'    shape="junction");
#' ggplot(grlJunc2df2, aes(x=x, y=y, group=gr_name, fill=gr_name)) +
#'    geom_diagonal_wide_arc(alpha=0.7) +
#'    colorjam::scale_fill_jam() +
#'    colorjam::theme_jam() +
#'    ggtitle("Junctions stacked at boundaries")
#'
#'
#' @export
stackJunctions <- function
(gr,
 scoreColname="score",
 scoreFactor=1,
 matchFrom=NULL,
 matchTo=NULL,
 strandedScore=TRUE,
 baseline=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to stack junctions by score (width) so they do
   ## not overlap at the start or end of each junction.
   if (!scoreColname %in% colnames(values(gr))) {
      stop("The scoreColname must be present in colnames(values(gr))");
   }

   ## Validate matchFrom, matchTo
   if (length(matchFrom) == 0) {
      if ("nameFrom" %in% colnames(values(gr))) {
         matchFrom <- "nameFrom";
      } else {
         matchFrom <- c("seqnames", "start");
      }
   }
   if (length(matchTo) == 0) {
      if ("nameTo" %in% colnames(values(gr))) {
         matchTo <- "nameTo";
      } else {
         matchTo <- c("seqnames", "end");
      }
   }
   ##
   grValueColnames <- colnames(values(gr));
   grValues <- c("seqnames", "start", "end", "width", "strand",
      grValueColnames);
   matchFrom <- intersect(matchFrom, grValues);
   if (length(matchFrom) == 0) {
      stop("matchFrom must match colnames(as.data.frame(gr))");
   }
   matchTo <- intersect(matchTo, grValues);
   if (length(matchTo) == 0) {
      stop("matchTo must match colnames(as.data.frame(gr))");
   }
   if (verbose) {
      printDebug("stackJunctions(): ",
         "matchFrom:", matchFrom,
         ", matchTo:", matchTo);
   }
   ## Optionally enforce strandedness
   if (strandedScore) {
      if (verbose) {
         printDebug("stackJunctions(): ",
            "Adjusting stranded score.");
      }
      matchFrom <- c(matchFrom, "strand");
      matchTo <- c(matchTo, "strand");
      scoreV <- scoreFactor *
         ifelse(as.character(GenomicRanges::strand(gr)) %in% "-", -1, 1) *
         abs(values(gr)[[scoreColname]]);
   } else {
      if (verbose) {
         printDebug("stackJunctions(): ",
            "Adjusting unstranded score.");
      }
      scoreV <- scoreFactor *
         values(gr)[[scoreColname]];
   }
   values(gr)[[scoreColname]] <- scoreV;
   # Combine value colnames and non-value colnames, using only
   # the required colnames since as.data.frame(gr) will fail if
   # any columns are not coercible to data.frame, for example list.
   exonsFrom <- pasteByRow(sep="_",
      as.data.frame(gr[,intersect(matchFrom, grValueColnames)])[,matchFrom,drop=FALSE])
   exonsTo <- pasteByRow(sep="_",
      as.data.frame(gr[,intersect(matchTo, grValueColnames)])[,matchTo,drop=FALSE])

   exonsFrom <- gsub("[.][-]*[0-9]+$", "",
      exonsFrom);
   exonsTo <- gsub("[.][-]*[0-9]+$", "",
      exonsTo);
   if (verbose) {
      printDebug("stackJunctions(): ",
         "head(exonsFrom):",
         head(exonsFrom, 10));
      printDebug("stackJunctions(): ",
         "head(exonsTo):",
         head(exonsTo, 10));
   }

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
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
   ## group by exonsFrom which allows the spliceBuffer to work
   order1 <- do.call(order, list(exonsFrom, width(gr)));
   values(gr)[order1,"yStart"] <- shrinkMatrix(
      scoreV[order1],
      #groupBy=start(gr[order1]),
      groupBy=exonsFrom[order1],
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))})$x +
      baselineV[exonsFrom[order1]];

   ## End position
   ## group by exonsTo which allows the spliceBuffer to work
   order2 <- do.call(order, list(exonsTo, width(gr)));
   values(gr)[order2,"yEnd"] <- shrinkMatrix(
      scoreV[order2],
      #groupBy=end(gr[order2]),
      groupBy=exonsTo[order2],
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))})$x +
      baselineV[exonsTo[order2]];

   ## Bonus points
   ## rank junctions by score at the exonsFrom and exonsTo position
   ## so the rank can be used to colorize dominant junctions
   shrink_junc <- function(i){
      if (jamba::igrepHas("numeric|integer", class(i))) {
         (rank(-abs(i)) == 1) * 1;
      } else {
         i;
      }
   }
   juncRankFrom <- shrinkMatrix(
      values(gr)[,c("score","nameFromToSample")],
      pasteByRow(values(gr)[,c("sample_id","nameFrom")]),
      shrinkFunc=shrink_junc);
   juncRankTo <- shrinkMatrix(
      values(gr)[,c("score","nameFromToSample")],
      pasteByRow(values(gr)[,c("sample_id","nameTo")]),
      #shrinkFunc=c)
      shrinkFunc=shrink_junc);
   juncRankDF <- data.frame(
         juncRankFrom[match(values(gr)[,"nameFromToSample"], juncRankFrom[,"nameFromToSample"]),],
         juncRankTo[match(values(gr)[,"nameFromToSample"], juncRankTo[,"nameFromToSample"]),]
   );
   juncRank <- (1 +
         juncRankFrom[match(values(gr)[,"nameFromToSample"], juncRankFrom[,"nameFromToSample"]),"score"] +
         juncRankTo[match(values(gr)[,"nameFromToSample"], juncRankTo[,"nameFromToSample"]),"score"]
   );
   values(gr)[,"junction_rank"] <- factor(juncRank);
   return(gr);
}

#' Jam Sashimi plot
#'
#' Jam Sashimi plot
#'
#' This function uses Sashimi data prepared by `prepareSashimi()`
#' and creates a ggplot graphical object ready for visualization.
#' As a result, this function provides several arguments to
#' customize the visualization.
#'
#' @family jam plot functions
#' @family jam ggplot2 functions
#' @family splicejam core functions
#'
#' @param sashimi Sashimi data prepared by `prepareSashimi()` which
#'    is a `list` with `covDF` coverage data in data.frame format,
#'    `juncDF` junction data in data.frame format,
#'    `juncLabelDF` junction label coordinates in data.frame format,
#'    `exonLabelDF` exon label coordinates per coverage polygon in
#'    data.frame format,
#'    `ref2c` list output from `make_ref2compressed()` to
#'    transform genomic coordinates.
#' @param show character vector of Sashimi plot features to include:
#'    `"coverage"` sequence read coverage data;
#'    `"junction"` splice junction read data.
#' @param coord_method character value indicating the type of
#'    coordinate scaling to use:
#'    `"scale"` uses `ggplot2::scale_x_continuous()`;
#'    `"coord"` uses `ggplot2::coord_trans()`;
#'    `"none"` does not compress genomic coordinates.
#' @param exonsGrl GRangesList object with one or more gene or
#'    transcript exon models, where exons are disjoint (not
#'    overlapping.)
#' @param junc_color,junc_fill character string with valid R color,
#'    used for junction outline, and fill, for the junction arc
#'    polygon. Alpha transparency is recommended for `junc_fill`
#'    so overlapping junction arcs are visible.
#' @param fill_scheme character string for how the exon coverages
#'    will be color-filled: `"exon"` will define colors for each
#'    distinct exon, using the GRanges names from `flatExonsByGene`;
#'    `"sample_id"` to color all exons the same by sample_id.
#' @param color_sub optional character vector of R compatible colors
#'    or hex strings, whose names are used to color or fill features
#'    in the ggplot object. For example, if `fill_sheme="sample_id"`
#'    the `color_sub` should have names for each `"sample_id"` value.
#'    If any values are missing, they will be filled in using
#'    `colorjam::rainbowJam()`.
#' @param use_jam_themes logical indicating whether to apply
#'    `colorjam::theme_jam()`, by default for the ggplot theme.
#' @param apply_facet logical indicating whether to apply
#'    `ggplot2::facet_wrap()` with `"~sample_id"` defining each
#'    panel.
#' @param facet_scales character value used as `"scales"` argument in
#'    `ggplot2::facet_wrap()` when `apply_facet=TRUE`.
#' @param ref2c optional output from `make_ref2compressed()` used to
#'    compress axis coordinates during junction arc calculations.
#' @param verbose logical indicating whether to print verbose output.
#' @param ... additional arguments are sent to `grl2df()`.
#'
#' @examples
#' suppressPackageStartupMessages(library(GenomicRanges));
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
plotSashimi <- function
(sashimi,
 show=c("coverage", "junction",
    "junctionLabels"),
 coord_method=c("scale", "coord", "none"),
 exonsGrl=NULL,
 junc_color=alpha2col("goldenrod2", 0.3),
 junc_fill=alpha2col("goldenrod2", 0.9),
 fill_scheme=c("sample_id", "exon"),
 color_sub=NULL,
 use_jam_themes=TRUE,
 apply_facet=TRUE,
 facet_scales="free_y",
 ref2c=NULL,
 do_highlight=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take prepared Sashimi data, and return ggplot
   coord_method <- match.arg(coord_method);
   fill_scheme <- match.arg(fill_scheme);
   if (length(ref2c) > 0) {
      sashimi$ref2c <- ref2c;
   }
   if (!"ref2c" %in% names(sashimi)) {
      coord_method <- "none";
   }
   ggCov <- NULL;
   #color_sub <- NULL;
   if ("coverage" %in% show && "covDF" %in% names(sashimi)) {
      #if (do_highlight) {
      #   covDF <- highlight_key(sashimi$covDF, key=~gr);
      #}
      # Prepare data.frame to be merged with juncDF
      covDF <- sashimi$covDF;
      covDF$type <- "coverage";
      covDF$name <- pasteByRow(covDF[,c("gr", "cov", "sample_id")], sep=" ");
      covDF$feature <- covDF$gr;
      covDF$row <- seq_len(nrow(covDF));
      if ("exon" %in% fill_scheme) {
         covDF$color_by <- covDF$gr;
      } else {
         covDF$color_by <- covDF$sample_id;
      }
      # Define colors
      color_sub_cov <- colorjam::group2colors(unique(covDF$color_by));
      use_names <- setdiff(names(color_sub_cov), names(color_sub));
      color_sub[use_names] <- color_sub_cov[use_names];

      if ("exonLabels" %in% show && "exonLabelDF" %in% names(sashimi)) {
         exonLabelDF <- sashimi$exonLabelDF;
         exonLabelDF$type <- "exon_label";
         exonLabelDF$name <- pasteByRow(exonLabelDF[,c("gr","sample_id")], sep=" ");
         exonLabelDF$feature <- exonLabelDF$gr;
         exonLabelDF$row <- seq_len(nrow(exonLabelDF));
         exonLabelDF$color_by <- NA;

      }
   }
   ## Junction data
   if ("junction" %in% show && "juncDF" %in% names(sashimi)) {
      color_sub_junc <- NULL;
      juncDF <- sashimi$juncDF;
      if (!"junction_rank" %in% colnames(juncDF)) {
         juncDF$junction_rank <- 3;
      }
      juncDF$type <- "junction";
      juncDF$name <- pasteByRow(juncDF[,c("nameFromTo", "sample_id")], sep=" ");
      juncDF$feature <- juncDF$nameFromTo;
      juncDF$row <- seq_len(nrow(juncDF));
      if (fill_scheme %in% "sample_id") {
         juncDF$color_by <- pasteByRow(juncDF[,c("sample_id","junction_rank")], sep=".");
      } else {
         #juncDF$color_by <- juncDF$nameFromToSample;
         juncDF$color_by <- as.character(juncDF$junction_rank);
      }

      # order so the lower rank, lesser junctions, are drawn on top
      juncDF <- juncDF[order(-juncDF$junction_rank),,drop=FALSE];

      ## gradually transparent white to shade by junction_rank
      junc_blank_3 <- nameVector(
         alpha2col(rep("#DDDDDD", 3), alpha=c(0.3, 0.67, 0.8)),
         c(3, 2, 1));
      if (fill_scheme %in% "sample_id") {
         if (all(unique(juncDF$sample_id) %in% names(color_sub))) {
            color_sub_samples <- color_sub[intersect(juncDF$sample_id, names(color_sub))];
         } else {
            color_sub_samples <- group2colors(unique(juncDF$sample_id));
         }
         # blend with gradually transparent white
         color_sub_junc_3 <- colorspace::hex(
            colorspace::mixcolor(
               hex2RGB(rgb2col(col2rgb(rep(color_sub_samples, each=3)))),
               hex2RGB(junc_blank_3),
               alpha=col2alpha(junc_blank_3)
            )
         );
         names(color_sub_junc_3) <- paste(names(color_sub_junc_3),
            names(junc_blank_3),
            sep=".");
      } else {
         # Otherwise blend junc_fill with white
         color_sub_junc_3 <- colorspace::hex(
            colorspace::mixcolor(
               hex2RGB(rgb2col(col2rgb(rep(junc_fill[1], each=3)))),
               hex2RGB(junc_blank_3),
               alpha=col2alpha(junc_blank_3)
            )
         );
         names(color_sub_junc_3) <- names(junc_blank_3);
      }
      use_names <- setdiff(names(color_sub_junc_3), names(color_sub));
      color_sub[use_names] <- color_sub_junc_3[use_names];
      junc_alpha <- mean(col2alpha(junc_fill));

      if ("junctionLabels" %in% show && "juncLabelDF" %in% names(sashimi)) {
         juncLabelDF <- sashimi$juncLabelDF;
         juncLabelDF$type <- "junction_label";
         juncLabelDF$name <- pasteByRow(juncLabelDF[,c("nameFromTo", "sample_id")], sep=" ");
         juncLabelDF$feature <- juncLabelDF$nameFromTo;
         juncLabelDF$row <- seq_len(nrow(juncLabelDF));
         if (!"junction_rank" %in% colnames(juncLabelDF)) {
            juncLabelDF$junction_rank <- "3";
         }
         if (fill_scheme %in% "sample_id") {
            juncLabelDF$color_by <- pasteByRow(juncLabelDF[,c("sample_id","junction_rank")], sep=".");
         } else {
            #juncDF$color_by <- juncDF$nameFromToSample;
            juncLabelDF$color_by <- as.character(juncLabelDF$junction_rank);
         }
         juncLabelDF$label <- scales::comma(round(juncLabelDF$score));
      }
   }

   # Assemble one data.frame in order to keep ggplot2 data in sync
   cjL <- list();
   if ("coverage" %in% show) {
      covDF$text <- paste0(
         "score:", scales::comma(covDF$y),
         "<br>coord:", scales::comma(covDF$x),
         "<br>feature:", as.character(covDF$gr),
         "<br>sample_id:", covDF$sample_id,
         "<br>track:", as.character(covDF$cov));
      cjL$coverage <- covDF;
   }
   if ("junction" %in% show) {
      juncDF$text <- paste0(
         "score:", scales::comma(juncDF$score),
         "<br>nameFrom:", juncDF$nameFrom,
         "<br>nameTo:", juncDF$nameTo,
         "<br>sample_id:", juncDF$sample_id
      );
      cjL$junction <- juncDF;
   }
   if ("junctionLabels" %in% show) {
      cjL$junction_label <- juncLabelDF;
   }
   if (length(cjL) == 1) {
      cjDF <- cjL[[1]];
   } else {
      cjDF <- jamba::mergeAllXY(cjL);
   }
   cjDF <- mixedSortDF(cjDF,
      byCols=c("type","row"));
   # order columns by presence of NA values
   na_ct <- apply(cjDF, 2, function(i){
      sum(is.na(i))
   });
   cjDF <- cjDF[,order(na_ct),drop=FALSE];

   # Create ggplot2 piece by piece
   if (do_highlight) {
      cjDFh <- highlight_key(cjDF,
         key=~feature);
      gg_sashimi <- ggplot(
         data=cjDFh,
         aes(
            x=x,
            y=y,
            group=name,
            color=color_by,
            fill=color_by
         ));
   } else {
      gg_sashimi <- ggplot(
         data=cjDF,
         aes(
            x=x,
            y=y,
            group=name,
            color=color_by,
            fill=color_by
         ));
   }

   # Add coverage layer
   color_sub_d <- makeColorDarker(color_sub, darkFactor=1.2);
   if (all(c("coverage","junction") %in% cjDF$type)) {
      gg_sashimi <- gg_sashimi +
         geom_polygon(
            data=. %>% filter(type %in% "coverage"),
            show.legend=FALSE) +
         geom_diagonal_wide_arc(
            data=. %>% filter(type %in% "junction"),
            show.legend=FALSE);
   } else if ("coverage" %in% cjDF$type) {
      gg_sashimi <- gg_sashimi +
         ggforce::geom_shape(
            data=. %>% filter(type %in% "coverage"),
            show.legend=FALSE
         );
   } else if ("junction" %in% cjDF$type) {
      gg_sashimi <- gg_sashimi +
         geom_diagonal_wide_arc(
            data=. %>% filter(type %in% "junction"),
            show.legend=FALSE,
            alpha=junc_alpha,
            strength=0.4
         );
   }

   # Add junction_labels
   if ("junction_label" %in% cjDF$type) {
      gg_sashimi <- gg_sashimi +
         ggrepel::geom_text_repel(
            data=. %>% filter(type %in% "junction_label"),
            angle=90,
            vjust=0.5,
            direction="y",
            point.padding=0,
            color="black",
            fill="transparent",
            aes(
               text=NULL,
               label=label
            )
         );
   }

   if ("scale" %in% coord_method) {
      gg_sashimi <- gg_sashimi +
         scale_x_continuous(trans=sashimi$ref2c$trans_grc);
   } else if ("coord" %in% coord_method) {
      gg_sashimi <- gg_sashimi +
         coord_trans(x=sashimi$ref2c$trans_grc);
   }
   if (use_jam_themes) {
      gg_sashimi <- gg_sashimi +
         colorjam::theme_jam(
            panel.grid.major.colour="grey80",
            panel.grid.minor.colour="grey90");
   }
   # Apply colors
   gg_sashimi <- gg_sashimi +
      scale_fill_manual(values=color_sub) +
      scale_color_manual(values=makeColorDarker(color_sub, darkFactor=1.2));

   # Apply facet
   if (apply_facet) {
      gg_sashimi <- gg_sashimi +
         facet_grid(sample_id~.,
         scales=facet_scales);
   }
   return(gg_sashimi);
}

#' Support plotly for GeomShape
#'
#' Support plotly for GeomShape
#'
#' This function is a helper function intended to provide a simple
#' but not fully-equivalent link between `ggforce::geom_shape()` and
#' `plotly::ggplotly()`. This function converts data to
#' `ggplot2::geom_polygon` which is sufficient for the splicejam
#' package, but which does not provide the extra effects from
#' `ggforce::geom_shape()` such as rounded corners or resizing.
#' This function may be extended then submitted to
#' the `ggforce` package to assist plotly support.
#'
#' This function is not exported, as it is only used by plotly
#' specifically during conversion of a ggplot object to plotly object.
#'
#' @param data,prestats_data,layout,params,p,... arguments provided
#'    by plotly during rendering, after stats have been applied.
#'    Currently only `data` is passed to `ggplot2::GeomPolygon`.
#'
to_basic.GeomShape <- function
(data, prestats_data, layout, params, p, ...)
{
   plotly:::prefix_class(data, "GeomPolygon");
}

