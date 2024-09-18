
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
#' @returns `data.frame` with `x,y` coordinates, and `id` which is used
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
#' @param grl `GRangesList`, or `GRanges` which will be converted to a
#'    GRangesList of length=1.
#' @param keepGRvalues,keepGRLvalues `logical` indicating whether the
#'    output data.frame should include column values from GRangesList and
#'    GRanges, if available.
#' @param addGaps `logical` indicating whether to add gap GRanges between
#'    same-strand GRanges features within each GRangesList element.
#'    When `TRUE` the gaps will be drawn between each GRanges rectangle.
#' @param width `numeric` value of default width for features when
#'    `shape="rectangle"`.
#' @param widthV `numeric` vector whose names are column values, using values
#'    from the first available colname from `width_colname`. Some common
#'    defaults are provided. Values are suggested to be between 0 and 1,
#'    since each GRangesList element is separated by 1 y-axis unit.
#' @param width_colname when `widthV` is used to determine width based
#'    upon a column value, the first matching colname of `width_colname`
#'    is used.
#' @param shape `character` string indicating whether input data should
#'    be processed as rectangular segments or splice junction arcs.
#' @param scoreColname,scoreArcMinimum,scoreFactor,scoreArcFactor `numeric`
#'    values used to determine junction ribbon height, the minimum
#'    height of the arc above the starting y-axis values based upon
#'    the score, the scaling factor for score values, and the
#'    relative height of the arc above the starting y-axis values
#'    multiplied by the score.
#' @param sampleColname `character` string indicating the column
#'    containing biological sample identifier. This column is only
#'    used when `type="junction"`, and when `doStackJunctions=TRUE`.
#'    It is used to ensure the junctions are only stacked within
#'    each sample. When `sampleColname` is not present in
#'    `colnames(GenomicRanges::values(grl@unlistData))` then all junctions are
#'    stacked.
#' @param doStackJunctions `logical` indicating whether to stack junctions
#'    at the start and end of junctions sharing the same coordinate,
#'    in order of shortest to longest junction width.
#' @param ref2c optional output from `make_ref2compressed()` used to
#'    compress axis coordinates during junction arc calculations.
#' @param verbose `logical` indicating whether to print verbose output.
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
#' gr <- GenomicRanges::GRanges(seqnames=rep(c("chr1"), 7),
#'    ranges=IRanges::IRanges(start=c(50, 100, 1300, 2500, 23750, 24900, 25000),
#'       end=c(100, 150, 1450, 2600, 23800, 25000, 25200)),
#'    strand=rep("+", 7),
#'    feature_type=rep(c("noncds", "cds", "noncds"), c(1,5,1)));
#' names(gr) <- jamba::makeNames(rep("exon", 7));
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
#'    ggplot2::coord_trans(x=ref2c$trans_grc) + colorjam::theme_jam();
#' print(gg3);
#'
#' ## An example showing splice junction data
#' data(test_junc_wide_gr);
#' junc_wide_df <- grl2df(test_junc_wide_gr, shape="junction");
#' ggWide1 <- ggplot2::ggplot(junc_wide_df,
#'    ggplot2::aes(x=x, y=y, group=gr_name, fill=gr_name, color=gr_name)) +
#'   splicejam::geom_diagonal_wide_arc() +
#'   colorjam::theme_jam() +
#'   colorjam::scale_fill_jam(alpha=0.7) +
#'   colorjam::scale_color_jam() +
#'   ggplot2::xlab("chr1") +
#'   ggplot2::ggtitle("junctions (full intron width)")
#' print(ggWide1);
#'
#' @export
grl2df <- function
(grl,
 keepGRvalues=TRUE,
 keepGRLvalues=FALSE,
 addGaps=TRUE,
 width=0.6,
 widthV=c(exon=0.6, cds=0.6, noncds=0.3, intron=0.01, gap=0.01, `NA`=0.5),
 width_colname=c("subclass", "feature_type"),
 shape=c("rectangle", "junction"),
 baseline=NULL,
 scoreColname="score",
 sampleColname="sample_id",
 scoreArcMinimum=200,
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
         jamba::printDebug("grl2df(): ",
            "converting input to GRangesList");
      }
      grl <- GenomicRanges::GRangesList(list(gr=grl));
   }
   if (addGaps && !"junction" %in% shape) {
      if (verbose) {
         jamba::printDebug("grl2df(): ",
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
         jamba::printDebug("grl2df(): ",
            "rectangle processing");
      }
      xCoords <- as.vector(rbind(
         GenomicRanges::start(grl@unlistData),
         GenomicRanges::end(grl@unlistData),
         GenomicRanges::end(grl@unlistData),
         GenomicRanges::start(grl@unlistData)
      ));
      width <- rep(width, length.out=length(grl@unlistData));
      width_colname_use <- head(
         intersect(width_colname,
            colnames(GenomicRanges::values(grl@unlistData))),
         1);
      if (length(width_colname_use) > 0) {
         if (verbose) {
            jamba::printDebug("grl2df(): ",
               "width_colname_use:",
               width_colname_use);
         }
         wFound <- match(jamba::rmNA(naValue="NA",
            GenomicRanges::values(grl@unlistData)[[width_colname_use]]),
            names(widthV));
         width[!is.na(wFound)] <- widthV[wFound[!is.na(wFound)]];
      }
      offset <- rep(width, each=4) / 2;
      if (length(baseline) > 0) {
         if (length(baseline) == length(grl)) {
            yBaseline <- rep(
               rep(baseline,
                  S4Vectors::elementNROWS(grl)),
               each=4);
         } else {
            yBaseline <- rep(
               rep(baseline,
                  length.out=length(grl@unlistData)),
               each=4);
         }
      } else {
         yBase <- rep(seq_along(grl) - 1,
            S4Vectors::elementNROWS(grl));
         ySeqnames <- as.character(unlist(GenomicRanges::seqnames(grl)));
         abs_rank <- function(x){
            xu <- jamba::nameVector(unique(x));
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
         seqnames=rep(GenomicRanges::seqnames(grl@unlistData), each=4));
      if (length(names(grl)) > 0) {
         grlNames <- factor(
            rep(
               rep(names(grl),
                  S4Vectors::elementNROWS(grl)),
               each=4),
            levels=names(grl));
         df$grl_name <- grlNames;
      }
      if (length(names(grl@unlistData)) > 0) {
         if (any(duplicated(names(grl@unlistData)))) {
            names(grl@unlistData) <- jamba::makeNames(jamba::rmNA(naValue="",
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
         for (iCol in colnames(GenomicRanges::values(grl@unlistData))) {
            df[,iCol] <- rep(GenomicRanges::values(grl@unlistData)[,iCol], each=4);
         }
      }
      if (keepGRLvalues) {
         for (iCol in colnames(GenomicRanges::values(grl))) {
            if (!iCol %in% colnames(df)) {
               df[,iCol] <- rep(
                  rep(GenomicRanges::values(grl)[,iCol],
                     S4Vectors::elementNROWS(grl)),
                  each=4);
            }
         }
      }
   } else if ("junction" %in% shape) {
      #############################################################
      ## optionally stack junction ends so they do not overlap
      if (verbose) {
         jamba::printDebug("grl2df(): ",
            "junction processing");
      }
      if (!all(scoreFactor == 1)) {
         if (verbose) {
            jamba::printDebug("grl2df(): ",
               "scoreFactor:",
               scoreFactor);
         }
         GenomicRanges::values(grl@unlistData)[[scoreColname]] <- (scoreFactor *
               GenomicRanges::values(grl@unlistData)[[scoreColname]])
      }
      if (doStackJunctions) {
         if (verbose) {
            jamba::printDebug("grl2df(): ",
               "stackJunctions()");
         }
         if (length(baseline) == 0) {
            baseline <- 0;
         }
         sampleColname <- intersect(sampleColname,
            colnames(GenomicRanges::values(grl@unlistData)));
         grlNew <- GenomicRanges::GRangesList(
            lapply(grl, function(iGR){
               stackJunctions(gr=iGR,
                  scoreColname=scoreColname,
                  sampleColname=sampleColname,
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
      xStart <- GenomicRanges::start(grl@unlistData);
      xEnd <- GenomicRanges::end(grl@unlistData);
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
      if ("yStart" %in% colnames(GenomicRanges::values(grl@unlistData))) {
         yStart <- GenomicRanges::values(grl@unlistData)$yStart;
      } else {
         yStart <- rep(0, length(grl@unlistData));
      }
      if ("yEnd" %in% colnames(GenomicRanges::values(grl@unlistData))) {
         yEnd <- GenomicRanges::values(grl@unlistData)$yEnd;
      } else {
         yEnd <- rep(0, length(grl@unlistData));
      }

      ## Calculate the height of each junction arc
      yScore <- GenomicRanges::values(grl@unlistData)[[scoreColname]] * 1;
      yMaxStartEnd <- pmax(abs(yStart), abs(yEnd));
      if (verbose) {
         jamba::printDebug("grl2df(): ",
            "calling internal_junc_score()");
         if (verbose > 1) {
            print(head(as.data.frame(grl@unlistData), 20));
            print(tail(as.data.frame(grl@unlistData), 20));
         }
      }
      ## sampleColname should be empty if there is no sample_id
      sampleColname <- intersect("sample_id",
         colnames(GenomicRanges::values(grl@unlistData)));

      internalMaxScore <- internal_junc_score(grl@unlistData,
         scoreColname=scoreColname,
         sampleColname=sampleColname,
         verbose=verbose,
         ...);
      if (verbose) {
         jamba::printDebug("grl2df(): ",
            "internalMaxScore:", internalMaxScore);
      }

      #yHeight <- pmax(yMaxStartEnd, internalMaxScore) - internalMaxScore
      arcBaseHeight <- pmax(abs(yScore), scoreArcMinimum) * scoreArcFactor;
      internalArcBaseHeight <- abs(internalMaxScore) * (1.1+scoreArcFactor);
      #yBaseHeight <- noiseFloor(internalArcBaseHeight - yMaxStartEnd,
      #   minimum=scoreArcMinimum*scoreArcFactor);
      yBaseHeight <- jamba::noiseFloor(internalArcBaseHeight - yStart,
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
         jamba::makeNames(
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
                  S4Vectors::elementNROWS(grl)),
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
         for (iCol in colnames(GenomicRanges::values(grl@unlistData))) {
            df[,iCol] <- rep(GenomicRanges::values(grl@unlistData)[[iCol]],
               each=8);
         }
      }
      if (keepGRLvalues) {
         for (iCol in colnames(GenomicRanges::values(grl))) {
            df[,iCol] <- rep(
               rep(GenomicRanges::values(grl)[,iCol],
                  S4Vectors::elementNROWS(grl)),
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
#' @param gene `character` string of the gene to plot, compared
#'    with `names(flatExonsByGene)` and `values(flatExonsByTx)$gene_name`.
#' @param tx `character` vector of the transcripts to plot, useful
#'    when specifying specific transcripts. Values are matched with
#'    `names(flatExonsByTx)`.
#' @param flatExonsByGene,flatExonsByTx `GRangesList` objects, named
#'    by `"gene_name"` or `"transcript_id"` respectively, containing
#'    disjoint (non-overlapping) exons within each GRangesList
#'    element. The data is expected to be in the form provided by
#'    `flattenExonsBy()`.
#' @param geneColor `character` color used as the base color for
#'    exons, where the color is varied for each feature type or
#'    subclass.
#' @param labelExons `logical` indicating whether to print text
#'    labels beneath each exon, using the values in colname
#'    `"gene_nameExon"`. Typically the gene and transcripts are
#'    named using consistent names, in which case one exon label
#'    is placed at the bottom of the lowest transcript for each
#'    unique exon label.
#' @param exonLabelAngle `numeric` angle in degrees (0 to 360)
#'    indicating how to rotate exon labels, where `90` is
#'    vertical, and `0` is horizontal.
#' @param exonLabelSize `numeric` value or `unit` object from `grid::unit()`.
#'    Numeric values are assumed to have unit `"pt"` which refers to
#'    font point size. Used to size exon labels when `labelExons=TRUE`.
#' @param newValues argument passed to `addGRLgaps()` to fill
#'    column values for newly created gap entries. It is useful
#'    to have `feature_type="gap"` so gaps have a different value
#'    than exons. It is also useful to have `subclass="gap"`
#'    when there are `"cds"` and `"noncds"` entries in the
#'    provided `flatExonsByGene` data.
#' @param gene_order `character` value indicating whether the
#'    flattened gene model should be plotted `"first"` above the
#'    transcript exon models, or `"last"` and below the
#'    transcript exon models.
#' @param return_type `character` value indicating whether to return
#'    the ggplot graphic object `"grob"`, or the data.frame
#'    `"df"` used to create the ggplot object.
#' @param ref2c `list` output from `make_ref2compressed()` which
#'    contains among other things, the `trans_grc` data of
#'    class `trans` or `transform` depending upon the versions
#'    of `scales` and `ggplot2` packages. It is used by
#'    `ggplot2::coord_trans()` or `ggplot2::scale_x_continuous()`.
#'    Note: The use of `trans` or `transform` object types should be
#'    consistent with the version of `scales` and `ggplot2`,
#'    for example an older version from cached data cannot be
#'    used with newer version of `ggplot2`. In that case the
#'    remedy is to delete the cache and start anew. Specifically,
#'    delete `sashimi_memoise`.
#' @param hjust,vjust `numeric` value to position exon labels
#'    passed to `ggrepel::geom_text_repel()`.
#' @param direction `character` string passed to `ggrepel::geom_text_repel()`
#'    to restrict placement of labels to one axis direction.
#' @param compressGaps `logical` indicating whether to compress gaps
#'    between exons. When `ref2c` is supplied, this argument is
#'    ignored and the supplied `ref2c` is used directly.
#' @param tx2geneDF `data.frame` or NULL, optionally used to help
#'    identify matching transcripts for the requested `gene` value,
#'    used when `"gene_name"` is not present in `values(flatExonsByTx)`.
#' @param label_coords `numeric` vector length 2, optional range of
#'    genomic coordinates to restrict labels, so labels are not
#'    arranged by `ggrepel::geom_text_repel()` even when `coord_cartesian()`
#'    is used to zoom into a specific x-axis range.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to relevant functions
#'    as needed, including `make_ref2compressed()`.
#'
#' @examples
#' ## Assume we start with flattened gene exons
#' data(test_exon_wide_gr);
#' test_flatExonsByGene <- GenomicRanges::split(test_exon_wide_gr,
#'    GenomicRanges::values(test_exon_wide_gr)[,"gene_name"]);
#'
#' # The most basic plot of exons
#' gene2gg(gene="TestGene1", flatExonsByGene=test_flatExonsByGene);
#'
#' # You can be fancy and number the exons
#' test_flatExonsByGene <- assignGRLexonNames(test_flatExonsByGene,
#'    geneSymbolColname="gene_name");
#' gene2gg(gene="TestGene1", flatExonsByGene=test_flatExonsByGene);
#'
#' # Or the exon labels can be hidden
#' gene2gg(gene="TestGene1", flatExonsByGene=test_flatExonsByGene, labelExons=FALSE)
#'
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
 exonLabelSize=8,
 geneSymbolColname="gene_name",
 newValues=list(feature_type="gap",
    subclass="gap",
    gene_nameExon="gap"),
 gene_order=c("first",
    "last"),
 return_type=c("grob",
    "df"),
 ref2c=NULL,
 hjust=0.5,
 vjust=0.5,
 direction=c("both",
    "x",
    "y"),
 compressGaps=TRUE,
 tx2geneDF=NULL,
 label_coords=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is a lightweight wrapper around grl2df() specifically intended
   ## for gene and transcript exon structure
   direction <- match.arg(direction);
   if (!suppressPackageStartupMessages(jamba::check_pkg_installed("ggplot2"))) {
      stop("gene2gg() requires ggplot2.");
   }
   if (!suppressPackageStartupMessages(jamba::check_pkg_installed("ggforce"))) {
      stop("gene2gg() requires ggforce.");
   }
   gene_order <- match.arg(gene_order);
   return_type <- match.arg(return_type);

   ## Convert exonLabelSize to have proper units, default "pt"
   if (!grid::is.unit(exonLabelSize)) {
      exonLabelSize <- grid::unit(exonLabelSize, "pt");
   }
   exonLabelMm <- grid::convertUnit(exonLabelSize, "mm", valueOnly=TRUE);

   if (length(flatExonsByGene) > 0 && length(gene) > 0) {
      if (!any(gene %in% names(flatExonsByGene))) {
         stop("gene was not found in names(flatExonsByGene)");
      }
      if (verbose) {
         jamba::printDebug("gene2gg(): ",
            "flatExonsByGene[gene]");
      }
      grl1a <- flatExonsByGene[names(flatExonsByGene) %in% gene];
   } else {
      grl1a <- NULL;
   }
   if (length(flatExonsByTx) > 0 && jamba::igrepHas("GRanges", class(flatExonsByTx))) {
      grl1 <- NULL;
      if (length(gene) > 0) {
         if (geneSymbolColname %in% colnames(GenomicRanges::values(flatExonsByTx))) {
            if (verbose) {
               jamba::printDebug("gene2gg(): ",
                  "values(flatExonsByTx)$gene_name %in% gene");
            }
            grl1 <- subset(flatExonsByTx, gene_name %in% gene);
         } else if (length(tx2geneDF) > 0 &&
            geneSymbolColname %in% colnames(tx2geneDF)) {
            if (verbose) {
               jamba::printDebug("gene2gg(): ",
                  "subset(tx2geneDF, ", geneSymbolColname, " %in% gene)$transcript_id");
            }
            tx <- unique(c(tx,
               subset(tx2geneDF, tx2geneDF[[geneSymbolColname]] %in% gene)$transcript_id));
         }
      }
      if (length(tx) > 0) {
         grl1 <- GenomicRanges::GRangesList(c(grl1,
            flatExonsByTx[names(flatExonsByTx) %in% tx]))
         GenomicRanges::values(grl1)[,geneSymbolColname] <- tx2geneDF[match(names(grl1),
            tx2geneDF$transcript_id),geneSymbolColname];
         GenomicRanges::values(grl1@unlistData)[,geneSymbolColname] <- rep(
            GenomicRanges::values(grl1)[,geneSymbolColname],
            S4Vectors::elementNROWS(grl1));
         GenomicRanges::values(grl1)$transcript_id <- names(grl1);
         GenomicRanges::values(grl1@unlistData)$transcript_id <- rep(
            GenomicRanges::values(grl1)$transcript_id,
            S4Vectors::elementNROWS(grl1));
      }
   } else {
      grl1 <- NULL;
   }
   if (length(grl1) == 0 && length(grl1a) == 0) {
      return(NULL);
   }
   if (length(grl1) == 0) {
      grl1a1 <- grl1a;
   } else if (length(grl1a) == 0) {
      grl1a1 <- rev(grl1)
   } else if ("first" %in% gene_order) {
      grl1a1 <- GenomicRanges::GRangesList(
         jamba::rmNULL(
            c(rev(grl1),
               grl1a)));
   } else {
      grl1a1 <- GenomicRanges::GRangesList(c(
         grl1a,
         rev(grl1)));
   }
   if (verbose) {
      jamba::printDebug("class(grl1a1):", class(grl1a1));
   }
   if (length(grl1a1) == 0) {
      stop("no exon models found for the gene and tx arguments given.");
   }
   ## Convert to tall data.frame
   grl1a1df <- grl2df(grl1a1,
      newValues=newValues,
      ...);

   ## Refactor to colorization
   ## - apply geneColor as a vector to each GRangesList named entry
   ## - apply gradient by subclass within each GRangesList named color
   if (length(geneColor) == 2) {
      geneColor <- jamba::nameVector(
         rep(geneColor,
            c(length(grl1a), length(grl1))),
         c(names(grl1a), names(grl1)));
   }
   gc <- jamba::nameVector(rep(geneColor,
      length.out=length(grl1a1)),
      names(grl1a1));
   if (length(names(geneColor)) > 0) {
      gc[names(geneColor)] <- geneColor;
   }
   if (verbose) {
      jamba::printDebug("grl2df(): ",
         names(gc),
         fgText=list("orange", setTextContrastColor(gc)),
         bgText=list(NA, gc));
   }
   colorColname <- head(
      intersect(c("subclass", "feature_type"),
         colnames(grl1a1df)),
      1);
   ## For now, color by grl_name and subclass
   ## In future, define colors for each gr_name, which will allow
   ## highlighting individual features if needed.
   if (length(colorColname) > 0) {
      grl1a1df$color_by <- jamba::pasteByRow(grl1a1df[,c("grl_name", colorColname)]);
      subclassV <- jamba::provigrep(
         c("noncds", "utr", "cds", "exon", "intron", "gap", "^NA$", "."),
         jamba::rmNA(naValue="NA",
            unique(grl1a1df[[colorColname]]))
      )
      colorSubV <- unlist(lapply(names(grl1a1), function(iname){
         jamba::color2gradient(
            jamba::nameVector(
               rep(gc[iname],
                  length.out=length(subclassV)),
               paste0(iname, "_", subclassV)));
      }));
   } else {
      ## If for some reason there is no subclass,feature_type
      ## we just color by the gene itself
      grl1a1df$color_by <- grl1a1df$grl_name;
      colorSubV <- gc;
   }

   ## Make a data.frame to label each exon
   exonLabelDF <- NULL;
   if (labelExons) {
      exonLabelDF <- jamba::renameColumn(
         from="groupBy",
         to="id",
         shrinkMatrix(grl1a1df[,c("x","y")],
            groupBy=grl1a1df[,"id"]));
      exonColname <- intersect(paste0(geneSymbolColname, "Exon"),
         colnames(grl1a1df));
      if (length(exonColname) == 0) {
         jamba::printDebug("gene2gg(): ",
            "Warning: exonColname not found in grl1a1df, skipping exon labels.",
            fgText=c("darkorange", "red"))
         labelExons <- FALSE;
      } else {
         exonLabelDF[,exonColname] <- grl1a1df[match(exonLabelDF[,"id"],
            grl1a1df[,"id"]), exonColname];

         exonLabelDF$y <- shrinkMatrix(grl1a1df[,c("x","y")],
            groupBy=grl1a1df[,"id"],
            shrinkFunc=min)$y;
         # Remove gap labels
         exonLabelDF <- subset(exonLabelDF,
            !grepl("^gap$|,", gene_nameExon));
         # Optionally remove labels outside the label_coords range
         if (length(label_coords) > 0 && nrow(exonLabelDF) > 0) {
            exonLabelDF <- subset(exonLabelDF,
               x >= min(label_coords) &
                  x <= max(label_coords));
         }
      }
   }
   if (length(ref2c) == 0) {
      if (compressGaps) {
         if (verbose) {
            jamba::printDebug("gene2gg(): ",
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
   ## based upon the number of transcripts being displayed.
   ## Todo: Clean up this code - non-overlapping labels are still
   ## a jumbled mess very often.
   ymin <- (-0.5 +
         -1 * (exonLabelMm/3) *
         ((labelExons*1) * length(grl1a1))/2);
   ymin <- -0.5;

   ## Put it together
   grl1a1gg <- ggplot2::ggplot(grl1a1df,
         ggplot2::aes(x=x,
            y=y,
            fill=color_by,
            color=color_by,
            text=paste0(sub("_", "<br>", gr_name),
               "<br>",subclass,
               "<br>",
               seqnames,
               ":",
               paste(scales::comma(range(x), accuracy=1),
                  collapse="-")
            ),
            group=id)) +
      ggforce::geom_shape(show.legend=FALSE) +
      colorjam::theme_jam() +
      ggplot2::ylab("") +
      ggplot2::scale_color_manual(
         na.value=jamba::makeColorDarker(tail(colorSubV, 1), darkFactor=1.3),
         values=jamba::makeColorDarker(colorSubV, darkFactor=1.3)) +
      ggplot2::scale_fill_manual(
         na.value=jamba::alpha2col(tail(colorSubV, 1), alpha=1),
         values=jamba::alpha2col(colorSubV, alpha=1)) +
      ggplot2::scale_y_continuous(breaks=seq_along(grl1a1)-1,
         limits=c(ymin, length(grl1a1)-0.5),
         expand=ggplot2::expansion(mult=c(labelExons * 2, 0)),
         labels=names(grl1a1));
   if (length(ref2c) > 0) {
      grl1a1gg <- grl1a1gg +
         ggplot2::scale_x_continuous(trans=ref2c$trans_grc);
   }
   if (labelExons && length(exonLabelDF) > 0 && nrow(exonLabelDF) > 0) {
      if (1 || length(ref2c) == 0) {
         grl1a1gg <- grl1a1gg +
            ggrepel::geom_text_repel(
               inherit.aes=FALSE,
               data=exonLabelDF,
               ggplot2::aes_(x=~x,
                  y=~min(y),
                  #text=NULL,
                  label=as.name(exonColname)),
               angle=exonLabelAngle,
               hjust=vjust,
               vjust=hjust,
               nudge_y=-0.2,
               ylim=c(NA, ymin),
               segment.color="grey35",
               min.segment.length=0,
               #fill="white",
               size=exonLabelMm,
               direction=direction);
      } else {
         grl1a1gg <- grl1a1gg +
            ggplot2::geom_text(
               inherit.aes=FALSE,
               data=exonLabelDF,
               ggplot2::aes(x=x,
                  y=min(y),
                  #text=NULL,
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
         ggplot2::ggtitle(paste0(jamba::cPaste(gene), " exon model"));
   }
   if ("df" %in% return_type) {
      return(grl1a1df);
   }
   if (length(ref2c) > 0) {
      attr(grl1a1gg, "ref2c") <- ref2c;
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
#' When the input data does not contain columns `"nameFrom"` and
#' and `"nameTo"`, the junctions are by default stacked by
#' coordinates.
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
#' @param scoreColname character string matching one of `colnames(GenomicRanges::values(gr))`
#'    that contains a numeric value representing the abundance of
#'    each splice junction observed.
#' @param sampleColname character string with the column or columns
#'    that contain biological sample identifier, used to ensure junctions
#'    are only stacked within a sample, and not across samples. When
#'    `sampleColname` is `NULL`, all junctions are stacked.
#' @param scoreFactor numeric value multiplied by the value in `scoreColname`
#'    to allow scaling the junctions across samples. Note that
#'    `scoreFactor` can be a vector, which would be applied to the
#'    vector of scores.
#' @param matchFrom,matchTo optional colnames to use when grouping
#'    junctions at the start and end positions. By default `"nameFrom"`
#'    and `"nameTo"` are used, as output from `closestExonToJunctions()`,
#'    which has the benefit of grouping junctions within the
#'    `spliceBuffer` distance from exon boundaries. If those values
#'    are not present `colnames(GenomicRanges::values(gr))` then the new
#'    default `c("seqnames", "start", "strand")` is used for
#'    `matchFrom`, and `c("seqnames", "end", "strand")` is used
#'    for `matchTo`. That said, if `matchFrom` or `matchTo` are supplied,
#'    those colnames are used from `as.data.frame(gr)`. Multiple colnames
#'    are allowed.
#'    Note also that `sampleColname` is appended to `matchFrom` and `matchTo`
#'    to ensure that matching is only performed within each
#'    `sampleColname` value.
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
#'    ranges=IRanges::IRanges(
#'       start=c(100, 300, 500,  900),
#'         end=c(200, 400, 750, 1000)),
#'    strand=rep("+", 4));
#' names(grExons) <- jamba::makeNames(rep("exon", length(grExons)),
#'    suffix="");
#'
#' grJunc <- GRanges(seqnames=rep("chr1", 6),
#'    ranges=IRanges::IRanges(start=c(200, 200, 400, 400, 750, 750),
#'       end=c(300, 500, 500, 900, 900, 1200)),
#'    strand=rep("+", 6),
#'    score=c(200, 50, 160, 40, 210, 10));
#' names(grJunc) <- jamba::makeNames(rep("junc", length(grJunc)));
#'
#' # quick plot showing exons and junctions using rectangles
#' grl <- c(
#'    GenomicRanges::GRangesList(exons=grExons),
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
#' grJunc <- closestExonToJunctions(grJunc, grExons, spliceBuffer=5)$spliceGRgene;
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
#'    ggtitle("Junctions stacked at boundaries");
#'
#' ## Another view showing the junction_rank
#' ## based upon max reads entering and exiting each exon edge
#' ggplot(grlJunc2df2, aes(x=x, y=y, group=gr_name)) +
#'    geom_diagonal_wide_arc(aes(alpha=junction_rank), fill="orange") +
#'    scale_alpha_manual(values=c(`1`=0.4, `2`=0.6, `3`=0.7)) +
#'    colorjam::scale_fill_jam() +
#'    colorjam::theme_jam() +
#'    ggtitle("Junctions stacked at boundaries")
#'
#' ## Last example showing how two samples are kept separate
#' grJunc_samples <- c(grJunc, grJunc);
#' values(grJunc_samples)[,"sample_id"] <- rep(c("SampleA","SampleB"),
#'    each=length(grJunc));
#' names(grJunc_samples) <- jamba::makeNames(GenomicRanges::values(grJunc_samples)[,"sample_id"]);
#' grlJunc2df_samples <- grl2df(grJunc_samples,
#'    scoreArcMinimum=20,
#'    shape="junction");
#' ggplot(grlJunc2df_samples, aes(x=x, y=y, group=gr_name, fill=gr_name)) +
#'    geom_diagonal_wide_arc(alpha=0.7,
#'       show.legend=FALSE) +
#'    colorjam::scale_fill_jam() +
#'    colorjam::theme_jam() +
#'    ggtitle("Junctions stacked at boundaries") +
#'    facet_wrap(~sample_id)
#'
#' @export
stackJunctions <- function
(gr,
 scoreColname="score",
 sampleColname="sample_id",
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
   if (!scoreColname %in% colnames(GenomicRanges::values(gr))) {
      stop("The scoreColname must be present in colnames(GenomicRanges::values(gr))");
   }

   ## Validate matchFrom, matchTo
   if (length(matchFrom) == 0) {
      if ("nameFrom" %in% colnames(GenomicRanges::values(gr))) {
         matchFrom <- "nameFrom";
      } else {
         matchFrom <- c("seqnames", "start");
      }
   }
   matchFrom <- unique(c(matchFrom, sampleColname));
   if (length(matchTo) == 0) {
      if ("nameTo" %in% colnames(GenomicRanges::values(gr))) {
         matchTo <- "nameTo";
      } else {
         matchTo <- c("seqnames", "end");
      }
   }
   matchTo <- unique(c(matchTo, sampleColname));
   ##
   grValueColnames <- colnames(GenomicRanges::values(gr));
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

   ## Optionally enforce strandedness
   if (strandedScore) {
      if (verbose) {
         jamba::printDebug("stackJunctions(): ",
            "Adjusting score by strand, using scoreFactor:",
            scoreFactor);
      }
      matchFrom <- unique(c(matchFrom, "strand"));
      matchTo <- unique(c(matchTo, "strand"));
      scoreV <- scoreFactor *
         ifelse(as.character(GenomicRanges::strand(gr)) %in% "-", -1, 1) *
         abs(GenomicRanges::values(gr)[[scoreColname]]);
   } else {
      if (verbose) {
         jamba::printDebug("stackJunctions(): ",
            "Adjusting unstranded score, using scoreFactor:",
            scoreFactor);
      }
      scoreV <- scoreFactor *
         GenomicRanges::values(gr)[[scoreColname]];
   }
   GenomicRanges::values(gr)[[scoreColname]] <- scoreV;

   if (verbose > 1) {
      jamba::printDebug("stackJunctions(): ",
         "matchFrom:", matchFrom,
         ", matchTo:", matchTo);
   }

   # Combine value colnames and non-value colnames, using only
   # the required colnames since as.data.frame(gr) will fail if
   # any columns are not coercible to data.frame, for example list.
   exonsFrom <- jamba::pasteByRow(sep="_",
      as.data.frame(gr[,intersect(matchFrom, grValueColnames)])[,matchFrom,drop=FALSE])
   exonsTo <- jamba::pasteByRow(sep="_",
      as.data.frame(gr[,intersect(matchTo, grValueColnames)])[,matchTo,drop=FALSE])

   exonsFrom <- gsub("[.][-]*[0-9]+$", "",
      exonsFrom);
   exonsTo <- gsub("[.][-]*[0-9]+$", "",
      exonsTo);
   if (verbose > 1) {
      jamba::printDebug("stackJunctions(): ",
         "head(exonsFrom):",
         head(exonsFrom, 10));
      jamba::printDebug("stackJunctions(): ",
         "head(exonsTo):",
         head(exonsTo, 10));
   }

   ## Extend baseline to length of gr, so the baseline applies
   ## to each exon
   allExons <- jamba::mixedSort(unique(
      c(exonsFrom, exonsTo)));
   baselineV <- jamba::nameVector(
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
   if (verbose) {
      jamba::printDebug("stackJunctions(): ",
         "Stacking exonsFrom");
   }
   if (length(jamba::tcount(names(gr), minCount=2)) > 0) {
      names(gr) <- jamba::makeNames(names(gr));
   }
   order1 <- do.call(order, list(exonsFrom, GenomicRanges::width(gr)));
   #values(gr)[order1,"yStart"] <- shrinkMatrix(
   yStart_df <- shrinkMatrix(
      scoreV[order1],
      groupBy=exonsFrom[order1],
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))});
   yRow_df <- shrinkMatrix(
      names(gr)[order1],
      groupBy=exonsFrom[order1],
      shrinkFunc=c);
   if (verbose > 1) {
      jamba::printDebug("head(yStart_df):");
      print(head(yStart_df, 20));
      print(gr[yRow_df$x]);
      jamba::printDebug("exonsFrom:", exonsFrom);
      jamba::printDebug("baselineV:");
      print(baselineV);
      jamba::printDebug("baselineV[exonsFrom]:", baselineV[exonsFrom]);
   }
   order1rev <- match(names(gr), yRow_df$x);
   GenomicRanges::values(gr)[,"yStart"] <- yStart_df$x[order1rev] + baselineV[exonsFrom];

   ## End position
   ## group by exonsTo which allows the spliceBuffer to work
   if (verbose) {
      jamba::printDebug("stackJunctions(): ",
         "Stacking exonsTo");
   }
   order2 <- do.call(order, list(exonsTo, GenomicRanges::width(gr)));
   yEnd_df <- shrinkMatrix(
      scoreV[order2],
      groupBy=exonsTo[order2],
      shrinkFunc=function(x){cumsum(head(c(0,x), length(x)))});
   yRow_df <- shrinkMatrix(
      names(gr)[order2],
      groupBy=exonsTo[order2],
      shrinkFunc=c);
   order2rev <- match(names(gr), yRow_df$x);
   GenomicRanges::values(gr)[,"yEnd"] <- yEnd_df$x[order2rev] + baselineV[exonsTo];

   ## Experimental "fix" for negative strand using flipped stacking
   if (any(as.character(GenomicRanges::strand(gr)) %in% "-")) {
      is_neg <- (as.character(GenomicRanges::strand(gr)) %in% "-");
      GenomicRanges::values(gr)[is_neg, c("yEnd", "yStart")] <- (
         GenomicRanges::values(gr)[is_neg, c("yStart", "yEnd")]);
   }

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
   ## Ensure "nameFromTo" exists
   if (all(c("nameFrom", "nameTo") %in% colnames(GenomicRanges::values(gr)))) {
      if (!"nameFromTo" %in% colnames(GenomicRanges::values(gr))) {
         GenomicRanges::values(gr)[,"nameFromTo"] <- jamba::pasteByRow(GenomicRanges::values(gr)[,c("nameFrom", "nameTo")],
            sep=" ");
      }
      sampleColname <- intersect(sampleColname, colnames(GenomicRanges::values(gr)));
      if (length(sampleColname) == 0) {
         ## Stack without using sample_id
         value_colnames <- c(scoreColname, "nameFromTo");
         from_colnames <- c("nameFrom");
         to_colnames <- c("nameTo");
         nfts_colname <- "nameFromTo";
      } else {
         ## Stack within sample_id
         if (!"nameFromToSample" %in% colnames(GenomicRanges::values(gr))) {
            GenomicRanges::values(gr)[,"nameFromToSample"] <- jamba::pasteByRow(GenomicRanges::values(gr)[,c("nameFromTo", sampleColname)],
               sep=":!:");
         }
         value_colnames <- c(scoreColname, "nameFromToSample");
         from_colnames <- c(sampleColname, "nameFrom");
         to_colnames <- c(sampleColname, "nameTo");
         nfts_colname <- "nameFromToSample";
      }
      if (verbose) {
         jamba::printDebug("stackJunctions(): ",
            "Calculating junction ranks with sampleColname:",
            sampleColname);
      }
      juncRankFrom <- shrinkMatrix(
         GenomicRanges::values(gr)[,value_colnames],
         jamba::pasteByRow(GenomicRanges::values(gr)[,from_colnames]),
         shrinkFunc=shrink_junc);
      juncRankTo <- shrinkMatrix(
         GenomicRanges::values(gr)[,value_colnames],
         jamba::pasteByRow(GenomicRanges::values(gr)[,to_colnames]),
         shrinkFunc=shrink_junc);
      juncRankDF <- data.frame(
            juncRankFrom[match(GenomicRanges::values(gr)[,nfts_colname], juncRankFrom[,nfts_colname]),],
            juncRankTo[match(GenomicRanges::values(gr)[,nfts_colname], juncRankTo[,nfts_colname]),]
      );
      juncRank <- (1 +
            juncRankFrom[match(GenomicRanges::values(gr)[,nfts_colname], juncRankFrom[,nfts_colname]),scoreColname] +
            juncRankTo[match(GenomicRanges::values(gr)[,nfts_colname], juncRankTo[,nfts_colname]),scoreColname]
      );
      GenomicRanges::values(gr)[,"junction_rank"] <- factor(juncRank);
   }
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
#' @param junc_alpha numeric value between 0 and 1, to define the
#'    alpha transparency used for junction colors, where 0 is
#'    fully transparent, and 1 is completely non-transparent.
#' @param junc_nudge_pct `numeric` value to nudge junction labels
#'    by the percent of the maximum y-axis junction label position.
#'    * The default `junc_nudge_pct=0.05` nudges labels up by 5%,
#'    which makes them consistently appear above the top edge of
#'    the junction ribbon.
#'    * Negative values will place labels just below the top edge of
#'    the junction ribbon.
#'    * A vector can be supplied to nudge each junction label
#'    individually, applied in order the labels appear from
#'    left to right.
#'    * Note that when the distance from label to the top edge
#'    of the ribbon exceeds a threshold, a line segment is drawn
#'    from the label to the top edge of the junction ribbon.
#'    This threshold is controlled by
#'    `ggrepel::geom_text_repel(..., min.segment.length=0.5)` and
#'    is not configurable at this time.
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
#' @param ylabel character string used as the y-axis label, by default
#'    `"score"` reflects the coverage score and junction score,
#'    respectively for coverage and junction data. Scores are
#'    also adjusted using the `scale_factor` value for each
#'    `sample_id` as defined in the `filesDF`. Set to `NULL` to
#'    hide the y-axis label completely.
#' @param xlabel character string used to define the x-axis name,
#'    which takes priority over argument `xlabel_ref`. When
#'    `xlabel` is `NULL` and `xlabel_ref` is `FALSE`, then the
#'    x-axis name is `""`, which displays no x-axis label.
#' @param xlabel_ref logical indicating whether the x-axis name
#'    should be determined by the reference (chromosome).
#' @param use_jam_themes logical indicating whether to apply
#'    `colorjam::theme_jam()`, by default for the ggplot theme.
#' @param apply_facet logical indicating whether to apply
#'    `ggplot2::facet_wrap()` with `"~sample_id"` defining each
#'    panel.
#' @param facet_scales character value used as `"scales"` argument in
#'    `ggplot2::facet_wrap()` when `apply_facet=TRUE`.
#' @param ref2c optional output from `make_ref2compressed()` used to
#'    compress axis coordinates during junction arc calculations.
#' @param label_coords numeric vector length 2, optional range of
#'    genomic coordinates to restrict labels, so labels are not
#'    arranged by `ggrepel::geom_text_repel()` even when `coord_cartesian()`
#'    is used to zoom into a specific x-axis range.
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
 junc_color=jamba::alpha2col("goldenrod2", 0.3),
 junc_fill=jamba::alpha2col("goldenrod2", 0.9),
 junc_alpha=0.8,
 junc_accuracy=1,
 junc_nudge_pct=0.05,
 fill_scheme=c("sample_id", "exon"),
 color_sub=NULL,
 ylabel="read depth",
 xlabel=NULL,
 xlabel_ref=TRUE,
 use_jam_themes=TRUE,
 apply_facet=TRUE,
 facet_scales="free_y",
 ref2c=NULL,
 label_coords=NULL,
 do_highlight=FALSE,
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
      if ("ref2c" %in% names(attributes(sashimi))) {
         sashimi$ref2c <- attr(sashimi$df, "ref2c");
      }
   }
   if (!"ref2c" %in% names(sashimi)) {
      coord_method <- "none";
   }
   ggCov <- NULL;
   cov_rows <- (sashimi$df$type %in% "coverage");
   if ("coverage" %in% show && any(cov_rows)) {
      if ("exon" %in% fill_scheme) {
         sashimi$df$color_by[cov_rows] <- as.character(sashimi$df$gr[cov_rows]);
      } else {
         sashimi$df$color_by[cov_rows] <- as.character(sashimi$df$sample_id[cov_rows]);
      }
      # Define colors
      color_sub_cov <- colorjam::group2colors(unique(sashimi$df$color_by[cov_rows]));
      use_names <- setdiff(names(color_sub_cov), names(color_sub));
      color_sub[use_names] <- color_sub_cov[use_names];

      exonlabel_rows <- (sashimi$df$type %in% "exon_label");
   }

   ## Junction data
   junc_rows <- (sashimi$df$type %in% "junction");
   #if ("junction" %in% show && "juncDF" %in% names(sashimi)) {
   if ("junction" %in% show && any(junc_rows)) {
      color_sub_junc <- NULL;
      if (!"junction_rank" %in% colnames(sashimi$df)) {
         sashimi$df$junction_rank <- 3;
      }
      if (fill_scheme %in% "sample_id") {
         sashimi$df$color_by[junc_rows] <- jamba::pasteByRow(
            sashimi$df[junc_rows,c("sample_id","junction_rank"), drop=FALSE],
            sep=".");
      } else {
         sashimi$df$color_by[junc_rows] <- as.character(sashimi$df$junction_rank[junc_rows]);
      }

      # order so the lower rank, lesser junctions, are drawn on top
      # disabled in version 46
      #juncDF <- juncDF[order(-juncDF$junction_rank),,drop=FALSE];

      ## gradually transparent white to shade by junction_rank
      junc_blank_3 <- jamba::nameVector(
         jamba::alpha2col(rep("#DDDDDD", 3), alpha=c(0.3, 0.67, 0.8)),
         c(3, 2, 1));
      if (fill_scheme %in% "sample_id") {
         if (all(unique(sashimi$df$sample_id[junc_rows]) %in% names(color_sub))) {
            color_sub_samples <- color_sub[intersect(sashimi$df$sample_id[junc_rows], names(color_sub))];
         } else {
            color_sub_samples <- colorjam::group2colors(unique(sashimi$df$sample_id[junc_rows]));
         }
         # blend with gradually transparent white
         color_sub_junc_3 <- colorspace::hex(
            colorspace::mixcolor(
               colorspace::hex2RGB(jamba::rgb2col(col2rgb(
                  rep(color_sub_samples, each=3)))),
               colorspace::hex2RGB(junc_blank_3),
               alpha=jamba::col2alpha(junc_blank_3)
            )
         );
         names(color_sub_junc_3) <- paste(names(color_sub_junc_3),
            names(junc_blank_3),
            sep=".");
      } else {
         # Otherwise blend junc_fill with white
         color_sub_junc_3 <- colorspace::hex(
            colorspace::mixcolor(
               colorspace::hex2RGB(jamba::rgb2col(col2rgb(rep(junc_fill[1], each=3)))),
               colorspace::hex2RGB(junc_blank_3),
               alpha=jamba::col2alpha(junc_blank_3)
            )
         );
         names(color_sub_junc_3) <- names(junc_blank_3);
      }
      use_names <- setdiff(names(color_sub_junc_3), names(color_sub));
      color_sub[use_names] <- color_sub_junc_3[use_names];
      color_sub[names(color_sub_junc_3)] <- jamba::alpha2col(
         color_sub[names(color_sub_junc_3)],
         alpha=junc_alpha);

      ## junction labels
      junclabel_rows <- (sashimi$df$type %in% "junction_label");
      if (any(junclabel_rows)) {
         if (fill_scheme %in% "sample_id") {
            sashimi$df$color_by[junclabel_rows] <- jamba::pasteByRow(sashimi$df[junclabel_rows, c("sample_id","junction_rank"), drop=FALSE],
               sep=".");
         } else {
            sashimi$df$color_by[junclabel_rows] <- as.character(sashimi$df$junction_rank[junclabel_rows]);
         }
      }
   }

   ## Pull out the data.frame
   cjDF <- sashimi$df;
   if (verbose) {
      jamba::printDebug("plotSashimi(): ",
         "head(df):");
      print(head(sashimi$df));
      jamba::printDebug("plotSashimi(): ",
         "table(sashimi$df$color_by):");
      print(table(sashimi$df$color_by));
   }

   # Create ggplot2 piece by piece
   if (do_highlight) {
      cjDFh <- plotly::highlight_key(cjDF,
         key=~feature);
      gg_sashimi <- ggplot2::ggplot(
         data=cjDFh,
         ggplot2::aes(
            x=x,
            y=y,
            group=name,
            color=color_by,
            fill=color_by,
            text=name
         ));
   } else {
      gg_sashimi <- ggplot2::ggplot(
         data=cjDF,
         ggplot2::aes(
            x=x,
            y=y,
            group=name,
            color=color_by,
            fill=color_by
         ));
   }
   # comma-delimit y-axis values
   # and remove the y-axis label
   gg_sashimi <- gg_sashimi +
      ggplot2::scale_y_continuous(
         labels=scales::comma,
         name=ylabel) +
      ggplot2::ylab(ylabel);

   # Add coverage layer
   color_sub_d <- jamba::makeColorDarker(color_sub, darkFactor=1.2);
   if (all(c("coverage","junction") %in% cjDF$type)) {
      if (verbose) {
         jamba::printDebug("plotSashimi(): ",
            "Including coverage and splice junctions.");
      }
      gg_sashimi <- gg_sashimi +
         #geom_polygon(
         ggforce::geom_shape(
            data=. %>% dplyr::filter(type %in% "coverage"),
            show.legend=FALSE) +
         geom_diagonal_wide_arc(
            data=. %>% dplyr::filter(type %in% "junction"),
            show.legend=FALSE,
            alpha=junc_alpha,
            strength=0.4);
   } else if ("coverage" %in% cjDF$type) {
      if (verbose) {
         jamba::printDebug("plotSashimi(): ",
            "Including coverage without splice junctions.");
      }
      gg_sashimi <- gg_sashimi +
         #geom_polygon(
         ggforce::geom_shape(
            data=. %>% dplyr::filter(type %in% "coverage"),
            show.legend=FALSE
         );
   } else if ("junction" %in% cjDF$type) {
      if (verbose) {
         jamba::printDebug("plotSashimi(): ",
            "Including splice junctions without coverage.");
      }
      gg_sashimi <- gg_sashimi +
         geom_diagonal_wide_arc(
            data=. %>% dplyr::filter(type %in% "junction"),
            show.legend=FALSE,
            alpha=junc_alpha,
            strength=0.4
         );
   }

   # Add junction_labels only for non-plotly
   if (jamba::igrepHas("junctionlabel|junction.label", show) &&
         "junction_label" %in% cjDF$type) {
      if (do_highlight) {
         if (verbose) {
            jamba::printDebug("plotSashimi(): ",
               "Disabled junction labels because do_highlight=",
               do_highlight);
         }
      } else {
         if (verbose) {
            jamba::printDebug("plotSashimi(): ",
               "Adding junction labels.");
         }
         if (length(label_coords) == 0) {
            label_coords <- range(
               subset(cjDF,type %in% "junction_label")$x,
               na.rm=TRUE);
         }
         ## Todo: nudge_y to adjust labels consistently above the ribbon
         max_junc_y <- max(subset(cjDF, type %in% "junction_label")$y,
            na.rm=TRUE);
         gg_sashimi <- gg_sashimi +
            ggrepel::geom_text_repel(
               data=. %>% dplyr::filter(type %in% "junction_label" &
                     x >= min(label_coords) & x <= max(label_coords)),
               angle=90,
               vjust=0.5,
               direction="y",
               point.padding=0,
               color="black",
               nudge_y=(max_junc_y * junc_nudge_pct),
               #fill="transparent",
               ggplot2::aes(
                  label=scales::comma(score,
                     accuracy=junc_accuracy)
               )
            );
      }
   }

   ## Determine an appropriate x-axis label
   if (length(xlabel) == 0) {
      if (length(xlabel_ref) > 0 && xlabel_ref) {
         xlabel <- as.character(unique(seqnames(attr(sashimi$ref2c, "gr"))));
      } else {
         xlabel <- "";
      }
   }

   if ("scale" %in% coord_method) {
      gg_sashimi <- gg_sashimi +
         ggplot2::scale_x_continuous(trans=sashimi$ref2c$trans_grc,
            name=xlabel);
   } else if ("coord" %in% coord_method) {
      gg_sashimi <- gg_sashimi +
         ggplot2::coord_trans(x=sashimi$ref2c$trans_grc) +
         #coord_cartesian(expand=FALSE) +
         ggplot2::xlab(xlabel);
   }

   if (use_jam_themes) {
      gg_sashimi <- gg_sashimi +
         colorjam::theme_jam(
            panel.grid.major.colour="grey80",
            panel.grid.minor.colour="grey90");
   }
   # Apply colors
   gg_sashimi <- gg_sashimi +
      ggplot2::scale_fill_manual(values=color_sub) +
      ggplot2::scale_color_manual(values=jamba::makeColorDarker(color_sub, darkFactor=1.2));

   # Apply facet
   if (apply_facet) {
      gg_sashimi <- gg_sashimi +
         ggplot2::facet_grid(sample_id~.,
            scales=facet_scales);
   }
   return(gg_sashimi);
}

# #' @importFrom plotly to_basic
# #' @export
# plotly::to_basic

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
#' This function is exported, since otherwise it caused problems
#' when invoking `plotly::ggplotly()` and this function was not
#' consistently used in that process for some reason.
#'
#' Currently only `data` is passed to `ggplot2::GeomPolygon`.
#'
#' @inheritParams plotly::to_basic
#'
#' @family jam ggplot2 functions
#'
#' @importFrom plotly to_basic
#' @export
to_basic.GeomShape <- function
(data,
 prestats_data,
 layout,
 params,
 p,
 ...)
{
   plotly:::prefix_class(data, "GeomPolygon");
}

