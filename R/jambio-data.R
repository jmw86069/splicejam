
#' Sample junction data GRangesList
#'
#' Sample junction data GRangesList
#'
#' This dataset contains RNA-seq splice junction data
#' stored as a GRangesList.
#'
#' @format GRangesList where each GRangesList item represents
#'    a distinct biological sample. Each GRanges item represents
#'    one splice junction, whose score represents the abundance
#'    of splice junction reads observed. The start of each splice
#'    junction should be one base after the end of the corresponding
#'    exon, using 1-based coordinates. Therefore, an exon spanning
#'    1-10 covers 10 bases, the corresponding junction would
#'    begin at position 11. Similarly, if the connected exon
#'    begins at position 100, the junction would end at position
#'    99.
#'
#' @family splicejam data
#'
#' @examples
#' # The code below is used to create the junction test data
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' test_junc_gr <- GRanges(seqnames=rep("chr1", 5),
#'    ranges=IRanges::IRanges(
#'       start=c(200, 200, 400, 400, 750),
#'         end=c(299, 499, 499, 899, 899)),
#'    strand=rep("+", 5),
#'    score=c(200, 50, 120, 80, 170),
#'    sample_id=rep("sample_A", 5));
#' names(test_junc_gr) <- jamba::makeNames(
#'    rep("junc", length(test_junc_gr)),
#'    suffix="");
#' test_junc_gr;
#'
#' # To plot junctions, use grl2df(..., shape="junction")
#' junc_df <- grl2df(test_junc_gr, shape="junction")
#' gg1 <- ggplot(junc_df, aes(x=x, y=y, group=id, fill=gr_name)) +
#'    ggforce::geom_diagonal_wide(alpha=0.7) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam()
#' print(gg1);
"test_junc_gr"

#' Sample junction data GRangesList with wide introns
#'
#' Sample junction data GRangesList with wide introns
#'
#' This dataset contains RNA-seq splice junction data
#' stored as a GRangesList.
#'
#' Intron and exon sizes are more consistent with
#' mammalian gene structure, and are intended to
#' demonstrate the challenge with visualizing exon
#' coverage data on a genomic scale. See examples
#' for steps to compress the intron sizes.
#'
#' @format GRangesList where each GRangesList item represents
#'    a distinct biological sample. Each GRanges item represents
#'    one splice junction, whose score represents the abundance
#'    of splice junction reads observed. The start of each splice
#'    junction should be one base after the end of the corresponding
#'    exon, using 1-based coordinates. Therefore, an exon spanning
#'    1-10 covers 10 bases, the corresponding junction would
#'    begin at position 11. Similarly, if the connected exon
#'    begins at position 100, the junction would end at position
#'    99.
#'
#' @family splicejam data
#'
#' @examples
#' # The code below is used to create the junction test data
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' data(test_junc_gr);
#' test_junc_wide_gr <- test_junc_gr;
#' xshift <- c(0, 10000, 20000, 39000);
#' end(test_junc_wide_gr) <- end(test_junc_gr) + xshift[c(2,3,3,4,4)];
#' start(test_junc_wide_gr) <- start(test_junc_gr) + xshift[c(1,1,2,2,3)];
#' names(test_junc_wide_gr) <- jamba::makeNames(
#'    rep("juncwide", length(test_junc_gr)),
#'    suffix="");
#' test_junc_wide_gr;
#'
#' # To plot junctions, use grl2df(..., shape="junction")
#' junc_wide_df <-grl2df(test_junc_wide_gr, shape="junction")
#' ggWide1 <- ggplot(junc_wide_df, aes(x=x, y=y, group=id, fill=gr_name)) +
#'    ggforce::geom_diagonal_wide(alpha=0.7) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    xlab("chr1") +
#'    ggtitle("junctions (full intron width)")
#' print(ggWide1);
#'
#' # The exons are required to define compressed ranges
#' # Note: DO NOT USE THESE STEPS
#' # this plot shows the distortion of arcs by compressed x-axis
#' # and we will fix it in the next example
#' data(test_exon_wide_gr);
#' ref2c <- make_ref2compressed(test_exon_wide_gr,
#'    nBreaks=10);
#' ggWide1c <- ggWide1 +
#'    scale_x_continuous(trans=ref2c$trans_grc) +
#'    xlab("chr1 (compressed introns)") +
#'    ggtitle("junctions (compressed introns, distorted)");
#' print(ggWide1c);
#'
#' # to fix the arc shapes, supply the transform to grl2df()
#' # Note: USE THESE STEPS
#' junc_wide_c_df <-grl2df(test_junc_wide_gr,
#'    shape="junction",
#'    ref2c=ref2c);
#' ggWide1c2 <- ggplot(junc_wide_c_df, aes(x=x, y=y, group=id, fill=gr_name)) +
#'    ggforce::geom_diagonal_wide(alpha=0.7) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    scale_x_continuous(trans=ref2c$trans_grc) +
#'    xlab("chr1 (compressed introns)") +
#'    ggtitle("junctions (compressed introns)");
#' print(ggWide1c2);
#'
"test_junc_wide_gr"

#' Sample exon data GRanges
#'
#' Sample exon data GRanges
#'
#' This dataset contains RNA-seq splice junction data
#' stored as a GRangesList.
#'
#' @format GRanges object where each segment represents one
#'    exon for an arbitrary gene. It has one column of values,
#'    `"gene_name"` used for Sashimi plot preparation.
#'
#' @family splicejam data
#'
#' @examples
#' # The code below is used to create the exon test data
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' test_exon_gr <- GRanges(seqnames=rep("chr1", 4),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 300, 500, 900),
#'         end=c(199, 399, 749, 999)),
#'    strand=rep("+", 4),
#'    gene_name=rep("TestGene1", 4));
#' names(test_exon_gr) <- jamba::makeNames(rep("exon", length(test_exon_gr)),
#'    suffix="");
#'
#' # To plot a simple GRanges object
#' ggplot(grl2df(test_exon_gr), aes(x=x, y=y, group=id, fill=feature_type)) +
#'    geom_polygon() +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam()
"test_exon_gr"

#' Sample exon data GRanges with wide introns
#'
#' Sample exon data GRanges with wide introns
#'
#' This dataset contains RNA-seq splice junction data
#' stored as a GRangesList.
#'
#' Intron and exon sizes are more consistent with
#' mammalian gene structure, and are intended to
#' demonstrate the challenge with visualizing exon
#' coverage data on a genomic scale. See examples
#' for steps to compress the intron sizes.
#'
#' @format GRanges object where each segment represents one
#'    exon for an arbitrary gene. It has one column of values,
#'    `"gene_name"` used for Sashimi plot preparation.
#'
#' @family splicejam data
#'
#' @examples
#' # The code below is used to create the exon test data
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' test_exon_wide_gr <- GRanges(seqnames=rep("chr1", 4),
#'    ranges=IRanges::IRanges(
#'       start=c(100, 10300, 20500,  39900),
#'         end=c(200, 10400, 20750, 40000)),
#'    strand=rep("+", 4),
#'    gene_name=rep("TestGene1", 4));
#' names(test_exon_wide_gr) <- jamba::makeNames(rep("wide", length(test_exon_wide_gr)),
#'    suffix="");
#' test_exon_wide_gr;
#'
#' # To plot a simple GRanges object
#' widedf <- grl2df(test_exon_wide_gr);
#' ggWide <- ggplot(widedf, aes(x=x, y=y, group=id, fill=feature_type)) +
#'    geom_polygon() +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    xlab("chr1") +
#'    ggtitle("exons (full introns)")
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
"test_exon_wide_gr"

#' Sample exon coverage data GRanges
#'
#' Sample exon coverage data GRanges
#'
#' This dataset represents exon GRanges with an additional column
#' with NumericList values representing RNA-seq read coverage
#' across these exons.
#'
#' @family splicejam data
#'
#' @examples
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#' data(test_exon_gr);
#'
#' # The steps below demonstrate how to create coverage data manually
#' test_cov_gr <- test_exon_gr;
#' set.seed(123);
#' cov_values <- c(250, 200, 180, 250);
#' values(test_cov_gr)[,"sample_A"] <- NumericList(
#'    lapply(seq_along(test_exon_gr), function(i){
#'       cov_values[i] +
#'          round(rnorm(width(test_exon_gr)[i])*7)
#'    }));
#' test_cov_gr;
#'
#' # You can plot coverage using exoncov2polygon()
#' exondf <- exoncov2polygon(test_cov_gr, covNames="sample_A");
#' gg3 <- ggplot(exondf,
#'       aes(x=x, y=y, group=gr, fill=gr, color=gr)) +
#'    ggforce::geom_shape(alpha=0.8) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    colorjam::scale_color_jam();
#' print(gg3);
#'
#' # you can add junctions to exons in one plot
#' junc_df <- grl2df(test_junc_gr, shape="junction")
#' gg3 +
#'    ggforce::geom_diagonal_wide(data=junc_df,
#'       inherit.aes=FALSE,
#'       show.legend=FALSE,
#'       fill="orange",
#'       aes(x=x, y=y, group=id),
#'       alpha=0.6)
"test_cov_gr"

#' Sample exon coverage data GRanges with wide introns
#'
#' Sample exon coverage data GRanges with wide introns
#'
#' This dataset represents exon GRanges with an additional column
#' with NumericList values representing RNA-seq read coverage
#' across these exons.
#'
#' Intron and exon sizes are more consistent with
#' mammalian gene structure, and are intended to
#' demonstrate the challenge with visualizing exon
#' coverage data on a genomic scale. See examples
#' for steps to compress the intron sizes.
#'
#' @family splicejam data
#'
#' @examples
#' # The steps below demonstrate how to create coverage data manually
#' suppressPackageStartupMessages(library(GenomicRanges));
#' suppressPackageStartupMessages(library(ggplot2));
#'
#' data(test_cov_gr);
#' data(test_exon_wide_gr);
#' test_cov_wide_gr <- test_cov_gr;
#' ranges(test_cov_wide_gr) <- ranges(test_exon_wide_gr);
#'
#' # You can plot coverage using exoncov2polygon()
#' widecovdf <- exoncov2polygon(test_cov_wide_gr, covNames="sample_A");
#' ggWide3 <- ggplot(widecovdf,
#'       aes(x=x, y=y, group=gr, fill=gr, color=gr)) +
#'    ggforce::geom_shape(alpha=0.7) +
#'    colorjam::theme_jam() +
#'    colorjam::scale_fill_jam() +
#'    colorjam::scale_color_jam();
#' print(ggWide3);
#'
#' # Now compress the introns keeping axis labels
#' ref2c <- make_ref2compressed(test_cov_wide_gr,
#'    nBreaks=10);
#' ggWide3c <- ggWide3 +
#'    scale_x_continuous(trans=ref2c$trans_grc) +
#'    xlab("chr1 (compressed introns)") +
#'    ggtitle("exons (compressed introns)");
#' print(ggWide3c);
#'
#' # you can add junctions to exons in one plot
#' junc_wide_c_df <-grl2df(test_junc_wide_gr,
#'    shape="junction",
#'    ref2c=ref2c);
#' ggWide3c +
#'    ggforce::geom_diagonal_wide(data=junc_wide_c_df,
#'       inherit.aes=FALSE,
#'       show.legend=FALSE,
#'       fill="orange",
#'       aes(x=x, y=y, group=id),
#'       alpha=0.6)
"test_cov_wide_gr"

