
# sashimi wrapper function
# - perform prepareSashimi(), plotSashimi()
# - apply logic similar to launchSashimiApp() for convenience


#' Make Splicejam Sashimi Figure
#' 
#' Make Splicejam Sashimi Figure, with optional columns, gene model,
#' subset region by coordinates or exon range.
#' 
#' This function is intended to mimic the functional steps
#' performed by `launchSashimiApp()` which are aimed to provide
#' a user-friendly interface to the many internal details.
#' 
#' Notable features:
#' 
#' * 'exon_range': new argument allows zooming into a genome
#' coordinate range using exon names.
#' * gene model is shown and aligned by default, using cowplot.
#' * 'layout_ncol': allows multiple columns with splicejam and
#' gene/transcript/exon model.
#' 
#' ## Todo:
#' 
#' * Consider providing interactive output via `ggplotly()`.
#' * Expand the approach to permit multiple genes.
#' * Apply memoise to plot data, as in `launchSashimiApp()`
#' 
#' @family splicejam core functions
#' 
#' @param sjenv `environment` containing data produced by
#'    `sashimiDataConstants()`, specifically:
#'    * 'flatExonsByGene'
#'    * 'flatExonsByTx'
#'    * 'detectedTx'
#'    * 'detectedGenes'
#'    * tx2geneDF': `data.frame` with columns: 'gene_name', 
#'    'transcript_id'.
#' @param gene `character` string with gene(s) to include in the
#'    flattened exon model.
#' @param gene_sd `list` optional sashimi data output from
#'    `prepareSashimi()`.
#' @param filesDF `data.frame` with columns 'url', 'sample_id', 'type'.
#'    See `prepareSashimi()` for more information.
#' @param sample_id `character` default NULL, optionally used to
#'    subset data to show only data for samples 'sample_id'.
#' @param exon_range `character` vector with first and last exons
#'    to include in the genome coordinate range. Typically in the
#'    form: 'genesymbol_exon1', 'genesymbol_exon8b'
#' @param display_coords `numeric` default NULL, with optional genome
#'    coordinates to use for the x-axis range.
#' @param facet_scales `character` string passed to
#'    `ggplot2::facet_wrap()` argument 'scales'. Use 'fixed' to share
#'    the same y-axis range in each panel.
#' @param minJunctionScore `numeric` default 10, passed to
#'    `prepareSashimi()` to require junctions with at least this score
#'    in order to be displayed. It is included here as a common
#'    argument to customize via `prepareSashimi()`.
#' @param scoreArcFactor `numeric` default 0.6 with default arc height
#'    added to the mean ribbon height, relative to its width.
#'    Default 0.6 ensures the ribbon rises at least 60% its width.
#'    Note argument 'scoreArcMinimum=100' can be adjusted through
#'    '...', and is not modified otherwise.
#' @param use_ylim `list` of `numeric` y-axis limits, applied to
#'    each 'sample_id' panel. It uses `ggh4x::scale_y_facet()` and
#'    is currently experimental, working out the enquo mechanics.
#' @param base_size `numeric` with base font size, default 18.
#' @param exonLabelSize `numeric` with exon label size, default 14.
#' @param label_junctions `logical` whether to label junctions with
#'    count/score, default TRUE.
#' @param junction_alpha `numeric` default 0.8, default alpha for
#'    junction ribbon arcs. Higher values are more opaque, and less
#'    transparent.
#' @param layout_ncol `integer` number of columns for output, default 1.
#' @param gene_panel_height `numeric` relative panel height for
#'    the gene/transcript/exon model, relative to each sashimi panel
#'    with relative height 1. Default is 1.
#' @param do_plot `logical` default TRUE, whether to render the plot.
#' @param use_memoise `logical` default FALSE, whether to enable
#'    data cache methods with memoise.
#' @param verbose `logical` whether to print verbose output, default TRUE.
#'    * `verbose == TRUE`: prints major steps and elapsed time (default)
#'    * `verbose == 2`: as above, adding `prepareSashimi()` major steps
#'    * `verbose == 3`: as above, adding `plotSashimi()` and `use_ylim`
#'    in detail. Probably too verbose at that point.
#' @param ... additional arguments are passed to relevant functions:
#'    `prepareSashimi()`, `gene2gg()`, `plotSashimi()`.
#'    Some common examples:
#'    * show: passed to `plotSashimi()` to customize which features
#'    are displayed: 'coverage', 'junctions', 'junctionLabels'
#'    * fill_scheme: passed to `plotSashimi()` to define coloring by
#'    'sample_id' (default) or 'exon'. When 'exon' the junctions are
#'    colored using 'junc_color' and 'junc_fill'.
#' 
#' @returns `list` with data elements:
#'    * cp_sashimi: the ggplot sashimi object
#'    * cp_sashimi_list: `list` of ggplot sashimi objects, per panel
#'    * cp_gene: the ggplot gene/transcript/exon object
#'    * cp_genes: `list` of ggplot gene/transcript/exon objects
#'    * gene_sd: the sashimi data from `prepareSashimi()`
#'    * layout_ncol: the layout ncol as used
#'    * use_ylim: the `list` of ylim as used
#'    * cp: the final assembled ggplot figure
#' 
#' @examples
#' sjenv <- sashimiDataConstants()
#' 
#' sjenv$filesDF <- farrisdata::farris_sashimi_files_df
#' sjfig <- splicejamFigure(sjenv=sjenv, gene="Gria1", sample_id=c("CA2_CB", "CA2_DE"))
#' 
#' sjfig <- splicejamFigure(sjenv=sjenv, gene="Gria1", sample_id=c("CA2_CB", "CA2_DE"),
#'    exon_range=c("Gria1_exon13", "Gria1_exon16"), scoreArcFactor=0.4)
#' 
#' sjfig2 <- splicejamFigure(sjenv=sjenv, gene="Gria1", sample_id=c("CA2_CB", "CA2_DE"),
#'    layout_ncol=2,
#'    exon_range=c("Gria1_exon13", "Gria1_exon16"), scoreArcFactor=0.4)
#' 
#' @export
splicejamFigure <- function
(sjenv=NULL,
 gene="Myom1",
 gene_sd=NULL,
 filesDF=NULL,
 sample_id=NULL,
 exon_range=NULL,
 display_coords=NULL,
 facet_scales="free_y",
 minJunctionScore=10,
 scoreArcFactor=0.6,
 use_ylim=NULL,
 base_size=18,
 exonLabelSize=14,
 label_junctions=TRUE,
 junction_alpha=0.8,
 layout_ncol=1,
 gene_panel_height=1,
 do_plot=TRUE,
 use_memoise=FALSE,
 panel_method=c(
    "patchwork",
    "cowplot"),
 verbose=TRUE,
 ...)
{
   # arglist
   arglist <- list(...)
   panel_method <- match.arg(panel_method);
   # exon_range
   # Todo: Make it work for multiple (neighboring) genes
   # - probably by creating new flattened exons
   if (length(exon_range) == 0) {
      # use overall first and last
      exon_range <- jamba::middle(n=2,
         names(sort(sjenv$flatExonsByGene[gene]@unlistData)))
   } else if (all(exon_range %in% names(sjenv$flatExonsByGene[gene]@unlistData))) {
      exon_range <- jamba::middle(n=2,
         names(sort(sjenv$flatExonsByGene[gene]@unlistData[exon_range])))
   } else {
      exon_range <- jamba::provigrep(exon_range,
         names(sjenv$flatExonsByGene[gene]@unlistData));
      if (length(exon_range) > 0) {
         exon_range <- jamba::middle(n=2,
            names(sort(sjenv$flatExonsByGene[gene]@unlistData[exon_range])))
      } else {
         exon_range <- jamba::middle(n=2,
            names(sort(sjenv$flatExonsByGene[gene]@unlistData)))   
      }
   }
   # genome coordinates
   if (length(display_coords) == 0) {
      display_coords <- range(as.data.frame(range(
         subset(sjenv$flatExonsByGene[gene]@unlistData,
            gene_nameExon %in% exon_range)
         ))[,c("start", "end")]);
   }

   # filesDF
   if (length(filesDF) > 0) {
      use_filesDF <- filesDF;
   } else {
      use_filesDF <- sjenv$filesDF;
   }
   if (length(use_filesDF) == 0 || nrow(use_filesDF) == 0) {
      stop("filesDF contains no data.");
   }

   # sample_id
   if (length(sample_id) == 0) {
      sample_id <- unique(use_filesDF$sample_id);
   } else {
      sample_id <- intersect(sample_id, use_filesDF$sample_id);
      if (length(sample_id) == 0) {
         stop("sample_id does not match filesDF$sample_id.");
      }
   }

   # prepareSashimi()
   st0 <- 0;
   st1 <- list(elapsed=0);
   if (length(gene_sd) == 0) {
      if (verbose) jamba::printDebug("splicejamFigure(): ", "prepareSashimi()");
      st1 <- system.time({
         gene_sd <- prepareSashimi(
            gene=gene,
            sample_id=sample_id,
            flatExonsByGene=sjenv$flatExonsByGene[gene],
            scoreArcFactor=scoreArcFactor,
            # flatExonsByTx=sjenv$flatExonsByTx,
            minJunctionScore=minJunctionScore,
            filesDF=use_filesDF,
            use_memoise=use_memoise,
            verbose=verbose > 1,
            ...)
      })
   } else {
      # Todo: Check output and subset for matching sample_id.
      use_sample_id <- sample_id;
      if (any(!gene_sd$df$sample_id %in% use_sample_id)) {
         if (verbose) jamba::printDebug("splicejamFigure(): ", "subset gene_sd by sample_id");
         gene_sd$df <- subset(gene_sd$df, sample_id %in% use_sample_id)
      }
   }
   asSeconds <- function(x, digits=2, ...){
      format(difftime(as.POSIXct(x), 0), digits=digits, ...)
   }
   if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st1["elapsed"]));
   
   # Compressed introns
   if ('ref2c' %in% names(gene_sd)) {
      ref2c <- gene_sd$ref2c;
   } else {
      ref2c <- sjenv$ref2c;
   }

   # layout_ncol <- 1;
   # base_size <- 18;
   xlim <- display_coords;
   ref_name <- as.character(head(
      GenomicRanges::seqnames(sjenv$flatExonsByGene[gene]@unlistData), 1));

   # plot panel heights, layout_ncol
   num_samples <- length(unique(sample_id))
   panel_height <- 1;
   # gene_panel_height <- 3;
   if (num_samples < layout_ncol) {
      layout_ncol <- num_samples;
   }
   plot_heights <- c(
      rep(panel_height, ceiling(num_samples / layout_ncol)),
      gene_panel_height)

   # plotSashimi()
   if (verbose) jamba::printDebug("splicejamFigure(): ", "plotSashimi()");
   use_show <- jamba::rmNA(c("coverage",
         "junction",
         ifelse(label_junctions, "junctionLabels", NA)))
   if ("show" %in% names(arglist)) {
      use_show <- arglist$show;
   }
   st2 <- list(elapsed=0);
   st2 <- system.time({
      gg_sashimi <- jamba::call_fn_ellipsis(plotSashimi,
         sashimi=gene_sd,
         show=use_show,
         color_sub=sjenv$color_sub,
         do_highlight=FALSE,
         facet_scales=facet_scales,
         junc_alpha=junction_alpha,
         label_coords=display_coords,
         # ref2c=ref2c,
         apply_facet=FALSE,
         # fill_scheme="sample_id",
         verbose=verbose > 2,
         ...);
   })
   if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st2["elapsed"]));
   
   # apply facet_wrap
   if (layout_ncol == 1) {
      strip.position <- "right";
   } else {
      strip.position <- "top";
   }
   gg_sashimi <- gg_sashimi +
      # ggplot2::facet_wrap(~sample_id,
      ggh4x::facet_wrap2(~sample_id,
         ncol=layout_ncol,
         remove_labels="y",
         strip.position=strip.position,
         scales=facet_scales);
      # Optional arguments:
      # strip.position 

   # detectedTx
   detectedTx <- NULL;
   if ('detectedTx' %in% names(arglist)) {
      detectedTx <- arglist[["detectedTx"]];
   } else if ('detectedTx' %in% names(sjenv)) {
      detectedTx <- sjenv$detectedTx;
   }
   gene_tx2geneDF <- subset(sjenv$tx2geneDF, gene_name %in% gene);
   detectedTx <- intersect(detectedTx,
      unique(gene_tx2geneDF$transcript_id))
   if (length(detectedTx) == 0) {
      if (verbose) {
         jamba::printDebug("splicejamFigure(): ",
            "using all transcript_id for 'gene'.");
      }
      detectedTx <- unique(gene_tx2geneDF$transcript_id);
   }

   ## Gene-exon model
   if (verbose) jamba::printDebug("splicejamFigure(): ", "gene2gg()");
   if (length(sjenv$flatExonsByTx) == 0) {
      txMatch <- 0;
   } else {
      txMatch <- match(detectedTx,
         names(sjenv$flatExonsByTx))
      if (any(is.na(txMatch))) {
         stop("Not all detectedTx were in names(sjenv$flatExonsByTx).");
      }
   }
   st3 <- list(elapsed=0);
   st3 <- system.time({
      gg_gene <- gene2gg(gene=gene,
         flatExonsByGene=sjenv$flatExonsByGene,
         flatExonsByTx=sjenv$flatExonsByTx[txMatch],
         label_coords=display_coords,
         ref2c=ref2c,
         # layout_ncol=layout_ncol,
         layout_ncol=1,
         exonLabelSize=exonLabelSize,
         ...);
   })
   if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st3["elapsed"]));

   ## Prepare for multi-panel figure
   if (verbose) jamba::printDebug("splicejamFigure(): ", paste(panel_method, "prep"));
   # apply x-axis limits
   cp_sashimi <- gg_sashimi +
      colorjam::theme_jam(base_size=base_size) +
      ggplot2::theme(axis.text.x=ggplot2::element_blank()) +
      ggplot2::xlab(NULL) +
      ggplot2::coord_cartesian(
         xlim=xlim,
         ylim=NULL);
   
   ## Apply multi-panel ylim
   panel_names <- names(table(gene_sd$df$sample_id))
   if (length(use_ylim) > 0) {
      if (!is.list(use_ylim)) {
         if (length(use_ylim) == 1) {
            use_ylim <- range(c(0, use_ylim))
         }
         use_ylim <- list(use_ylim);
         use_ylim <- rep(use_ylim, length.out=num_samples);
      }
      if (length(names(use_ylim)) > 0) {
         use_ylim <- use_ylim[panel_names];
         if (any(lengths(use_ylim) == 0)) {
            stop("use_ylim[sample_id] have empty entries.");
         }
      } else {
         use_ylim <- rep(use_ylim, length.out=num_samples);
         names(use_ylim) <- panel_names;
      }
   }

   if (FALSE && length(use_ylim) > 0) {
      ## Not active because it does not work well with facets,
      ## something about quoted values in the filter. Ugh.
      if (verbose) jamba::printDebug("splicejamFigure(): ", "use_ylim");
      if (!is.list(use_ylim)) {
         if (length(use_ylim) == 1) {
            use_ylim <- range(c(0, use_ylim))
         }
         use_ylim <- list(use_ylim);
         use_ylim <- rep(use_ylim, length.out=num_samples);
      }
      panel_names <- names(table(gene_sd$df$sample_id));
      if (length(names(use_ylim)) > 0) {
         use_ylim <- use_ylim[panel_names];
         if (any(lengths(use_ylim) == 0)) {
            stop("use_ylim[sample_id] have empty entries.");
         }
      } else {
         use_ylim <- rep(use_ylim, length.out=num_samples);
         names(use_ylim) <- panel_names;
      }
      
      # panel_names <- seq_len(num_samples);
      for (ipanel in panel_names) {
         ipanel1 <- match(ipanel, panel_names);
         if (verbose > 2) {
            jamba::printDebug("splicejamFigure(): ",
               "Applying ylim: ", use_ylim[[ipanel]],
               " to panel ", ipanel);
         }
         # not sure how to pass variable value
         # and not some enquosure thing.
         # enquo(variable)??
         suppressMessages({
            cp_sashimi <- cp_sashimi +
               ggh4x::scale_y_facet(
                  sample_id %in% ipanel | PANEL %in% ipanel,
                  limits=use_ylim[[ipanel1]])
         })
      }
   }

   # apply gene x-limits
   cp_gene <- gg_gene +
      colorjam::theme_jam(base_size=base_size) +
      ggplot2::ggtitle(NULL) +
      ggplot2::xlab(ref_name) +
      ggplot2::coord_cartesian(
         xlim=xlim);
   
   ## concise way to blank y-axis labels:
   # gg_plot + ggplot2::theme(axis.y = ggplot2::element_blank())
   blank_yaxis_labels <- function(...) {
      ggplot2::theme(
         # axis.ticks.y=ggplot2::element_blank(),
         axis.title.y=ggplot2::element_blank(),
         axis.text.y=ggplot2::element_blank()
      )
   }

   # split facet ggplot2 into list of ggplots
   split_gg_facets <- function
   (ggfacet, facetby="sample_id", ...)
   {
      facet_vals <- unique(ggfacet$data[[facetby]])
      # split facets into list
      lapply(jamba::nameVector(facet_vals), function(i){
         ggnew <- ggfacet;
         ggnew$data <- subset(ggnew$data, ggnew$data[[facetby]] %in% i);
         ggnew
      })
   }
   # split
   cp_sashimi_list <- split_gg_facets(cp_sashimi);
   cp_genes <- lapply(seq_len(layout_ncol), function(i){
      if (i == 1) {
         cp_gene
      } else {
         cp_gene + blank_yaxis_labels()
      }
   })

   # now apply y-axis limit to each panel
   if (length(use_ylim) > 0) {
      if (verbose) jamba::printDebug("splicejamFigure(): ", "use_ylim");
      panel_names <- names(table(gene_sd$df$sample_id))
      st4a <- system.time({
         for (ipanel in panel_names) {
            suppressMessages({
               cp_sashimi_list[[ipanel]] <- cp_sashimi_list[[ipanel]] +
                  ggplot2::scale_y_continuous(
                     limits=use_ylim[[ipanel]],
                     labels=scales::comma,
                     name="read depth")
            })
         }
      })
      if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st4a["elapsed"]));
   }

   # Create plot spacer in case needed
   cp_spacers <- NULL;
   if (layout_ncol > 1) {
      num_spacers <- (num_samples %% layout_ncol)
      cp_spacers <- lapply(seq_len(num_spacers), function(i){
         patchwork::plot_spacer()
      })
   }
   # Todo: fill blank panes when nsamples %% layout_ncol != 0
   layout_nrow <- ceiling(num_samples / layout_ncol) + 1;

   st4 <- list(elapsed=0);
   if (verbose) jamba::printDebug("splicejamFigure(): ", panel_method);
   st4 <- system.time({
      if ("cowplot" %in% panel_method) {
         # Cowplot multi-panel figure
         cp <- cowplot::plot_grid(
            plotlist=c(cp_sashimi_list,
               cp_spacers,
               cp_genes),
            ncol=layout_ncol,
            nrow=layout_nrow,
            align="v",
            axis="lr",
            rel_heights=plot_heights)
         if (isTRUE(do_plot)) plot(cp);
      } else {
         cp <- patchwork::wrap_plots(
            c(cp_sashimi_list,
               cp_spacers,
               cp_genes),
            ncol=layout_ncol,
            nrow=layout_nrow,
            heights=plot_heights) +
            patchwork::plot_layout(axis_titles="collect")
         #
         if (isTRUE(do_plot)) plot(cp);
      }
   })
   if (verbose) jamba::printDebug("", "elapsed ", indent=19, asSeconds(st4["elapsed"]));
   return(invisible(list(
      cp_sashimi=cp_sashimi,
      cp_sashimi_list=cp_sashimi_list,
      cp_gene=cp_gene,
      cp_genes=cp_genes,
      gene_sd=gene_sd,
      layout_ncol=layout_ncol,
      use_ylim=use_ylim,
      cp=cp
   )))
   }
