
#' Sashimi Shiny app server
#'
#' Sashimi Shiny app server
#'
#' This function contains the server logic for the R-shiny
#' Splicejam Sashimi viewer.
#'
#' The R-shiny app is started by `launchSashimiApp()`, which
#' calls `shiny::shinyApp()`, using arguments `server`, `ui`,
#' `onStart`, and `options`. This function fulfills the
#' argument `server`.
#'
#'
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#'
#' @family splicejam R-shiny functions
#'
#' @import jamba
#' @import dplyr
#' @import ggplot2
#' @import plotly
#' @import GenomicRanges
#'
#' @export
sashimiAppServer <- function
(input,
 output,
 session)
{
   #
   options("warn"=-1);
   output$sashimiplot_guide <- renderUI(sashimiplot_guide);
   output$sashimiplotviz_guide <- renderUI(sashimiplotviz_guide);

   ## show sessionInfo()
   output$sessionInfo <- renderPrint({
      capture.output(sessionInfo())
   })

   ## Update gene_choices
   get_gene_choices <- reactive({
      search_genelist <- input$search_genelist;
      if (length(search_genelist) > 0 && jamba::igrepHas("detected", search_genelist)) {
         gene_choices <- detectedGenes;
      } else {
         gene_choices <- jamba::mixedSort(unique(
            tx2geneDF$gene_name));
      }
      return(gene_choices);
   });

   get_default_gene <- function() {
      gene_choices <- tryCatch({
         get_gene_choices();
      }, error=function(e){
         detectedGenes;
      });
      jamba::printDebug("sashimiAppServer(): ",
         "length(gene_choices):",
         jamba::formatInt(length(gene_choices)));
      default_gene <- head(
         jamba::provigrep(c("Gria1",
            "Ntrk2",
            "Actb",
            "Gapd",
            "^[A-Z][a-z]{3}", "."),
            gene_choices),
         1);
      jamba::printDebug("sashimiAppServer(): ",
         c("default_gene:",
            default_gene), sep="");
      if (exists("gene")) {
         new_default_gene <- intersect(get("gene"),
            gene_choices);
         if (length(new_default_gene) > 0) {
            default_gene <- head(new_default_gene, 1);
         }
      }
      jamba::printDebug("sashimiAppServer(): ",
         c("default_gene:",
            default_gene), sep="");
      return(default_gene);
   }
   observe({
      search_genelist <- input$search_genelist;
      gene_choices <- get_gene_choices();
      default_gene <- get_default_gene();
      if (verbose) {
         jamba::printDebug("sashimiAppServer(): ",
            "updateSelectizeInput default_gene:",
            default_gene);
      }
      updateSelectizeInput(session,
         "gene",
         choices=gene_choices,
         selected=default_gene,
         server=TRUE);
   });

   ## server-side selectize gene list
   jamba::printDebug("sashimiAppServer(): ",
      "length(detectedGenes):",
      jamba::formatInt(length(detectedGenes)));
   default_gene <- get_default_gene();
   updateSelectizeInput(session,
      "gene",
      choices=detectedGenes,
      selected=default_gene,
      #choices=get_gene_choices(),
      #selected=head(
      #   jamba::provigrep(c("Gria1", "Ntrk2", "^[A-Z][a-z]{3}", "."),
      #      detectedGenes), 1),
      server=TRUE);

   output$gene_coords_label <- renderText({"Genome coordinate range"});

   # Define verbose flag
   if (!exists("verbose")) {
      verbose <- FALSE;
      jamba::printDebug("sashimiAppServer(): ",
         "verbose:", verbose);
   }

   # use debounce to slow down rendering a plot while changing parameters
   debounce_ms <- 1000;
   junction_alpha <- reactive({
      if (length(input$junction_alpha) == 0) {
         return(0.7);
      }
      input$junction_alpha;
   });
   junction_alpha_d <- debounce(junction_alpha, debounce_ms);
   share_y_axis <- reactive({
      if (length(input$share_y_axis) == 0) {
         return(TRUE);
      }
      return(input$share_y_axis);
   });
   share_y_axis_d <- debounce(share_y_axis, debounce_ms);
   show_gene_model <- reactive({
      if (length(input$show_gene_model) == 0) {
         return(TRUE)
      }
      input$show_gene_model
   });
   show_gene_model_d <- debounce(show_gene_model, debounce_ms);
   show_tx_model <- reactive({
      if (length(input$show_tx_model) == 0) {
         return(TRUE)
      }
      input$show_tx_model
   });
   show_tx_model_d <- debounce(show_tx_model, debounce_ms);
   show_detected_tx <- reactive({
      if (length(input$show_detected_tx) == 0) {
         return(TRUE)
      }
      input$show_detected_tx;
   });
   show_detected_tx_d <- debounce(show_detected_tx, debounce_ms);
   exon_label_size <- reactive({
      font_sizing <- input$exon_label_size;
      if (length(font_sizing) == 0) {
         font_sizing <- "Default"
      }
      base_font_size <- 1.0 * dplyr::case_when(
         igrepHas("-4", font_sizing) ~ -4,
         igrepHas("-3", font_sizing) ~ -3,
         igrepHas("-2", font_sizing) ~ -2,
         igrepHas("-1", font_sizing) ~ -1,
         igrepHas("Default", font_sizing) ~ 0,
         igrepHas("+1", font_sizing) ~ 1,
         igrepHas("+2", font_sizing) ~ 2,
         igrepHas("+3", font_sizing) ~ 3,
         igrepHas("+4", font_sizing) ~ 4,
         igrepHas(".", font_sizing) ~ 0
      );
      base_font_size;
   });
   exon_label_size_d <- debounce(exon_label_size, debounce_ms);
   panel_height <- reactive({
      if (length(input$panel_height) == 0) {
         return(200)
      }
      input$panel_height;
   });
   panel_height_d <- debounce(panel_height, debounce_ms);
   font_sizing <- reactive({
      font_sizing <- input$font_sizing;
      if (length(font_sizing) == 0) {
         font_sizing <- "Default"
      }
      base_font_size <- 1.0 * dplyr::case_when(
         igrepHas("-4", font_sizing) ~ -4,
         igrepHas("-3", font_sizing) ~ -3,
         igrepHas("-2", font_sizing) ~ -2,
         igrepHas("-1", font_sizing) ~ -1,
         igrepHas("Default", font_sizing) ~ 0,
         igrepHas("+1", font_sizing) ~ 1,
         igrepHas("+2", font_sizing) ~ 2,
         igrepHas("+3", font_sizing) ~ 3,
         igrepHas("+4", font_sizing) ~ 4,
         igrepHas(".", font_sizing) ~ 0
      );
      base_font_size * 1.5;
   });
   font_sizing_d <- debounce(font_sizing, debounce_ms);
   do_plotly <- reactive({
      if (length(input$do_plotly) == 0) {
         return(FALSE)
      }
      input$do_plotly;
   });
   do_plotly_d <- debounce(do_plotly, debounce_ms);
   enable_highlights <- reactive({
      if (length(input$enable_highlights) == 0) {
         return(FALSE)
      }
      input$enable_highlights;
   });
   enable_highlights_d <- debounce(enable_highlights, debounce_ms);


   # update the "Update" button when something has changed
   observe({
      gene <- input$gene;
      sample_order <- input$selectionto_order;
      min_junction_reads <- input$min_junction_reads;
      include_strand <- input$include_strand;
      layout_ncol <- input$layout_ncol;
      exon_range <- input$exon_range;
      gene_coords <- input$gene_coords;
      shinyjs::enable("calc_gene_params");
   });

   # Update the slider bar for each gene, or when slider type is changed
   observe({
      gene <- input$gene;
      if (verbose) {
         jamba::printDebug("sashimiAppServer(): ",
            "observe input$gene");
         jamba::printDebug("sashimiAppServer(): ",
            "updateInputSlider gene:",
            gene);
      }
      if (length(gene) > 0 && nchar(gene) > 0) {
         ## Handle "All Genes" where it is not present in flatExonsByGene
         flatExonsByGene1 <- get_flat_gene_exons();
         jamba::printDebug("sashimiAppServer(): ",
            "update gene slider bar, gene:", gene,
            ", length(flatExonsByGene1):", length(flatExonsByGene1),
            ", length(flatExonsByGene1[[gene]]):", length(flatExonsByGene1[[gene]]));
         chr_range <- as.data.frame(range(flatExonsByGene1[[gene]]))[,c("start", "end")];
         coords_label <- paste0("Coordinate range for ",
            gene,
            " on ",
            as.character(seqnames(head(flatExonsByGene1[[gene]], 1))),
            " (", as.character(strand(head(flatExonsByGene1[[gene]], 1))),
            "):");
         if ("gene_nameExon" %in% colnames(values(flatExonsByGene1[[gene]]))) {
            exon_names <- jamba::mixedSort(
               values(flatExonsByGene1[[gene]])$gene_nameExon);
            jamba::printDebug("sashimiAppServer(): ",
               "exon_names:",
               exon_names);
         } else {
            exon_names <- NULL;
         }
         if (length(chr_range) > 0) {
            ## Update slider text label
            output$gene_coords_label <- renderText({
               coords_label
            });
            ## Update sliderInput for gene_coords
            updateSliderInput(session,
               "gene_coords",
               min=min(chr_range[["start"]]),
               max=max(chr_range[["end"]]),
               value=range(c(chr_range[["start"]], chr_range[["end"]]))
            );
         }
         ## update exon name range slider
         if (length(exon_names) > 0) {
            updateSliderTextInput(session,
               "exon_range",
               choices=exon_names,
               selected=c(
                  head(exon_names, 1),
                  tail(exon_names, 1))
            );
         }
      }
   });

   get_active_gene <- reactive({
      input$calc_gene_params;
      gene <- isolate(input$gene);
      gene;
   });

   get_flat_gene_exons <- reactive({
      gene <- input$gene;
      if (length(gene) > 0 && nchar(gene) > 0) {
         ## Handle "All Genes" where it is not present in flatExonsByGene
         if (!gene %in% names(flatExonsByGene)) {
            jamba::printDebug("sashimiAppServer(): ",
               "Creating flat exons for gene:",
               gene);
            flatExonsByGene1 <- flattenExonsBy(genes=gene,
               by="gene",
               exonsByTx=exonsByTx,
               cdsByTx=cdsByTx,
               tx2geneDF=tx2geneDF);
         } else {
            jamba::printDebug("sashimiAppServer(): ",
               "Re-using flat exons for gene:",
               gene);
            flatExonsByGene1 <- flatExonsByGene[gene];
         }
      } else {
         flatExonsByGene1 <- flatExonsByGene;
      }
      return(flatExonsByGene1);
   });

   get_flat_gene_exons_plot <- reactive({
      input$calc_gene_params;
      gene <- isolate(input$gene);
      if (length(gene) > 0 && nchar(gene) > 0) {
         ## Handle "All Genes" where it is not present in flatExonsByGene
         if (!gene %in% names(flatExonsByGene)) {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Creating flat exons for gene:",
                  gene);
            }
            flatExonsByGene1 <- flattenExonsBy(genes=gene,
               by="gene",
               exonsByTx=exonsByTx,
               cdsByTx=cdsByTx,
               tx2geneDF=tx2geneDF);
         } else {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Re-using flat exons for gene:",
                  gene);
            }
            flatExonsByGene1 <- flatExonsByGene[gene];
         }
      } else {
         flatExonsByGene1 <- flatExonsByGene;
      }
      return(flatExonsByGene1);
   });

   get_gene_coords <- reactive({
      input$calc_gene_params;
      gene <- isolate(input$gene);
      ## get gene coordinate range
      use_exon_names <- isolate(input$use_exon_names);
      flatExonsByGene1 <- get_flat_gene_exons_plot();
      if ("exon names" %in% use_exon_names) {
         exon_range <- isolate(input$exon_range);
         # Convert exon names to coordinates
         gene_coords <- range(as.data.frame(range(
            subset(flatExonsByGene1[[gene]], gene_nameExon %in% exon_range)
         ))[,c("start", "end")]);
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "gene_coords inferred from input$exon_range:",
               gene_coords);
         }
      } else {
         gene_coords <- isolate(input$gene_coords);
         if (length(gene_coords) == 0 || all(gene_coords %in% c(28,117))) {
            # default has not been initialized yet, so take full gene range
            gene_coords <- range(as.data.frame(range(
               flatExonsByGene1[[gene]]
            ))[,c("start", "end")]);
         }
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "gene_coords from input$gene_coords:",
               gene_coords);
         }
      }
      gene_coords;
   })

   get_sample_id <- reactive({
      input$calc_gene_params;
      ## Define sample_id from sample selection
      sample_order <- isolate(input$selectionto_order);
      if (!exists("sample_id")) {
         if (length(sample_order) == 0) {
            sample_id <- head(unique(filesDF$sample_id), 4);
         } else {
            sample_id <- sample_order;
         }
      } else if (length(sample_order) > 0) {
         sample_id <- sample_order;
      }
      if (verbose) {
         jamba::printDebug("sashimiAppServer(): ",
            c("get_sample_id('",
               sample_id,
               "')"),
            sep="");
      }
      sample_id;
   });

   get_sashimi_data <- reactive({
      input$calc_gene_params;
      shinyjs::disable("calc_gene_params");
      #gene <- get_active_gene();
      gene <- isolate(input$gene);
      flatExonsByGene1 <- get_flat_gene_exons_plot();
      if (length(gene) == 0) {
         return(NULL);
      }
      if (!exists("flatExonsByGene") ||
            !exists("filesDF")) {
         return(NULL);
      }

      ## Wrap the workflow in a progress bar
      if (!exists("prepareSashimi_m")) {
         prepareSashimi_m <- memoise::memoise(prepareSashimi,
            cache=memoise::cache_filesystem("sashimidata_memoise"));
      }

      ## Define sample_id from sample selection
      #sample_id <- get_sample_id();
      sample_id <- get_sample_id_dt();
      if (verbose) {
         jamba::printDebug("sashimiAppServer(): ",
            "Using sample_id:",
            sample_id);
         jamba::printDebug("sashimiAppServer(): ",
            "Using gene:",
            gene);
      }

      min_junction_reads <- isolate(input$min_junction_reads);
      include_strand <- isolate(input$include_strand);
      if (!exists("use_memoise")) {
         use_memoise <- TRUE;
      }
      withProgress(
         message="Preparing Sashimi data.",
         value=0,
         {
            if (1 || verbose) {
               gene_has_cache <- memoise::has_cache(prepareSashimi_m)(
                  gene=gene,
                  flatExonsByGene=flatExonsByGene1,
                  minJunctionScore=min_junction_reads,
                  sample_id=sample_id,
                  filesDF=filesDF,
                  include_strand=include_strand,
                  verbose=verbose,
                  use_memoise=use_memoise,
                  do_shiny_progress=TRUE);
               jamba::printDebug("sashimiAppServer(): ",
                  "gene:",
                  gene,
                  ", gene_has_cache:",
                  gene_has_cache);
            }
            sashimi_data <- prepareSashimi_m(
               gene=gene,
               flatExonsByGene=flatExonsByGene1,
               minJunctionScore=min_junction_reads,
               sample_id=sample_id,
               filesDF=filesDF,
               include_strand=include_strand,
               verbose=verbose,
               use_memoise=use_memoise,
               do_shiny_progress=TRUE);
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "sdim(sashimi_data):");
               print(jamba::sdim(sashimi_data));
            }
         }
      )
      ## Check for missing data
      some_null <- sashimi_data$some_null;
      if (length(some_null) && some_null) {
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "sashimi_data some_null:",
               some_null);
         }
         withProgress(
            message="Repeating Sashimi steps",
            value=0,
            {
               ## Some underlying data was NULL, therefore try to repair
               if (compareVersion(as.character(packageVersion("memoise")), "1.1.0.900") >= 0) {
                  ## memoise 1.1.0.900 has new function memoise::drop_cache()
                  if (1 || verbose) {
                     jamba::printDebug("sashimiAppServer(): ",
                        "Dropping memoise cache, then re-creating with memoise.");
                  }
                  memoise::drop_cache(prepareSashimi_m)(
                     gene=gene,
                     flatExonsByGene=flatExonsByGene1,
                     minJunctionScore=min_junction_reads,
                     sample_id=sample_id,
                     filesDF=filesDF,
                     include_strand=include_strand,
                     verbose=verbose,
                     use_memoise=use_memoise,
                     do_shiny_progress=TRUE);
                  ## Re-run prepareSashimi using memoise so it will
                  ## create a new cache file
                  sashimi_data <- prepareSashimi_m(
                     gene=gene,
                     flatExonsByGene=flatExonsByGene1,
                     minJunctionScore=min_junction_reads,
                     sample_id=sample_id,
                     filesDF=filesDF,
                     include_strand=include_strand,
                     verbose=verbose,
                     use_memoise=use_memoise,
                     do_shiny_progress=TRUE);
               } else {
                  ## Re-run prepareSashimi without using memoise
                  if (1 || verbose) {
                     jamba::printDebug("sashimiAppServer(): ",
                        "Invalid memoise cache file, re-running prepareSashimi() without memoise.",
                        fgText="red");
                  }
                  sashimi_data <- prepareSashimi(
                     gene=gene,
                     flatExonsByGene=flatExonsByGene1,
                     minJunctionScore=min_junction_reads,
                     sample_id=sample_id,
                     filesDF=filesDF,
                     include_strand=include_strand,
                     verbose=verbose,
                     use_memoise=use_memoise,
                     do_shiny_progress=TRUE);
               }
            }
         );
      }
      sashimi_data;
   });

   output$sashimiplot_output <- renderUI({
      sashimi_data <- get_sashimi_data();
      if ("ref2c" %in% names(attributes(sashimi_data))) {
         ref2c <- attr(sashimi_data, "ref2c");
      } else if ("ref2c" %in% names(sashimi_data)) {
         ref2c <- sashimi_data$ref2c;
      } else {
         ref2c <- NULL;
      }

      flatExonsByGene1 <- get_flat_gene_exons_plot();
      gene <- get_active_gene();

      # Update in parent environment
      #sashimi_data <- sashimi_data;
      if (length(sashimi_data) == 0) {
         ## Todo: make empty plot a minimum height in pixels
         jamba::printDebug("sashimiAppServer(): ",
            "Rendering blank panel for ",
            "sashimiplot_output");
         tagList(
            renderPlot(
               height=300,
               ggplot() + geom_blank()
            )
         );
      } else {
         if (share_y_axis_d()) {
            facet_scales <- "fixed";
         } else {
            facet_scales <- "free_y";
         }
         ## do_plotly <- do_plotly_d()
         #do_highlight <- (enable_highlights_d() && do_plotly_d());
         do_highlight <- (do_plotly_d());

         ## Obtain the baseline ggplot object
         gg_sashimi <- plotSashimi(sashimi_data,
            show=c("coverage", "junction", "junctionLabels"),
            color_sub=color_sub,
            do_highlight=do_highlight,
            facet_scales=facet_scales,
            junc_alpha=junction_alpha_d(),
            label_coords=get_gene_coords(),
            ref2c=ref2c,
            fill_scheme="sample_id");

         # Check layout_ncol
         layout_ncol <- isolate(input$layout_ncol);
         if (length(layout_ncol) > 0 && layout_ncol > 1) {
            # change from facet_grid to facet_wrap()
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Applying facet_wrap() with ncol:",
                  layout_ncol);
            }
            gg_sashimi <- gg_sashimi +
               facet_wrap(~sample_id,
                  ncol=layout_ncol,
                  scales=facet_scales);
         }
         #gg_sashimi <<- gg_sashimi;
         ## Optionally prepare gene-exon model
         if (show_gene_model_d()) {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Getting gene model with label_coords:",
                  get_gene_coords());
            }
            if (show_tx_model_d() && length(flatExonsByTx) > 0) {
               if (show_detected_tx_d()) {
                  jamba::printDebug("preparing to show detected transcripts:",
                     length(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]));
                  if (length(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]) == 0) {
                     jamba::printDebug("No transcripts in flatExonsByTx matched detectedTx. Showing head(names(flatExonsByTx)), then head(flatExonsByTx):");
                     print(head(names(flatExonsByTx)));
                     print(head(flatExonsByTx));
                  } else {
                     print(head(subset(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]@unlistData, gene_name %in% gene)));
                  }
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene1,
                     flatExonsByTx=flatExonsByTx[names(flatExonsByTx) %in% detectedTx],
                     label_coords=get_gene_coords(),
                     ref2c=ref2c,
                     exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
               } else {
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene1,
                     flatExonsByTx=flatExonsByTx,
                     label_coords=get_gene_coords(),
                     ref2c=ref2c,
                     exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
               }
            } else {
               gg_gene <- gene2gg(gene=gene,
                  flatExonsByGene=flatExonsByGene1,
                  label_coords=get_gene_coords(),
                  ref2c=ref2c,
                  exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
            }
            gg_gene <- gg_gene;
         }

         ## Prepare text label with genome coordinates
         ref_name <- as.character(head(seqnames(flatExonsByGene1[[gene]]), 1));
         coord_label <- paste0(
            ref_name,
            ":",
            paste(get_gene_coords(), collapse="-")
         );
         output$sashimitext_output <- renderText({
            HTML(paste0("Region ",
               coord_label))
         });

         ## Prepare plotly or ggplot2 output
         #sample_id <- unique(filesDF$sample_id);
         #num_samples <- max(c(length(sample_id), 1)) + 1;
         #num_samples <- 2;
         ## Define sample_id from sample selection
         sample_id <- get_sample_id_dt();

         num_samples <- length(unique(sample_id));
         panel_height <- panel_height_d();
         # adjust base_size for fonts using exponent based upon panel height
         font_exp <- 1/3;
         base_font_size <- 16 + 1.0 * font_sizing_d();
         base_size <- base_font_size * (panel_height^font_exp)/(250^font_exp);
         plot_height <- panel_height * (num_samples + show_gene_model_d());
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "num_samples:", num_samples,
               ", panel_height:", panel_height,
               ", base_size:", format(digits=1, base_size),
               ", base_font_size:", base_font_size,
               ", layout_ncol:", layout_ncol,
               ", plot_height:", plot_height);
         }
         if (do_plotly_d()) {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Preparing plotly output.");
            }
            if (show_gene_model_d()) {
               ## use plotly, showing gene model
               ggly1 <- plotly::ggplotly(
                  gg_sashimi +
                     colorjam::theme_jam(base_size=base_size) +
                     theme(axis.text.x=element_blank()) +
                     xlab(NULL) +
                     coord_cartesian(xlim=get_gene_coords()),
                  tooltip="text")# %>%
                  #plotly::style(
                  #   hoveron="fill"
                  #);
               ggly2 <- plotly::ggplotly(
                  gg_gene +
                     colorjam::theme_jam(base_size=base_size) +
                     ggtitle(NULL) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=get_gene_coords()),
                  tooltip="text");
               gg_ly <- suppressMessages(
                  plotly::subplot(
                     ggly1,
                     ggly2,
                     nrows=2,
                     heights=c(num_samples, 2 + show_tx_model_d())/(num_samples + 2 + show_tx_model_d()),
                     shareX=TRUE
                  ) %>% layout(height=plot_height)
               );
               if (enable_highlights_d()) {
                  gg_ly <- gg_ly %>%
                     plotly::highlight("plotly_hover",
                        opacityDim=0.8,
                        selected=attrs_selected(
                           line=list(color="#444444")));
               }
            } else {
               ## use plotly, showing gene model
               gg_ly <- plotly::ggplotly(
                  gg_sashimi +
                     colorjam::theme_jam(base_size=base_size) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=get_gene_coords()),
                  tooltip="text",
                  height=plot_height
               )
               # %>% style(hoveron="fill");
               if (enable_highlights_d()) {
                  gg_ly <- gg_ly %>%
                     highlight("plotly_hover",
                        opacityDim=0.8,
                        selected=attrs_selected(
                           line=list(color="#444444")));
               }
            }
            ## Remove the color legend (again)
            if (!input$plotly_legend) {
               gg_ly <- gg_ly %>%
                  layout(showlegend=FALSE);
            }
            # Try converting to plotlyOutput to define a fixed plot height
            output$plotly <- renderPlotly({
               gg_ly %>%
                  layout(margin=list(r=100, l=70, t=20, b=70))
            });
            plotly_out <- plotlyOutput(
               "plotly",
               height=plot_height
            );
            return(plotly_out);
         } else {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Preparing ggplot output.");
            }
            ## Non-plotly static plot output
            if (length(get_gene_coords()) > 0) {
               if (show_gene_model_d()) {
                  ## With gene-transcript-exon model
                  tagList(renderPlot(
                     height=plot_height,
                     suppressMessages(
                        cowplot::plot_grid(
                           gg_sashimi +
                              colorjam::theme_jam(base_size=base_size) +
                              theme(axis.text.x=element_blank()) +
                              xlab(NULL) +
                              coord_cartesian(xlim=get_gene_coords()),
                           gg_gene +
                              colorjam::theme_jam(base_size=base_size) +
                              ggtitle(NULL) +
                              xlab(ref_name) +
                              coord_cartesian(xlim=get_gene_coords()),
                           ncol=1,
                           align="v",
                           axis="lr",
                           rel_heights=c(num_samples, 2+show_tx_model_d()*1)
                        )
                     )
                  ));
               } else {
                  ## No gene-transcript-exon model
                  tagList(renderPlot(
                     height=plot_height,
                     suppressMessages(
                        gg_sashimi +
                           colorjam::theme_jam(base_size=base_size) +
                           xlab(ref_name) +
                           coord_cartesian(xlim=get_gene_coords())
                     )
                  ));
               }
            } else {
               ## Fallback when get_gene_coords() is empty
               ## sometimes happens at startup, unclear when/why
               tagList(renderPlot(
                  height=plot_height,
                  gg_sashimi +
                     colorjam::theme_jam(base_size=base_size)
               ));
            }
         }
      }
   });

   # Helper function to colorize a datatable
   colorize_DT <- function
   (dt,
    df,
    color_sub,
    ...)
   {
      # check each colname for matching colors
      for (iCol in colnames(df)) {
         if (all(unique(df[[iCol]]) %in% names(color_sub))) {
            dt <- dt %>%
               DT::formatStyle(
                  iCol,
                  backgroundColor=DT::styleEqual(
                     levels=names(color_sub),
                     values=rgb2col(col2rgb(color_sub))
                  ),
                  color=DT::styleEqual(
                     levels=names(color_sub),
                     values=setTextContrastColor(color_sub)
                  )
               );
         }
      }
      dt;
   }

   # filesDF table output
   output$files_df <- DT::renderDT({

      # add new colname to contain checkboxes
      shinyCheckboxInput <- function(FUN, len, id, ...) {
         inputs <- sapply(seq_len(len), function(i){
            as.character(
               FUN(
                  paste0(id, i),
                  label="selected?",
                  width="30px",
                  ...
               )
            )
         });
         inputs;
      }
      files_dt <- filesDF;
      #files_dt$selected <- shinyCheckboxInput(
      #   checkboxInput,
      #   nrow(files_dt),
      #   "cbox_");

      files_dt <- files_dt %>%
         dplyr::select(sample_id,
            everything(),
            -scale_factor, -type, -url,
            type, scale_factor, url) %>%
         dplyr::arrange(sample_id, type) %>%
         DT::datatable(
            editable=TRUE,
            rownames=FALSE,
            escape=TRUE,
            #height='15px',
            extensions=c("RowReorder"),
            selection="none",
            options=list(
               autoWidth=TRUE,
               deferRender=TRUE,
               pageLength=24,
               lengthMenu=c(12,24,48,120),
               rowReorder=TRUE
            )
         );
         # optionally colorize sample_id using color_sub
         if (!all(c("junction", "bw") %in% names(color_sub))) {
            color_sub["junction"] <- "slateblue1";
            color_sub["bw"] <- "darkslategray1";
            # update color_sub in the parent environment
            color_sub <<- color_sub;
         }

         # check each colname for matching colors
         files_dt <- colorize_DT(dt=files_dt,
            df=filesDF,
            color_sub=color_sub);

         files_dt;
      },
      options=list(
         lengthChange=FALSE,
         pageLength=24
      )
   );

   # Render samples_df for sample selection and ordering
   data <- reactiveValues(
      samples_data=unique(filesDF[,setdiff(colnames(filesDF), c("url", "type", "scale_factor")), drop=FALSE]),
      sample_id=if (exists("sample_id")) {
         sample_id
      } else if (all(c("CA1_CB", "CA1_DE", "CA2_CB", "CA2_DE") %in% filesDF$sample_id)) {
         c("CA1_CB", "CA1_DE", "CA2_CB", "CA2_DE")
      } else if (is.factor(filesDF$sample_id)) {
         levels(filesDF$sample_id)
      } else {
         unique(filesDF$sample_id)
      }
   );
   get_sample_id_dt <- reactive({
      # force update when Calculate Gene Params button is clicked
      # This step prevents updates from changing the sashimi plot
      # until the Calculate button is clicked.
      input$calc_gene_params;
      new_sample_id <- isolate(data$sample_id);
      if (verbose >= 0) {
         jamba::printDebug("get_sample_id_dt(): ",
            new_sample_id);
      }
      return(new_sample_id)
   })
   output$samplesdf <- DT::renderDataTable({
      isolate(data$samples_data)
   },
      options=list(
         pageLength=nrow(isolate(data$samples_data)),
         lengthMenu=nrow(isolate(data$samples_data)),
         dom='t',
         columnDefs=list(
            list(
               targets=seq_len(ncol(isolate(data$samples_data))),
               searchable=FALSE,
               sortable=FALSE))
      ),
      selection=list(
         mode='multiple',
         selected=match(get_sample_id_dt(),
            isolate(data$samples_data)$sample_id)
      )
   );
   # event of clicking the samples_df table
   observeEvent(input$samplesdf_cell_clicked, {
      selected_rows <- input$samplesdf_rows_selected;
      if (verbose > 1) {
         jamba::printDebug("selected_rows:", selected_rows);
      }
      # calculate new row order
      row_order <- order(
         seq_along(isolate(data$samples_data)$sample_id) %in% selected_rows,
         decreasing=TRUE
      )
      new_sample_id <- isolate(data$samples_data)$sample_id[selected_rows];
      if (verbose > 1) {
         jamba::printDebug("New sample_id values selected:",
            new_sample_id);
      }
      data$sample_id <- new_sample_id;
      data$samples_data <- isolate(data$samples_data)[row_order,,drop=FALSE];
      proxy <- DT::dataTableProxy('samplesdf')
      DT::replaceData(proxy, data$samples_data)
      # make sure to select the rows again
      DT::selectRows(proxy, seq_along(selected_rows))
      # enable the Calculate button so the sashimi plot can be updated
      shinyjs::enable("calc_gene_params");
   })

}
