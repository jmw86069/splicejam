
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
#' @importFrom plotly subplot
#' @import GenomicRanges
#' @importFrom magrittr %>%
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

   if (exists("verbose")) {
      verbose <- get("verbose");
   } else {
      verbose <- FALSE;
   }
   if (exists("dodebug")) {
      dodebug <- get("dodebug");
   } else {
      dodebug <- FALSE;
   }

   ##  gene_choices
   get_gene_choices <- reactive({
      search_genelist <- input$search_genelist;
      if (length(search_genelist) > 0 && jamba::igrepHas("detected", search_genelist)) {
         gene_choices <- detectedGenes;
      } else {
         gene_choices <- jamba::mixedSort(unique(
            tx2geneDF$gene_name));
      }
      gene_choices <- unique(c("blank", gene_choices));
      return(gene_choices);
   });

   get_default_gene <- function() {
      # return default_gene if defined
      if (exists("default_gene") && length(default_gene) > 0 && all(nchar(default_gene) > 0)) {
         default_gene <- head(default_gene[nchar(default_gene) > 0], 1);
         return(default_gene);
      }
      if (exists("detectedGenes")) {
         gene_choices <- detectedGenes;
      } else {
         gene_choices <- mixedSort(unique(tx2geneDF$gene_name));
      }
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
         if (length(new_default_gene) > 0 && !new_default_gene == default_gene) {
            default_gene <- head(new_default_gene, 1);
            jamba::printDebug("sashimiAppServer(): ",
               c("new_default_gene:",
                  default_gene), sep="");
         }
      }
      return(default_gene);
   }

   observe({
      #search_genelist <- input$search_genelist;
      gene_choices <- get_gene_choices();
      default_gene <- get_default_gene();

      jamba::printDebug("sashimiAppServer(): ",
         "updateSelectizeInput default_gene:",
         default_gene);

      shiny::updateSelectizeInput(session,
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
      choices=unique(c("blank", detectedGenes)),
      selected=default_gene,
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
         jamba::igrepHas("-4", font_sizing) ~ -4,
         jamba::igrepHas("-3", font_sizing) ~ -3,
         jamba::igrepHas("-2", font_sizing) ~ -2,
         jamba::igrepHas("-1", font_sizing) ~ -1,
         jamba::igrepHas("Default", font_sizing) ~ 0,
         jamba::igrepHas("+1", font_sizing) ~ 1,
         jamba::igrepHas("+2", font_sizing) ~ 2,
         jamba::igrepHas("+3", font_sizing) ~ 3,
         jamba::igrepHas("+4", font_sizing) ~ 4,
         jamba::igrepHas(".", font_sizing) ~ 0
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
   gene_panel_height <- reactive({
      if (length(input$gene_panel_height) == 0) {
         return(200)
      }
      input$gene_panel_height;
   });
   gene_panel_height_d <- debounce(gene_panel_height, debounce_ms);
   font_sizing <- reactive({
      font_sizing <- input$font_sizing;
      if (length(font_sizing) == 0) {
         font_sizing <- "Default"
      }
      base_font_size <- 1.0 * dplyr::case_when(
         jamba::igrepHas("-4", font_sizing) ~ -4,
         jamba::igrepHas("-3", font_sizing) ~ -3,
         jamba::igrepHas("-2", font_sizing) ~ -2,
         jamba::igrepHas("-1", font_sizing) ~ -1,
         jamba::igrepHas("Default", font_sizing) ~ 0,
         jamba::igrepHas("+1", font_sizing) ~ 1,
         jamba::igrepHas("+2", font_sizing) ~ 2,
         jamba::igrepHas("+3", font_sizing) ~ 3,
         jamba::igrepHas("+4", font_sizing) ~ 4,
         jamba::igrepHas(".", font_sizing) ~ 0
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
      if (!"blank" %in% gene) {
         shinyjs::enable("calc_gene_params");
      }
   });

   get_gene_coords_label <- reactive({
      flat_list <- get_flat_gene_exons();
      gene <- flat_list$gene;
      flatExonsByGene1 <- flat_list$flatExonsByGene;
      if (length(flatExonsByGene1) == 0) {
         coords_label <- "blank (no gene is displayed)";
      } else {
         chr_range <- as.data.frame(range(flatExonsByGene1[[gene]]))[,c("start", "end")];
         coords_label <- paste0("Coordinate range for ",
            gene,
            " on ",
            as.character(seqnames(head(flatExonsByGene1[[gene]], 1))),
            " (", as.character(strand(head(flatExonsByGene1[[gene]], 1))),
            "):");
      }
      coords_label;
   })

   output$gene_coords_label <- renderText({
      get_gene_coords_label()
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
         if (!"blank" %in% gene) {
            shinyjs::enable("calc_gene_params");
            shinyjs::enable("exon_range");
            shinyjs::enable("gene_coords");
            ## Handle "All Genes" where it is not present in flatExonsByGene
            flat_list <- get_flat_gene_exons();
            flatExonsByGene1 <- flat_list$flatExonsByGene;
            if (verbose > 1) {
               jamba::printDebug("sashimiAppServer(): ",
                  "update gene slider bar, gene:", gene,
                  ", length(flatExonsByGene1):", length(flatExonsByGene1),
                  ", length(flatExonsByGene1[[gene]]):", length(flatExonsByGene1[[gene]]));
            }
            chr_range <- as.data.frame(range(flatExonsByGene1[[gene]]))[,c("start", "end")];

            ## update exon name range slider
            if ("gene_nameExon" %in% colnames(values(flatExonsByGene1[[gene]]))) {
               exon_names <- jamba::mixedSort(
                  values(flatExonsByGene1[[gene]])$gene_nameExon);
               jamba::printDebug("sashimiAppServer(): ",
                  "exon_names:",
                  exon_names);
               shinyWidgets::updateSliderTextInput(session,
                  "exon_range",
                  choices=exon_names,
                  selected=c(
                     head(exon_names, 1),
                     tail(exon_names, 1))
               );
            } else {
               exon_names <- NULL;
               shinyjs::disable("exon_range");
            }
            if (length(chr_range) > 0) {
               ## Update sliderInput for gene_coords
               shiny::updateSliderInput(session,
                  "gene_coords",
                  min=min(chr_range[["start"]]),
                  max=max(chr_range[["end"]]),
                  value=range(c(chr_range[["start"]], chr_range[["end"]]))
               );
            } else {
               shinyjs::disable("gene_coords");
            }
         } else {
            # Update gene coord label and sliders for "blank" gene
            shinyjs::disable("calc_gene_params");
            shinyjs::disable("exon_range");
            shinyjs::disable("gene_coords");
            shinyWidgets::updateSliderTextInput(session,
               "gene_coords",
               choices=c("1", "2"),
               selected=c("1", "2"));
            shinyWidgets::updateSliderTextInput(session,
               "exon_range",
               choices=c("1", "2"),
               selected=c("1", "2"));
         }
      }
   });

   # assemble all gene query values into a list when "Update" button pressed
   get_active_gene <- reactive({
      input$calc_gene_params;
      # get isolated variable values
      gene <- isolate(input$gene);
      use_exon_names <- isolate(input$use_exon_names);
      gene_coords <- isolate(input$gene_coords);
      exon_range <- isolate(input$exon_range);
      min_junction_reads <- isolate(input$min_junction_reads);
      include_strand <- isolate(input$include_strand);
      gene_vals <- list(gene=gene,
         use_exon_names=use_exon_names,
         gene_coords=gene_coords,
         exon_range=exon_range,
         min_junction_reads=min_junction_reads,
         include_strand=include_strand);
      # disable the button once values are updated
      shinyjs::disable("calc_gene_params");
      gene_vals;
   });

   # This function reacts to changes in input$gene, the "Select Gene" widget
   get_flat_gene_exons <- reactive({
      gene <- input$gene;
      if (length(gene) > 0 && nchar(gene) > 0) {
         if ("blank" %in% gene) {
            return(NULL)
         }
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
         return(flatExonsByGene);
      }
      return(list(flatExonsByGene=flatExonsByGene1,
         gene=gene));
   });

   # Same function as get_flat_gene_exons() except it does not react to changes in input$gene
   # This function only reacts to "Update Sashimi Plots" input$calc_gene_params
   get_flat_gene_exons_plot <- reactive({
      gene_vals <- get_active_gene();
      gene <- gene_vals$gene;
      if (length(gene) > 0 && nchar(gene) > 0) {
         if ("blank" %in% gene) {
            return(NULL)
         }
         ## Handle "All Genes" where it is not present in flatExonsByGene
         if (!gene %in% names(flatExonsByGene)) {
            if (verbose) {
               jamba::printDebug("sashimiAppServer(): ",
                  "Creating flat exons for gene: ",
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
                  "Re-using flat exons for gene: ",
                  gene);
            }
            flatExonsByGene1 <- flatExonsByGene[gene];
         }
      } else {
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "No gene specified, returning full flatExonsByGene, gene: ",
               gene);
         }
         flatExonsByGene1 <- flatExonsByGene;
      }
      return(list(flatExonsByGene=flatExonsByGene1,
         gene_vals=gene_vals));
      #return(flatExonsByGene1);
   });

   get_display_coords <- reactive({
      ## 27jul2021: can this input$calc_gene_params be ignored?
      ## Instead rely upon get_flat_gene_exons_plot()?
      ## This change could force order of operation so these are not
      ## updating at the same time and out of sync.
      ##
      ## I am not sure how and when isolate(input$gene) gets updated
      #input$calc_gene_params;

      # move this reactive before the isolate() statements
      flat_list <- get_flat_gene_exons_plot();
      if (length(flat_list) == 0) {
         return(NULL)
      }
      flatExonsByGene1 <- flat_list$flatExonsByGene;
      gene_vals <- flat_list$gene_vals;
      gene <- gene_vals$gene;
      use_exon_names <- gene_vals$use_exon_names;
      exon_range <- gene_vals$exon_range;
      gene_coords <- gene_vals$gene_coords;

      ## get gene coordinate range
      if (verbose) {
         jamba::printDebug("get_display_coords(): ",
            "gene: ", gene,
            ", use_exon_names: ", use_exon_names,
            ", exon_range: ", exon_range,
            ", gene_coords: ", jamba::formatInt(gene_coords),
            sep="-");
      }
      if ("exon names" %in% use_exon_names) {
         # Convert exon names to coordinates
         gene_coords <- range(as.data.frame(range(
            subset(flatExonsByGene1[[gene]], gene_nameExon %in% exon_range)
         ))[,c("start", "end")]);
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "input$exon_range: ", exon_range,
               ", inferred gene_coords: ",
               jamba::formatInt(gene_coords),
               sep="-");
         }
      } else {
         if (length(gene_coords) == 0) {
            # default has not been initialized yet, so take full gene range
            gene_coords <- range(as.data.frame(range(
               flatExonsByGene1[[gene]]
            ))[,c("start", "end")]);
         }
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "input$gene_coords:",
               jamba::formatInt(gene_coords),
               sep="-");
         }
      }
      gene_coords;
   })

   # the main function to prepare sashimi data for display
   # when gene is empty, "" or "blank" it returns NULL
   get_sashimi_data <- reactive({
      ## Comment out input$calc_gene_params since we depend
      ## upon its reactive effects
      #input$calc_gene_params;
      #shinyjs::disable("calc_gene_params");
      if (!exists("flatExonsByGene") ||
            !exists("filesDF")) {
         return(NULL);
      }

      if (!exists("use_memoise")) {
         use_memoise <- TRUE;
      }

      ## Wrap the workflow in a progress bar
      if (!exists("prepareSashimi_m")) {
         jamba::printDebug("get_sashimi_data(): ",
            "create memoise function prepareSashimi_m()");
         prepareSashimi_m <- memoise::memoise(prepareSashimi,
            cache=memoise::cache_filesystem("sashimidata_memoise"));
      }

      # move down below the preparatory steps above
      flat_list <- get_flat_gene_exons_plot();
      flatExonsByGene1 <- flat_list$flatExonsByGene;
      gene_vals <- flat_list$gene_vals;
      gene <- gene_vals$gene;
      use_exon_names <- gene_vals$use_exon_names;
      exon_range <- gene_vals$exon_range;
      gene_coords <- gene_vals$gene_coords;
      min_junction_reads <- gene_vals$min_junction_reads;
      include_strand <- gene_vals$include_strand;

      ## Define sample_id from sample selection
      sample_id_list <- get_sample_id_dt();
      sample_id <- sample_id_list$sample_id;
      layout_ncol <- sample_id_list$layout_ncol;

      if (verbose) {
         jamba::printDebug("sashimiAppServer(): ",
            "Using sample_id:",
            sample_id);
         jamba::printDebug("sashimiAppServer(): ",
            "Using gene:",
            gene);
      }
      if (length(sample_id) == 0) {
         return(NULL)
      }

      # gene is assigned only after reactive get_flat_gene_exons_plot()
      if (length(gene) == 0 || nchar(gene) == 0 || "blank" %in% gene) {
         return(NULL);
      }

      withProgress(
         message="Preparing Sashimi data.",
         value=0,
         {
            if (verbose && use_memoise) {
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
            if (dodebug) assign("sashimi_data", value=sashimi_data, envir=globalenv());
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
                        sep=" ",
                        c("Invalid memoise cache file, re-running prepareSashimi() without memoise.",
                           "\nPlease update memoise to version 1.1.0.900 or higher,\n",
                           "for proper cache cleaning methods. See Github:",
                           " 'r-lib/memoise'"),
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
      return(list(
         sashimi_data=sashimi_data,
         sample_id_list=sample_id_list,
         flat_list=flat_list));
   });

   # main Sashimi plot render function
   output$sashimiplot_output <- renderUI({
      sashimi_data_list <- get_sashimi_data();

      render_blank_plot <- function() {
         jamba::printDebug("sashimiAppServer(): ",
            "Rendering blank panel for ",
            "sashimiplot_output");
         return(tagList(renderPlot(
            height=300,
            ggplot2::ggplot() + ggplot2::theme_void()
         )));
      }
      if (length(sashimi_data_list) == 0 || length(sashimi_data_list$sashimi_data) == 0) {
         return(render_blank_plot())
      }

      # get associated values
      sashimi_data <- sashimi_data_list$sashimi_data;
      sample_id_list <- sashimi_data_list$sample_id_list;
      sample_id <- sample_id_list$sample_id;
      layout_ncol <- sample_id_list$layout_ncol;
      flat_list <- sashimi_data_list$flat_list;
      flatExonsByGene1 <- flat_list$flatExonsByGene;
      gene_vals <- flat_list$gene_vals;
      gene <- gene_vals$gene;
      use_exon_names <- gene_vals$use_exon_names;
      exon_range <- gene_vals$exon_range;
      gene_coords <- gene_vals$gene_coords;
      min_junction_reads <- gene_vals$min_junction_reads;
      include_strand <- gene_vals$include_strand;
      if (length(gene) == 0 || length(gene_coords) == 0 || length(sample_id) == 0) {
         return(render_blank_plot())
      }

      if ("ref2c" %in% names(attributes(sashimi_data))) {
         ref2c <- attr(sashimi_data, "ref2c");
      } else if ("ref2c" %in% names(sashimi_data)) {
         ref2c <- sashimi_data$ref2c;
      } else {
         ref2c <- NULL;
      }

      ## main sashimi plot method
      if (share_y_axis_d()) {
         facet_scales <- "fixed";
      } else {
         facet_scales <- "free_y";
      }

      ## call plotSashimi() with do_highlight=TRUE so the plot data
      ## will be prepared once and rendering can change independently
      #
      # do_plotly <- do_plotly_d()
      # do_highlight <- (enable_highlights_d() && do_plotly_d());
      # do_highlight <- (do_plotly_d());

      ## Obtain the baseline ggplot object
      display_coords <- get_display_coords();
      gg_sashimi <- plotSashimi(sashimi_data,
         show=c("coverage", "junction", "junctionLabels"),
         color_sub=color_sub,
         do_highlight=do_plotly_d(),
         facet_scales=facet_scales,
         junc_alpha=junction_alpha_d(),
         label_coords=display_coords,
         ref2c=ref2c,
         fill_scheme="sample_id");
      if (dodebug) assign("gg_sashimi", value=gg_sashimi, envir=globalenv());

      # Check layout_ncol
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

      ## Optionally prepare gene-exon model
      if (show_gene_model_d()) {
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "Getting gene model with get_display_coords(): ",
               jamba::formatInt(display_coords), sep="-");
         }
         if (show_tx_model_d() && length(flatExonsByTx) > 0) {
            if (show_detected_tx_d()) {
               if (verbose) {
                  jamba::printDebug("preparing to show detected transcripts:",
                     length(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]));
                  if (length(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]) == 0) {
                     jamba::printDebug("No transcripts in flatExonsByTx matched detectedTx. Showing head(names(flatExonsByTx)), then head(flatExonsByTx):");
                     print(head(names(flatExonsByTx)));
                     print(head(flatExonsByTx));
                  } else {
                     print(head(subset(flatExonsByTx[names(flatExonsByTx) %in% detectedTx]@unlistData, gene_name %in% gene)));
                  }
               }
               gg_gene <- gene2gg(gene=gene,
                  flatExonsByGene=flatExonsByGene1,
                  flatExonsByTx=flatExonsByTx[names(flatExonsByTx) %in% detectedTx],
                  label_coords=display_coords,
                  ref2c=ref2c,
                  exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
            } else {
               # TODO: recalculate flat exons for the gene using all transcripts
               gg_gene <- gene2gg(gene=gene,
                  flatExonsByGene=flatExonsByGene1,
                  flatExonsByTx=flatExonsByTx,
                  label_coords=display_coords,
                  ref2c=ref2c,
                  exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
            }
         } else {
            gg_gene <- gene2gg(gene=gene,
               flatExonsByGene=flatExonsByGene1,
               label_coords=display_coords,
               ref2c=ref2c,
               exonLabelSize=14 + exon_label_size_d() + font_sizing_d());
         }
         # why is this statement here?
         gg_gene <- gg_gene;
      }

      ## Prepare text label with genome coordinates
      ref_name <- as.character(head(seqnames(flatExonsByGene1[[gene]]), 1));
      coord_label <- paste0(
         ref_name,
         ":",
         paste(jamba::formatInt(gene_coords), collapse="-")
      );
      output$sashimitext_output <- renderText({
         HTML(paste0(" Region ",
            coord_label))
      });

      ## Prepare plotly or ggplot2 output
      num_samples <- length(unique(sample_id));
      panel_height <- panel_height_d();
      gene_panel_height <- gene_panel_height_d();
      # adjust base_size for fonts using exponent based upon panel height
      font_exp <- 1/3;
      base_font_size <- 16 + 1.0 * font_sizing_d();
      base_size <- base_font_size * (panel_height^font_exp)/(250^font_exp);
      #plot_height <- panel_height * (num_samples + show_gene_model_d());
      plot_heights <- c(sashimi=panel_height * num_samples,
         gene=gene_panel_height * show_gene_model_d());
      plot_height <- sum(plot_heights);
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
         ##########################################################
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "Preparing plotly output.");
         }
         display_coords <- get_display_coords();
         if (show_gene_model_d()) {
            ## use plotly, showing gene model
            ggly1 <- plotly::ggplotly(
               gg_sashimi +
                  colorjam::theme_jam(base_size=base_size) +
                  theme(axis.text.x=element_blank()) +
                  xlab(NULL) +
                  coord_cartesian(xlim=display_coords),
               tooltip="text",
               plot_height=plot_heights[1])
               # %>% plotly::style(hoveron="fill")
            ggly2 <- plotly::ggplotly(
               gg_gene +
                  colorjam::theme_jam(base_size=base_size) +
                  ggtitle(NULL) +
                  xlab(ref_name) +
                  coord_cartesian(xlim=display_coords),
               tooltip="text",
               plot_height=plot_heights[2]);
            p_heights <- unname(plot_heights/sum(plot_heights));
            gg_ly <- tryCatch({
               plotly::subplot(
                  ggly1,
                  ggly2,
                  nrows=2,
                  shareX=TRUE,
                  heights=p_heights);
            }, error=function(e){
               jamba::printDebug("Error:");
               print(e);
               ggly1
            })
               # %>% plotly::layout(height=sum(plot_heights))
         } else {
            ## use plotly, without showing gene model
            gg_ly <- plotly::ggplotly(
               gg_sashimi +
                  colorjam::theme_jam(base_size=base_size) +
                  xlab(ref_name) +
                  coord_cartesian(xlim=get_display_coords()),
               tooltip="text",
               height=plot_height
            )
            # %>% plotly::style(hoveron="fill");
         }
         if (enable_highlights_d()) {
            gg_ly <- plotly::highlight(gg_ly,
                  "plotly_hover",
                  opacityDim=0.8,
                  selected=plotly::attrs_selected(
                     line=list(color="#444444")));
         }
         ## Remove the color legend (again)
         if (!input$plotly_legend) {
            gg_ly <- plotly::layout(gg_ly, showlegend=FALSE);
         }
         # Try converting to plotlyOutput to define a fixed plot height
         output$plotly <- plotly::renderPlotly({
            plotly::layout(gg_ly,
               margin=list(r=100, l=70, t=20, b=70))
         });
         plotly_out <- plotly::plotlyOutput(
            "plotly",
            height=plot_height
         );
         return(plotly_out);
      } else {
         ##########################################################
         if (verbose) {
            jamba::printDebug("sashimiAppServer(): ",
               "Preparing ggplot output.");
         }
         ## Non-plotly static plot output
         display_coords <- get_display_coords();
         if (verbose) {
            jamba::printDebug("output$sashimiplot_output: ",
               "display_coords: ",
               jamba::formatInt(display_coords), sep="-");
         }
         if (length(display_coords) > 0) {
            #if (length(reactive_gene_coords$xlim) == 0) {
               reactive_gene_coords$xlim <- display_coords;
            #}
            #jamba::printDebug("reactive_gene_coords$xlim: ", formatInt(reactive_gene_coords$xlim), sep="-");
            if (show_gene_model_d()) {
               ## With gene-transcript-exon model
               cp <- cowplot::plot_grid(
                  gg_sashimi +
                     colorjam::theme_jam(base_size=base_size) +
                     theme(axis.text.x=element_blank()) +
                     xlab(NULL) +
                     coord_cartesian(xlim=reactive_gene_coords$xlim),
                  gg_gene +
                     colorjam::theme_jam(base_size=base_size) +
                     ggtitle(NULL) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=reactive_gene_coords$xlim),
                  ncol=1,
                  align="v",
                  axis="lr",
                  rel_heights=plot_heights
               )
               output[["sashimi_plot"]] <- renderPlot(
                  height=plot_height,
                  cp
               )
               tagList(shiny::plotOutput("sashimi_plot",
                  height=paste0(plot_height, "px")))
            } else {
               ## No gene-transcript-exon model
               # slightly different method for dynamic zoom inside ggplot2
               output[["sashimi_plot"]] <- renderPlot(
                  height=plot_height,
                  #suppressMessages(
                     gg_sashimi +
                        colorjam::theme_jam(base_size=base_size) +
                        xlab(ref_name) +
                        coord_cartesian(xlim=reactive_gene_coords$xlim)
                  #)
               )
               observeEvent(input[["sashimi_plot_dblclick"]], {
                  brush <- input[["sashimi_plot_brush"]]
                  if (!is.null(brush)) {
                     xrange <- c(brush$xmin, brush$xmax);
                     coordrange <- ref2c$inverse(xrange);
                     reactive_gene_coords$xlim <- c(coordrange);
                  } else {
                     reactive_gene_coords$xlim <- gene_coords;
                  }
               })
               tagList(shiny::plotOutput("sashimi_plot",
                  height=paste0(plot_height, "px"),
                  dblclick="sashimi_plot_dblclick",
                  brush=brushOpts(
                     direction="x",
                     id="sashimi_plot_brush",
                     resetOnNew=TRUE)
               ))
            }
         } else {
            ## Fallback when get_display_coords() is empty
            ## sometimes happens at startup, unclear when/why
            tagList(renderPlot(
               height=plot_height,
               gg_sashimi +
                  colorjam::theme_jam(base_size=base_size)
            ));
         }
      }
   });

   reactive_gene_coords <- reactiveValues(xmin=NULL, xmax=NULL)

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
            dt <- DT::formatStyle(dt,
               iCol,
               backgroundColor=DT::styleEqual(
                  levels=names(color_sub),
                  values=jamba::rgb2col(col2rgb(color_sub))
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

      files_dt <- dplyr::select(files_dt,
         sample_id,
         everything(),
         -scale_factor, -type, -url,
         type, scale_factor, url);
      files_dt <- dplyr::arrange(files_dt,
         sample_id, type)
      files_dt <- DT::datatable(files_dt,
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
         ))
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
      layout_ncol <- isolate(input$layout_ncol);
      if (verbose >= 0) {
         jamba::printDebug("get_sample_id_dt(): ",
            new_sample_id);
      }
      return(list(sample_id=new_sample_id,
         layout_ncol=layout_ncol))
   })

   # sample table where selected rows are re-ordered on click
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
         selected=match(get_sample_id_dt()$sample_id,
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
