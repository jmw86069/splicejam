
#' Sashimi Shiny app server
#'
#' Sashimi Shiny app sserver
#'
#' @param input provided by shiny
#' @param output provided by shiny
#' @param stssion provided by shiny
#'
#' @family splicejam R-shiny functions
#'
#' @import jamba
#' @import dplyr
#' @import ggplot2
#' @import plotly
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

   ## server-side selectize gene list
   printDebug("length(detectedGenes):", length(detectedGenes));
   updateSelectizeInput(session,
      "gene",
      choices=detectedGenes,
      selected="Gria1",
      server=TRUE);

   output$gene_coords_label <- renderText({"Genome coordinate range"});


   # use debounce to slow down rendering a plot while changing parameters
   debounce_ms <- 1000;
   share_y_axis <- reactive({
      input$share_y_axis;
   });
   share_y_axis_d <- debounce(share_y_axis, debounce_ms);
   show_gene_model <- reactive({
      input$show_gene_model;
   });
   show_gene_model_d <- debounce(show_gene_model, debounce_ms);
   show_tx_model <- reactive({
      input$show_tx_model;
   });
   show_tx_model_d <- debounce(show_tx_model, debounce_ms);
   show_detected_tx <- reactive({
      input$show_detected_tx;
   });
   show_detected_tx_d <- debounce(show_detected_tx, debounce_ms);
   exon_label_size <- reactive({
      input$exon_label_size;
   });
   exon_label_size_d <- debounce(exon_label_size, debounce_ms);
   panel_height <- reactive({
      input$panel_height;
   });
   panel_height_d <- debounce(panel_height, debounce_ms);
   font_sizing <- reactive({
      input$font_sizing;
   });
   font_sizing_d <- debounce(font_sizing, debounce_ms);
   do_plotly <- reactive({
      input$do_plotly;
   });
   do_plotly_d <- debounce(do_plotly, debounce_ms);
   enable_highlights <- reactive({
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
      printDebug("updateInputSlider gene:", gene);
      if (length(gene) > 0 && nchar(gene) > 0) {
         chr_range <- as.data.frame(range(flatExonsByGene[[gene]]))[,c("start", "end")];
         if (length(chr_range) > 0) {
            ## Update slider text label
            output$gene_coords_label <- renderText({
               paste0("Coordinate range ",
                  as.character(seqnames(head(flatExonsByGene[[gene]], 1))),
                  ":")
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
         if ("gene_nameExon" %in% colnames(values(flatExonsByGene[[gene]]))) {
            exon_names <- jamba::mixedSort(
               values(flatExonsByGene[[gene]])$gene_nameExon);
            if (length(exon_names) > 0) {
               updateSliderTextInput(session,
                  "exon_range",
                  choices=exon_names,
                  selected=c(head(exon_names, 1), tail(exon_names, 1))
               );
            }
         }
      }
   })

   get_sashimi_data <- reactive({
      input$calc_gene_params;
      shinyjs::disable("calc_gene_params");
      gene <<- isolate(input$gene);
      if (!exists("flatExonsByGene") ||
            !exists("filesDF")) {
         return(NULL);
      }

      ## Wrap the workflow in a progress bar
      prepareSashimi_m <- memoise::memoise(prepareSashimi,
         cache=memoise::cache_filesystem("sashimidata_memoise"));

      ## Define sample_id from sample selection
      sample_order <- isolate(input$selectionto_order);
      if (!exists("sample_id")) {
         if (length(sample_order) == 0) {
            sample_id <<- head(unique(filesDF$sample_id), 3);
         } else {
            sample_id <<- sample_order;
         }
      } else if (length(sample_order) > 0) {
         sample_id <<- sample_order;
      }
      printDebug("Using sample_id:", sample_id);
      min_junction_reads <- isolate(input$min_junction_reads);
      include_strand <- isolate(input$include_strand);
      if (!exists("verbose")) {
         verbose <- FALSE;
      }
      if (!exists("use_memoise")) {
         use_memoise <- TRUE;
      }
      withProgress(
         message="Preparing Sashimi data.",
         value=0,
         {
            sashimi_data <- prepareSashimi_m(gene=gene,
               flatExonsByGene=flatExonsByGene,
               minJunctionScore=min_junction_reads,
               sample_id=sample_id,
               filesDF=filesDF,
               include_strand=include_strand,
               verbose=verbose,
               use_memoise=use_memoise,
               do_shiny_progress=TRUE);
            sashimi_data;
         }
      );
   });

   output$sashimiplot_output <- renderUI({
      sashimi_data <- get_sashimi_data();
      # Update in parent environment
      sashimi_data <<- sashimi_data;
      if (length(sashimi_data) == 0) {
         ## Todo: make empty plot a minimum height in pixels
         tagList(
            renderPlot(
               height=300,
               ggplot() + geom_blank()
            )
         );
      } else {
         if (!exists("gene") || length(gene) == 0) {
            gene <<- isolate(input$gene);
         }
         if (share_y_axis_d()) {
            facet_scales <- "fixed";
         } else {
            facet_scales <- "free_y";
         }
         if (!isolate(input$do_plotly)) {
            do_highlight <- FALSE;
         } else {
            do_highlight <- TRUE;
         }
         gg_sashimi <- plotSashimi(sashimi_data,
            color_sub=color_sub,
            do_highlight=do_highlight,
            facet_scales=facet_scales,
            fill_scheme="sample_id");
         #sashimi_data <- sashimi_data;

         # Check layout_ncol
         layout_ncol <- isolate(input$layout_ncol);
         if (length(layout_ncol) > 0 && layout_ncol > 1) {
            # change from facet_grid to facet_wrap()
            printDebug("Applying facet_wrap() with ncol:",
               layout_ncol);
            gg_sashimi <- gg_sashimi +
               facet_wrap(~sample_id,
                  ncol=layout_ncol,
                  scales=facet_scales);
         }
         gg_sashimi <<- gg_sashimi;
         ## Optionally prepare gene-exon model
         if (show_gene_model_d()) {
            if (show_tx_model_d() && length(flatExonsByTx) > 0) {
               if (show_detected_tx_d()) {
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene,
                     flatExonsByTx=flatExonsByTx[names(flatExonsByTx) %in% detectedTx],
                     exonLabelSize=exon_label_size_d());
               } else {
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene,
                     flatExonsByTx=flatExonsByTx,
                     exonLabelSize=exon_label_size_d());
               }
            } else {
               gg_gene <- gene2gg(gene=gene,
                  flatExonsByGene=flatExonsByGene,
                  exonLabelSize=exon_label_size_d());
            }
            gg_gene <<- gg_gene;
         }

         ## Optionally get gene coordinate range
         if (isolate(input$use_exon_names) %in% "exon names") {
            exon_range <- isolate(input$exon_range);
            # Convert exon names to coordinates
            gene_coords <- range(as.data.frame(range(
               subset(flatExonsByGene[[gene]], gene_nameExon %in% exon_range)
            ))[,c("start", "end")]);
         } else {
            gene_coords <- isolate(input$gene_coords);
         }

         ## Prepare text label with genome coordinates
         ref_name <- as.character(head(seqnames(flatExonsByGene[[gene]]), 1));
         coord_label <- paste0(
            ref_name,
            ":",
            paste(gene_coords, collapse="-")
         );
         output$sashimitext_output <- renderText({
            HTML(paste0("Region ",
               coord_label))
         });

         ## Prepare plotly or ggplot2 output
         #sample_id <- unique(filesDF$sample_id);
         #num_samples <- max(c(length(sample_id), 1)) + 1;
         #num_samples <- 2;
         num_samples <- length(unique(sample_id));
         panel_height <- panel_height_d();
         # adjust base_size for fonts using exponent based upon panel height
         font_exp <- 1/3;
         font_sizing <- font_sizing_d();
         base_font_size <- dplyr::case_when(
            igrepHas("-2", font_sizing) ~ 8,
            igrepHas("-1", font_sizing) ~ 8,
            igrepHas("Default", font_sizing) ~ 12,
            igrepHas("+1", font_sizing) ~ 14,
            igrepHas("+2", font_sizing) ~ 16,
            igrepHas(".", font_sizing) ~ 12
         )
         base_size <- base_font_size * (panel_height^font_exp)/(250^font_exp);
         plot_height <- panel_height * (num_samples + show_gene_model_d());
         printDebug("num_samples:", num_samples,
            ", panel_height:", panel_height,
            ", base_size:", format(digits=1, base_size),
            ", base_font_size:", base_font_size,
            ", font_sizing:", font_sizing,
            ", layout_ncol:", layout_ncol,
            ", plot_height:", plot_height);
         if (do_plotly_d()) {
            if (show_gene_model_d()) {
               ## use plotly, showing gene model
               ggly1 <<- plotly::ggplotly(
                  gg_sashimi +
                     theme_jam(base_size=base_size) +
                     theme(axis.text.x=element_blank()) +
                     xlab(NULL) +
                     coord_cartesian(xlim=gene_coords),
                  tooltip="name") %>%
                  plotly::style(
                     hoveron="fill"
                  );
               if (enable_highlights_d()) {
                  ggly1 <- ggly1 %>%
                     plotly::highlight(
                        on="plotly_hover",
                        off="plotly_doubleclick",
                        opacityDim=0.8,
                        selected=attrs_selected(
                           line=list(color="#444444")));
               }
               ggly2 <<- plotly::ggplotly(
                  gg_gene +
                     theme_jam(base_size=base_size) +
                     ggtitle(NULL) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=gene_coords),
                  tooltip="text");
               gg_ly <<- suppressMessages(
                  plotly::subplot(
                     ggly1,
                     ggly2,
                     nrows=2,
                     heights=c(num_samples, 1)/(num_samples + 1),
                     shareX=TRUE
                  ) %>% layout(height=plot_height)
               );
               gg_ly <- gg_ly %>%
                  plotly::highlight("plotly_hover",
                     opacityDim=0.8,
                     selected=attrs_selected(
                        line=list(color="#444444")));
            } else {
               ## use plotly, showing gene model
               gg_ly <- plotly::ggplotly(
                  gg_sashimi +
                     theme_jam(base_size=base_size) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=gene_coords),
                  tooltip="name",
                  height=plot_height
               ) %>% style(hoveron="fill");
               if (enable_highlights_d()) {
                  gg_ly <- gg_ly %>%
                     highlight("plotly_hover",
                        opacityDim=0.8,
                        selected=attrs_selected(line=list(color="#444444")));
               }
            }
            ## Remove the color legend (again)
            gg_ly <- gg_ly %>%
               layout(showlegend=FALSE);
            # Try converting to plotlyOutput to define a fixed plot height
            output$plotly <- renderPlotly({
               gg_ly %>%
                  layout(margin=list(r=100, l=70, t=20, b=70))
            });
            plotlyOutput(
               "plotly",
               height=plot_height
            );
            #tagList(
            #   htmltools::as.tags(gg_ly)
            #);
         } else {
            ## Non-plotly static plot output
            if (length(gene_coords) > 0) {
               if (show_gene_model_d()) {
                  tagList(renderPlot(
                     height=plot_height,
                     suppressMessages(
                        cowplot::plot_grid(
                           gg_sashimi +
                              theme_jam(base_size=base_size) +
                              theme(axis.text.x=element_blank()) +
                              xlab(NULL) +
                              coord_cartesian(xlim=gene_coords),
                           gg_gene +
                              theme_jam(base_size=base_size) +
                              ggtitle(NULL) +
                              xlab(ref_name) +
                              coord_cartesian(xlim=gene_coords),
                           ncol=1,
                           align="v",
                           axis="lr",
                           rel_heights=c(num_samples, 2+show_tx_model_d()*1)
                        )
                     )
                  ));
               } else {
                  tagList(renderPlot(
                     height=plot_height,
                     suppressMessages(
                        gg_sashimi +
                           theme_jam(base_size=base_size) +
                           xlab(ref_name) +
                           coord_cartesian(xlim=gene_coords)
                     )
                  ));
               }
            } else {
               tagList(renderPlot(
                  height=plot_height,
                  gg_sashimi +
                     theme_jam(base_size=base_size)
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

   # Render the orderInput elements
   output$selection_sort <- renderUI({
      if (exists("filesDF") && length(filesDF$sample_id) > 0) {
         sample_items <- unique(filesDF$sample_id);
      } else {
         sample_items <- NULL;
      }
      if (exists("sample_id") && length(sample_id) > 0) {
         display_items <- unique(sample_id);
         sample_items <- setdiff(sample_items, display_items);
      } else {
         display_items <- NULL;
      }
      tagList(
         shinyjqui::orderInput(
            inputId="selectionfrom",
            label="Hidden samples:",
            items=sample_items,
            item_class="primary",
            connect="selectionto",
            width="100%",
            placeholder="(Drag samples here to hide)"
         ),
         shinyjqui::orderInput(
            inputId="selectionto",
            label="Display samples:",
            items=display_items,
            item_class="primary",
            connect="selectionfrom",
            width="100%",
            placeholder="(Drag samples here to display)"
         )
      )
   });

}
