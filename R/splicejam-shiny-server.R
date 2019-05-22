
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

   ## server-side selectize gene list
   printDebug("length(detectedGenes):", length(detectedGenes));
   updateSelectizeInput(session,
      "gene",
      choices=detectedGenes,
      selected="Gria1",
      server=TRUE);

   output$gene_coords_label <- renderText({"Genome coordinate range"});

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
      gene <<- isolate(input$gene);
      if (!exists("flatExonsByGene") ||
            !exists("filesDF")) {
         return(NULL);
      }

      ## Wrap the workflow in a progress bar
      prepareSashimi_m <- memoise::memoise(prepareSashimi,
         cache=memoise::cache_filesystem("sashimidata_memoise"));
      ## Define sample_id for now
      if (!exists("sample_id")) {
         sample_id <<- head(unique(filesDF$sample_id), 3);
      }
      min_junction_reads <- isolate(input$min_junction_reads);
      include_strand <- isolate(input$include_strand);
      if (!exists("verbose")) {
         verbose <- FALSE;
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
               do_shiny_progress=TRUE);
            sashimi_data;
         }
      );
   });

#   output$sashimiplot_output <- renderPlot({
   output$sashimiplot_output <- renderUI({
      sashimi_data <- get_sashimi_data();
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
         if (input$share_y_axis) {
            facet_scales <- "fixed";
         } else {
            facet_scales <- "free_y";
         }
         gg_sashimi <- plotSashimi(sashimi_data,
            color_sub=color_sub,
            do_highlight=TRUE,
            facet_scales=facet_scales,
            fill_scheme="sample_id");
         sashimi_data <- sashimi_data;
         gg_sashimi <<- gg_sashimi;
         ## Optionally prepare gene-exon model
         if (input$show_gene_model) {
            if (input$show_tx_model && length(flatExonsByTx) > 0) {
               if (input$show_detected_tx) {
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene,
                     flatExonsByTx=flatExonsByTx[names(flatExonsByTx) %in% detectedTx],
                     exonLabelSize=input$exon_label_size);
               } else {
                  gg_gene <- gene2gg(gene=gene,
                     flatExonsByGene=flatExonsByGene,
                     flatExonsByTx=flatExonsByTx,
                     exonLabelSize=input$exon_label_size);
               }
            } else {
               gg_gene <- gene2gg(gene=gene,
                  flatExonsByGene=flatExonsByGene,
                  exonLabelSize=input$exon_label_size);
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
            paste0("  ",
               coord_label)
         });

         ## Prepare plotly or ggplot2 output
         #sample_id <- unique(filesDF$sample_id);
         #num_samples <- max(c(length(sample_id), 1)) + 1;
         #num_samples <- 2;
         num_samples <- length(unique(sample_id));
         panel_height <- input$panel_height;
         # adjust base_size for fonts using exponent based upon panel height
         font_exp <- 1/3;
         base_size <- 12 * (panel_height^font_exp)/(250^font_exp);
         plot_height <- panel_height * (num_samples + input$show_gene_model);
         printDebug("num_samples:", num_samples,
            ", panel_height:", panel_height,
            ", base_size:", format(digits=1, base_size),
            ", plot_height:", plot_height);
         if (input$do_plotly) {
            if (input$show_gene_model) {
               ## use plotly, showing gene model
               ggly1 <- plotly::ggplotly(
                  gg_sashimi +
                     theme_jam(base_size=base_size) +
                     theme(axis.text.x=element_blank()) +
                     xlab(NULL) +
                     coord_cartesian(xlim=gene_coords),
                  tooltip="name") %>%
                  plotly::style(
                     hoveron="fill"
                  );
               if (input$enable_highlights) {
                  ggly1 <- ggly1 %>%
                     plotly::highlight(
                        on="plotly_hover",
                        off="plotly_doubleclick",
                        opacityDim=0.8,
                        selected=attrs_selected(
                           line=list(color="#444444")));
               }
               ggly2 <- plotly::ggplotly(
                  gg_gene +
                     theme_jam(base_size=base_size) +
                     ggtitle(NULL) +
                     xlab(ref_name) +
                     coord_cartesian(xlim=gene_coords),
                  tooltip="text");
               gg_ly <- suppressMessages(
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
               if (input$enable_highlights) {
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
               if (input$show_gene_model) {
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
                           rel_heights=c(num_samples, num_samples+1)
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

   # filesDF table output
   output$files_df <- DT::renderDT({
      options("rowGroup.dataSrc"="sample_id");

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
      files_dt$selected <- shinyCheckboxInput(
         checkboxInput,
         nrow(files_dt),
         "cbox_");

      files_dt <- files_dt %>%
         dplyr::select(sample_id, selected, type, everything(), -url, url) %>%
         dplyr::arrange(sample_id, type) %>%
         DT::datatable(
            editable=TRUE,
            rownames=FALSE,
            escape=FALSE,
            height='15px',
            extensions=c("RowReorder"),
            options=list(
               autoWidth=TRUE,
               deferRender=TRUE,
               columnDefs=list(list(width='50px', targets=list(1,2,3))),
               pageLength=24,
               lengthMenu=c(12,24,48,120),
               rowReorder=TRUE,
               preDrawCallback=DT::JS('function() {
Shiny.unbindAll(this.api().table().node()); } '),
               drawCallback=DT::JS('function() {
Shiny.bindAll(this.api().table().node()); } ')
            )
         );
         # optionally colorize sample_id using color_sub
         color_sub["junction"] <- "slateblue1";
         color_sub["bw"] <- "darkslategray1";

         # check each colname for matching colors
         for (iCol in colnames(filesDF)) {
            if (all(unique(filesDF[[iCol]]) %in% names(color_sub))) {
               files_dt <- files_dt %>%
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
         files_dt;
      },
      options=list(
         lengthChange=FALSE,
         pageLength=24
      )
   );

}
