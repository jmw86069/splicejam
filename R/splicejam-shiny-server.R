
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


   observe({
      gene <- input$gene;
      printDebug("updateInputSlider gene:", gene);
      if (length(gene) > 0 && nchar(gene) > 0) {
         chr_range <- as.data.frame(range(flatExonsByGene[[gene]]))[,c("start", "end")];
         if (length(chr_range) > 0) {
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
      withProgress(
         message="Preparing Sashimi data.",
         value=0,
         {
            sashimi_data <- prepareSashimi_m(gene=gene,
               flatExonsByGene=flatExonsByGene,
               minJunctionScore=100,
               sample_id=c("CA1_CB", "CA2_CB"),
               filesDF=filesDF,
               verbose=TRUE,
               do_shiny_progress=FALSE);
            sashimi_data;
         }
      );
   });

#   output$sashimiplot_output <- renderPlot({
   output$sashimiplot_output <- renderUI({
      sashimi_data <- get_sashimi_data();
      if (length(sashimi_data) == 0) {
         tagList(renderPlot(
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
            facet_scales=facet_scales,
            fill_scheme="sample_id");
         ## Optionally prepare gene-exon model
         if (input$show_gene_model) {
            gg_gene <- gene2gg(gene=gene,
               flatExonsByGene=flatExonsByGene,
               exonLabelSize=input$exon_label_size);
         }

         ## Optionally get gene coordinate range
         if (isolate(input$use_exon_names)) {
            exon_range <- isolate(input$exon_range);
            # Convert exon names to coordinates
            gene_coords <- range(as.data.frame(range(
               subset(flatExonsByGene[[gene]], gene_nameExon %in% exon_range)
            ))[,c("start", "end")]);
         } else {
            gene_coords <- isolate(input$gene_coords);
         }

         ## Prepare plotly or ggplot2 output
         ref_name <- as.character(head(seqnames(flatExonsByGene[[gene]]), 1));
         sample_id <- unique(filesDF$sample_id);
         num_samples <- max(c(length(sample_id), 1)) + 1;
         num_samples <- 2;
         plot_height <- 250 * (num_samples +
               input$show_gene_model);
         printDebug("num_samples:", num_samples,
            ", plot_height:", plot_height);
         if (input$do_plotly) {
            if (input$show_gene_model) {
               gg_ly <- suppressMessages(
                  plotly::subplot(
                     plotly::ggplotly(
                        gg_sashimi +
                           theme(axis.text.x=element_blank()) +
                           xlab(NULL),
                        tooltip="text"),
                     plotly::ggplotly(
                        gg_gene +
                           ggtitle(NULL) +
                           xlab(ref_name),
                        tooltip="text"),
                     nrows=2,
                     heights=c(num_samples, 1)/(num_samples + 1),
                     shareX=TRUE
                  ) %>% layout(height=plot_height)
               );
            } else {
               gg_ly <- suppressMessages(
                  plotly::ggplotly(
                     gg_sashimi +
                        scale_y_continuous(labels=scales::comma) +
                        xlab(ref_name) +
                        coord_cartesian(xlim=gene_coords),
                     tooltip="text",
                     height=plot_height
                  )
               );
            }
            ## Remove the color legend (again)
            gg_ly <- gg_ly %>%
               layout(showlegend=FALSE);
            #%>%plotly::layout(autosize=TRUE)
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
                              scale_y_continuous(labels=scales::comma) +
                              theme(axis.text.x=element_blank()) +
                              xlab(NULL) +
                              coord_cartesian(xlim=gene_coords),
                           gg_gene +
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
                           xlab(ref_name) +
                           coord_cartesian(xlim=gene_coords)
                     )
                  ));
               }
            } else {
               tagList(renderPlot(
                  height=plot_height,
                  gg_sashimi
               ));
            }
         }
      }
   });

}
