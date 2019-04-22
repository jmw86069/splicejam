
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
#'
#' @export
sashimiAppServer <- function
(input,
 output,
 session)
{
   #
   printDebug("sashimiAppServer(): ",
      "class(sashimiplot_guide):",
      class(sashimiplot_guide));
   printDebug("sashimiAppServer(): ",
      "exists('flatExonsByGene'):",
      exists('flatExonsByGene'));
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
      gene <- isolate(input$gene);
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
               do_shiny_progress=TRUE);
            sashimi_data;
         }
      );
   });

#   output$sashimiplot_output <- renderPlot({
   output$sashimiplot_output <- renderUI({
      sashimi_data <- get_sashimi_data();
      if (length(sashimi_data) == 0) {
         tagList(renderPlot(
            ggplot() + geom_blank()
            )
         );
      } else {
         gg_sashimi <- plotSashimi(sashimi_data,
            junc_color=alpha2col("goldenrod1", 0.1),
            fill_scheme="sample_id");
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
         if (input$do_plotly) {
            gg_ly <- plotly::ggplotly(
               gg_sashimi +
                  scale_y_continuous(labels=scales::comma) +
                  coord_cartesian(xlim=gene_coords)
            );
            if (input$do_rangeslider) {
               gg_ly <- gg_ly %>%
                  plotly::rangeslider();
            }
            tagList(gg_ly);
         } else {
            if (length(gene_coords) > 0) {
               tagList(renderPlot(gg_sashimi +
                     coord_cartesian(xlim=gene_coords)
               ));
            } else {
               tagList(renderPlot(gg_sashimi));
            }
         }
      }
   });

}
