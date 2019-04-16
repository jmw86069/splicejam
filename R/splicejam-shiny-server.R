
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
            print(sdim(sashimi_data));
            sashimi_data;
         }
      );
   });

   output$sashimiplot_output <- renderPlot({
      sashimi_data <- get_sashimi_data();
      if (length(sashimi_data) == 0) {
         ggplot() + geom_blank();
      } else {
         gg_sashimi <- plotSashimi(sashimi_data,
            junc_color=alpha2col("goldenrod1", 0.1),
            fill_scheme="sample_id");
         print(gg_sashimi);
      }
   });

}
