
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
#'
#' @export
sashimiAppServer <- function
(input,
 output,
 session)
{
   #
   printDebug("class(sashimiplot_guide):",
      class(sashimiplot_guide));
   output$sashimiplot_guide <- renderUI(sashimiplot_guide);

   sashimi_data <- reactive({
      gene <- input$gene;
      sashimi_data <- prepareSashimi(gene=gene,
         flatExonsByGene=flatExonsByGeneMm10,
         minJunctionScore=100,
         sample_id=c("CA1_CB", "CA2_CB"),
         filesDF=filesDF)
   });

   output$sashimiplot_output <- renderPlot({
      # button "calc_gene_params"
      # Require non-empty numi_per_cell data.frame
      #req(numi_per_cell());
      #numi_per_cell_df <- numi_per_cell();
      #gg_rank_count_cell <- rank_count_plot(numi_per_cell_df[,2]);
      #gg_rank_count_cell;
      ggplot() + geom_blank();
   });

}
