
#' Sashimi Shiny app UI
#'
#' Sashimi Shiny app UI
#'
#' @family Sashimi Shiny functions
#'
#' @import shiny
#' @import shinydashboard
#'
#' @param ... additional arguments are ignored.
#'
#' @export
sashimiAppUI <- function
(...)
{
   # header
   header <- dashboardHeader(
      title=tagList("Splicejam Sashimi Viewer",
         icon("map"))
   );

   # sidebar
   sidebar <- dashboardSidebar(
      sidebarMenu(
         id="tabs",
         menuItem(
            text="Guides",
            tabName="guides",
            icon=icon("info")),
         menuItem(
            text="Sashimi Plot",
            tabName="sashimiplot",
            icon=icon("map"))
      )
   );

   # define Sashimi Plot tab
   sashimiplotTab <- fluidPage(
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            shinydashboard::box(
               title="Sashimi Plot Parameters",
               status="warning",
               solidHeader=TRUE,
               width=12,
               selectizeInput(
                  label="Select Gene",
                  inputId="gene",
                  selected="Gria1",
                  choices=c("Gria1", "Ntrk2"),
                  options=list(maxOptions=100),
                  multiple=FALSE
               ),
               actionButton("calc_gene_params",
                  label="Update Sashimi Plots")
            )
         )
      ),
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            shinydashboard::box(
               title="Sashimi Plot",
               status="primary",
               solidheader=FALSE,
               width=12,
               plotOutput("sashimiplot_output")
            )
         )
      )
   );

   ## Load guidesTab via salsaAppConstants()
   #sashimiAppConstants();

   # dashboard body
   body <- dashboardBody(
      tabItems(
         tabItem(tabName="guides", guidesTab),
         tabItem(tabName="sashimiplot", sashimiplotTab)
      )
   );

   printDebug("class(header):", class(header));
   printDebug("class(sidebar):", class(sidebar));
   printDebug("class(body):", class(body));
   dp <- dashboardPage(
      header,
      sidebar,
      body,
      skin="blue");
   printDebug("class(dp):", class(dp));
   dp;

}
