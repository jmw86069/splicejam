
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
   header <- dashboardHeaderPlus(
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
               collapsible=TRUE,
               width=12,
               selectizeInput(
                  label="Select Gene",
                  inputId="gene",
                  selected="Gria1",
                  choices=c("Gria1", "Ntrk2"),
                  options=list(maxOptions=100),
                  multiple=FALSE
               ),
               shinyWidgets::materialSwitch(
                  inputId="use_exon_names",
                  value=FALSE,
                  status="info",
                  label="Use Exon Names?"),
               conditionalPanel(
                  condition="input.use_exon_names == true",
                  sliderTextInput(
                     inputId="exon_range",
                     label="Gene exon range",
                     grid=TRUE,
                     force_edges=TRUE,
                     choices=c("exon1", "exon2", "exon3"),
                     selected=c("exon1", "exon3")
                  )
               ),
               conditionalPanel(
                  condition="input.use_exon_names == false",
                  sliderInput(
                     "gene_coords",
                     label="Genome coordinate range",
                     min=1,
                     max=2000,
                     value=c(2, 2000),
                     step=1,
                     round=TRUE
                  )
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
            shinydashboardPlus::boxPlus(
               title="Sashimi Plot",
               status="primary",
               solidheader=FALSE,
               closable=FALSE,
               enable_sidebar=TRUE,
               width=12,
               #plotOutput("sashimiplot_output"),
               uiOutput("sashimiplot_output"),
               sidebar_width=35,
               sidebar_start_open=FALSE,
               sidebar_content=tagList(
                  shinyWidgets::prettySwitch(
                     inputId="do_plotly",
                     value=FALSE,
                     slim=TRUE,
                     fill=TRUE,
                     status="info",
                     label="Interactive"),
                  conditionalPanel(
                     condition="input.do_plotly == true",
                     shinyWidgets::prettySwitch(
                        inputId="do_rangeslider",
                        value=FALSE,
                        slim=TRUE,
                        fill=TRUE,
                        status="info",
                        label="Dynamic range slider")
                  )
               )
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
