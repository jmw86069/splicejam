
#' Sashimi Shiny app UI
#'
#' Sashimi Shiny app UI
#'
#' @family Sashimi Shiny functions
#'
#' @import shiny
#' @import shinydashboard
#' @import shinydashboardPlus
#' @import shinyWidgets
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
   options("warn"=-1);

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
               fluidRow(
                  column(
                     width=4,
                     style="padding:0px",
                     shinyWidgets::prettySwitch(
                        inputId="use_exon_names",
                        value=FALSE,
                        slim=TRUE,
                        status="info",
                        width="90%",
                        label="Use exon names in the slider below?")
                  ),
                  column(
                     width=4,
                     style="padding:0px",
                     numericInput(
                        inputId="min_junction_reads",
                        label="Minimum junction reads",
                        value=100,
                        width="90%",
                        step=1,
                        min=1,
                        max=1000
                     )
                  ),
                  column(
                     width=4,
                     style="padding:0px",
                     selectInput(
                        inputId="include_strand",
                        label="Show coverage by strand:",
                        choices=c("both", "+", "-"),
                        selected="both",
                        multiple=FALSE,
                        width="90%"
                     )
                  )
               ),
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
                  htmltools::em(textOutput("gene_coords_label",
                     inline=TRUE)),
                  sliderInput(
                     "gene_coords",
                     label=NULL,#"Genome coordinate range",
                     min=1,
                     max=2000,
                     width="80%",
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
               solidheader=TRUE,
               closable=FALSE,
               collapsible=TRUE,
               enable_sidebar=TRUE,
               width=12,
               height="100%",
               fluidRow(
                  column(
                     width=12,
                     textOutput("sashimitext_output"),
                     uiOutput("sashimiplot_output")
                  )
               ),
               sidebar_width=25,
               sidebar_start_open=FALSE,
               sidebar_content=tagList(
                  shinyWidgets::prettySwitch(
                     inputId="do_plotly",
                     value=TRUE,
                     slim=TRUE,
                     status="info",
                     label="Interactive?"),
                  shinyWidgets::prettySwitch(
                     inputId="show_gene_model",
                     value=TRUE,
                     slim=TRUE,
                     status="info",
                     label="Show gene-exon model?"
                  ),
                  conditionalPanel(
                     condition="input.show_gene_model == true",
                     shinyWidgets::prettySwitch(
                        inputId="show_tx_model",
                        value=FALSE,
                        slim=TRUE,
                        status="info",
                        label="Include transcript isoforms?"
                     ),
                     conditionalPanel(
                        condition="input.show_tx_model == true",
                        shinyWidgets::prettySwitch(
                           inputId="show_detected_tx",
                           value=TRUE,
                           slim=TRUE,
                           status="info",
                           label="Show detected transcripts only?"
                        )
                     )
                  ),
                  shinyWidgets::prettySwitch(
                     inputId="share_y_axis",
                     value=TRUE,
                     slim=TRUE,
                     status="info",
                     label="Shared y-axis range?"
                  ),
                  sliderInput(
                     inputId="exon_label_size",
                     label="Exon label size",
                     min=2,
                     max=20,
                     value=6,
                     step=0.5
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
      setShadow(class="box"),
      setShadow(class="boxPlus"),
      tabItems(
         tabItem(tabName="guides", guidesTab),
         tabItem(tabName="sashimiplot", sashimiplotTab)
      )
   );

   dp <- dashboardPage(
      header,
      sidebar,
      body,
      skin="blue");
   printDebug("class(dp):", class(dp));
   dp;

}
