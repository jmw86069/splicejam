
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
   nbsp <- HTML("&nbsp;");
   nbsp3 <- HTML("&nbsp;&nbsp;&nbsp;");

   # sidebar
   sidebar <- dashboardSidebar(
      sidebarMenu(
         id="tabs",
         menuItem(
            text="Sashimi Plot",
            tabName="sashimiplot",
            icon=icon("map")),
         menuItem(
            text="Sample Selection",
            tabName="sampleselect",
            icon=icon("clipboard-list")),
         menuItem(
            text="Guides",
            tabName="guides",
            icon=icon("info"))
         #menuItem(
         #   text="Samples and Data",
         #   tabName="samplesdata",
         #   icon=icon("table"))
      )
   );

   # define Sashimi Plot tab
   sashimiplotTab <- fluidPage(
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            shinydashboard::box(
               title="Plot Parameters",
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
                     width=1,
                     style="padding:0px"
                  ),
                  column(
                     width=4,
                     style="padding:0px",
                     shinyWidgets::radioGroupButtons(
                        inputId="use_exon_names",
                        status="primary",
                        choices=c("coordinates", "exon names"),
                        selected="coordinates",
                        checkIcon = list(
                           yes=icon("ok", lib="glyphicon")
                        ),
                        #width="90%",
                        label="Slider bar measurement"
                     )
                  ),
                  column(
                     width=3,
                     style="padding:0px",
                     shinyWidgets::radioGroupButtons(
                        inputId="include_strand",
                        label="Show coverage by strand:",
                        choices=c("+", "-", "both"),
                        selected="both",
                        status="primary",
                        checkIcon=list(
                           yes=icon("ok", lib="glyphicon"))
                     )
                  ),
                  column(
                     width=3,
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
                     width=1,
                     style="padding:0px"
                  )
               ),
               conditionalPanel(
                  condition="input.use_exon_names == 'exon names'",
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
                  condition="input.use_exon_names == 'coordinates'",
                  tags$strong(textOutput("gene_coords_label",
                     inline=FALSE)),
                  sliderInput(
                     "gene_coords",
                     label=NULL,#"Genome coordinate range",
                     min=28,
                     max=117,
                     #width="80%",
                     value=c(28, 117),
                     step=1,
                     round=TRUE
                  )
               ),
               actionButton(
                  inputId="calc_gene_params",
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
               solidHeader=TRUE,
               closable=FALSE,
               collapsible=TRUE,
               enable_sidebar=TRUE,
               width=12,
               height="100%",
               fluidRow(
                  column(
                     width=12,
                     style="padding:15px",
                     textOutput("sashimitext_output"),
                     shinycssloaders::withSpinner(
                        uiOutput("sashimiplot_output"),
                        type=8
                     )
                  )
               ),
               sidebar_width=25,
               sidebar_start_open=FALSE,
               sidebar_content=tagList(
                  shinyWidgets::sliderTextInput(
                     inputId="panel_height",
                     label="Height per panel:",
                     choices=c(50,100,75,150,200,250,300,400,500),
                     selected=200,
                     grid=TRUE
                  ),
                  shinyWidgets::sliderTextInput(
                     inputId="font_sizing",
                     label="Font sizing:",
                     choices=c("-2 smaller", "-1 smaller", "Default", "+1 larger", "+2 larger"),
                     selected="Default",
                     grid=TRUE
                  ),
                  shiny::sliderInput(
                     inputId="junction_alpha",
                     label="Junction transparency:",
                     min=0.1,
                     max=1.0,
                     step=0.1,
                     value=1.0
                  ),
                  tags$b("Plot Style:"),
                  shinyWidgets::prettyCheckbox(
                     inputId="do_plotly",
                     value=FALSE,
                     icon=icon("check"),
                     status="primary",
                     label="Interactive plot"),
                  conditionalPanel(
                     condition="input.do_plotly == true",
                     nbsp3,
                     nbsp3,
                     shinyWidgets::prettyCheckbox(
                        inputId="enable_highlights",
                        inline=TRUE,
                        value=TRUE,
                        icon=icon("check"),
                        status="primary",
                        label="Enable highlighting"),
                     tags$br(),
                     nbsp3,
                     nbsp3,
                     shinyWidgets::prettyCheckbox(
                        inputId="plotly_legend",
                        inline=TRUE,
                        value=FALSE,
                        icon=icon("check"),
                        status="primary",
                        label="Display filter legend")
                  ),
                  tags$b("Exon Models:"),
                  shinyWidgets::prettyCheckbox(
                     inputId="show_gene_model",
                     value=TRUE,
                     icon=icon("check"),
                     status="warning",
                     label="Show gene-exon model"
                  ),
                  conditionalPanel(
                     condition="input.show_gene_model == true",
                     nbsp3,
                     shinyWidgets::prettyCheckbox(
                        inputId="show_tx_model",
                        inline=TRUE,
                        value=FALSE,
                        icon=icon("check"),
                        status="warning",
                        label="Show transcript-exon model"
                     ),
                     conditionalPanel(
                        condition="input.show_tx_model == true",
                        nbsp3,
                        nbsp3,
                        shinyWidgets::prettyCheckbox(
                           inputId="show_detected_tx",
                           inline=TRUE,
                           value=TRUE,
                           icon=icon("check"),
                           status="warning",
                           label="Detected transcripts only"
                        )
                     )
                  ),
                  br(),
                  tags$b("Axis Settings:"),
                  shinyWidgets::prettyCheckbox(
                     inputId="share_y_axis",
                     value=TRUE,
                     icon=icon("check"),
                     status="success",
                     label="Shared y-axis range"
                  ),
                  conditionalPanel(
                     condition="input.do_plotly == false",
                     sliderInput(
                        inputId="exon_label_size",
                        label="Exon label size",
                        min=2,
                        width="80%",
                        max=20,
                        value=5,
                        step=0.5
                     )
                  )
               )
            )
         )
      )
   );

   # define Samples and Data tab
   samplesdataTab <- fluidPage(
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            DT::DTOutput(
               outputId="files_df"
            )
         )
      )
   );

   # define Sample Selection tab
   sampleselectTab <- fluidPage(
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            shinydashboard::box(
               title="Sample Selection and Order",
               status="warning",
               solidHeader=TRUE,
               collapsible=TRUE,
               width=12,
               "Select samples by dragging them down to the Display Samples section.",
               fluidRow(
                  column(
                     width=7,
                     style="padding:15px",
                     uiOutput("selection_sort")
                  ),
                  column(
                     width=4,
                     style="padding:15px",
                     numericInput(
                        inputId="layout_ncol",
                        label="Number of columns for panel layout",
                        value=1,
                        width="90%",
                        step=1,
                        min=1,
                        max=8
                     )
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
      shinyjs::useShinyjs(),
      setShadow(class="box"),
      setShadow(class="boxPlus"),
      tabItems(
         tabItem(tabName="guides", guidesTab),
         tabItem(tabName="sashimiplot", sashimiplotTab),
         tabItem(tabName="sampleselect", sampleselectTab)
         #tabItem(tabName="samplesdata", samplesdataTab)
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
