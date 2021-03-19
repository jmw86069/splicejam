
#' Sashimi Shiny app UI
#'
#' Sashimi Shiny app UI
#'
#' This function contains the UI for the R-shiny
#' Splicejam Sashimi viewer.
#'
#' The R-shiny app is started by `launchSashimiApp()`, which
#' calls `shiny::shinyApp()`, using arguments `server`, `ui`,
#' `onStart`, and `options`. This function fulfills the
#' argument `ui`.
#'
#' @family splicejam R-shiny functions
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
   header <- dashboardHeader(
      title=tagList("Splicejam Sashimi Viewer",
         icon("map"))
   );
   options("warn"=-1);
   nbsp <- HTML("&nbsp;");
   nbsp3 <- HTML("&nbsp;&nbsp;&nbsp;");

   jam_get <- function
   (name, default, envir, verbose=FALSE, ...)
   {
      dotlist <- list(...);
      if (name %in% names(dotlist)) {
         if (verbose) {
            jamba::printDebug("jam_get(): ",
               "retrieved ", name, " from '", "...", "'");
         }
         return(dotlist[[name]]);
      }
      if (length(envir) > 0) {
         if (is.list(envir)) {
            for (ienv in envir) {
               if (exists(name, envir=ienv)) {
                  if (verbose) {
                     jamba::printDebug("jam_get(): ",
                        "retrieved ", name, " from '", "envir", "'");
                  }
                  return(get(name, envir=ienv));
               }
            }
         } else {
            if (exists(name, envir=envir)) {
               if (verbose) {
                  jamba::printDebug("jam_get(): ",
                     "retrieved ", name, " from '", "envir", "'");
               }
               return(get(name, envir=envir));
            }
         }
      }
      if (exists(name)) {
         if (verbose) {
            jamba::printDebug("jam_get(): ",
               "retrieved ", name, " from '", "get()", "'");
         }
         return(get(name));
      }
      return(default);
   }
   jam_as_logical <- function(x) {
      if (any(c("TRUE", "1", "true", "yes") %in% x)) {
         x <- TRUE;
      } else {
         x <- FALSE;
      }
      return(x);
   }

   min_junction_reads <- jam_get("min_junction_reads",
      100,
      verbose=TRUE,
      ...);
   exon_range_choices_default <- jamba::mixedSort(
      jamba::vigrep("_exon",
      #jamba::provigrep(unique(gsub("exon.*$", "exon", exon_range_selected)),
            values(flatExonsByGene@unlistData)$gene_nameExon));
   if (exists("gene") && length(gene) == 1 && nchar(gene) > 0) {
      exon_range_choices_default <- jamba::mixedSort(
         jamba::vigrep(paste0(gene, "_"),
            exon_range_choices_default));
   }
   exon_range_choices <- jam_get("exon_range_choices",
      exon_range_choices_default,
      verbose=TRUE,
      ...);
   exon_range_selected_default <- c(
      head(exon_range_choices, 1),
      tail(exon_range_choices, 1));
   exon_range_selected <- jam_get("exon_range",
      exon_range_selected_default,
      verbose=TRUE,
      ...);

   #cat(jamba::cPaste(c("exon_range_choices:", exon_range_choices), sep="\n"),
   #   file="debug_output.txt");

   layout_ncol <- jam_get("layout_ncol", 1, verbose=TRUE, ...);
   include_strand <- jam_get("layout_ncol", "both", verbose=TRUE, ...);
   use_exon_names <- jam_get("use_exon_names", "coordinates", verbose=TRUE, ...);
   share_y_axis <- jam_as_logical(jam_get("share_y_axis", TRUE, verbose=TRUE, ...));
   if (!any(c("coordinates", "exon names") %in% use_exon_names)) {
      use_exon_names <- "coordinates";
   }
   show_gene_model <- jam_as_logical(
      jam_get("show_gene_model", TRUE, verbose=TRUE, ...));
   show_tx_model <- jam_as_logical(
      jam_get("show_tx_model", TRUE, verbose=TRUE, ...));
   if (show_tx_model) {
      show_gene_model <- TRUE;
   }
   show_detected_tx <- jam_as_logical(
      jam_get("show_detected_tx", TRUE, verbose=TRUE, ...));

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
      )
   );

   # define Sashimi Plot tab
   sashimiplotTab <- fluidPage(
      fluidRow(
         column(
            width=12,
            style="padding:0px",
            shinydashboardPlus::box(
               title="Plot Parameters",
               status="warning",
               solidHeader=TRUE,
               collapsible=TRUE,
               width=12,
               fluidRow(
                  column(
                     width=1,
                     style="padding:0px"
                  ),
                  column(
                     width=5,
                     style="padding:0px",
                     selectizeInput(
                        label="Select Gene",
                        inputId="gene",
                        selected=default_gene,
                        choices=detectedGenes,
                        options=list(maxOptions=100),
                        multiple=FALSE
                     )
                  ),
                  column(
                     width=1,
                     style="padding:0px"
                  ),
                  column(
                     width=4,
                     style="padding:0px",
                     selectizeInput(
                        label="Genes to Search",
                        inputId="search_genelist",
                        selected="Detected",
                        choices=c("Detected", "All"),
                        multiple=FALSE
                     )
                  ),
                  column(
                     width=1,
                     style="padding:0px"
                  )
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
                        selected=use_exon_names,
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
                        label="Show coverage by strand",
                        choices=c("+", "-", "both"),
                        selected=include_strand,
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
                        value=min_junction_reads,
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
                  shinyWidgets::sliderTextInput(
                     inputId="exon_range",
                     label="Gene exon range",
                     grid=TRUE,
                     force_edges=TRUE,
                     choices=exon_range_choices,
                     selected=exon_range_selected
                     #selected=c("exon1", "exon3")
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
            shinydashboardPlus::box(
               title="Sashimi Plot",
               status="primary",
               solidHeader=TRUE,
               closable=FALSE,
               collapsible=TRUE,
               ## icons: https://adminlte.io/themes/AdminLTE/pages/UI/icons.html
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
                     ),
                     plotlyOutput("plotly_blank")
                  )
               ),
               #sidebar_icon="fa fa-cog",
               #sidebar_icon="fa fa-bar-chart",
               # glyphicon glyphicon-cog
               sidebar=boxSidebar(
                  id="sashimi_sidebar",
                  width=25,
                  startOpen=FALSE,
                  #background="#FFE5C8",
                  background="#AAAAAA",
                  icon=shiny::icon("gear"),
                  tagList(
                     shinyWidgets::sliderTextInput(
                        inputId="panel_height",
                        label="Height per panel:",
                        choices=c(50,75,100,150,200,250,300,400,500),
                        selected=200,
                        grid=TRUE
                     ),
                     shinyWidgets::sliderTextInput(
                        inputId="font_sizing",
                        label="Font sizing:",
                        choices=c("-4 smaller",
                           "-3 smaller",
                           "-2 smaller",
                           "-1 smaller",
                           "Default",
                           "+1 larger",
                           "+2 larger",
                           "+3 larger",
                           "+4 larger"),
                        selected="Default",
                        grid=TRUE
                     ),
                     shiny::sliderInput(
                        inputId="junction_alpha",
                        label="Junction transparency:",
                        min=0.1,
                        max=1.0,
                        step=0.1,
                        value=0.7
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
                           value=FALSE,
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
                     tags$b("Axis Settings:"),
                     shinyWidgets::prettyCheckbox(
                        inputId="share_y_axis",
                        value=share_y_axis,
                        icon=icon("check"),
                        status="success",
                        label="Shared y-axis range"
                     ),
                     tags$b("Exon Models:"),
                     shinyWidgets::prettyCheckbox(
                        inputId="show_gene_model",
                        value=show_gene_model,
                        icon=icon("check"),
                        status="warning",
                        label="Show gene model"
                     ),
                     conditionalPanel(
                        condition="input.show_gene_model == true",
                        nbsp3,
                        shinyWidgets::prettyCheckbox(
                           inputId="show_tx_model",
                           inline=TRUE,
                           value=show_tx_model,
                           icon=icon("check"),
                           status="warning",
                           label="Show transcripts"
                        ),
                        conditionalPanel(
                           condition="input.show_tx_model == true",
                           nbsp3,
                           nbsp3,
                           shinyWidgets::prettyCheckbox(
                              inputId="show_detected_tx",
                              inline=TRUE,
                              value=show_detected_tx,
                              icon=icon("check"),
                              status="warning",
                              label="Detected only"
                           )
                        )
                     ),
                     br(),
                     conditionalPanel(
                        condition="input.do_plotly == false",
                        shinyWidgets::sliderTextInput(
                           inputId="exon_label_size",
                           label="Exon font sizing:",
                           choices=c("-4 smaller",
                              "-3 smaller",
                              "-2 smaller",
                              "-1 smaller",
                              "Default",
                              "+1 larger",
                              "+2 larger",
                              "+3 larger",
                              "+4 larger"),
                           selected="Default",
                           grid=TRUE
                        )
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
               fluidRow(
                  column(
                     width=7,
                     style="padding:15px",
                     paste("Select samples in the order they should appear in the figure.",
                        "The table will re-order itself to place selected samples",
                        "at the top. De-select samples to remove them from the figure.",
                        "Click 'Update Sashimi Plots' to create a new figure."),
                     DT::dataTableOutput('samplesdf')
                  ),
                  column(
                     width=4,
                     style="padding:15px",
                     numericInput(
                        inputId="layout_ncol",
                        label="Number of columns for panel layout",
                        value=layout_ncol,
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

   dp <- shinydashboardPlus::dashboardPage(
      header,
      sidebar,
      body,
      skin="blue");
   dp;

}
