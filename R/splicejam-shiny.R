
#' Launch Sashimi R-shiny application
#'
#' Launch Sashimi R-shiny application
#'
#' This function launches the Sashimi visualization
#' R-shiny app.
#' 
#' The key input data should be supplied through argument
#' `envir` which should be prepared using
#' `sashimiDataConstants()`. It will prepare all the
#' necessary data to use.
#' 
#' ## Custom Parameters
#' 
#' The environment `envir` may have additional visual
#' customizations, which recognize the following names:
#' 
#' * 'panel_height', default `200` pixels
#' * 'gene_panel_height', default `400` pixels
#' * 'layout_ncol', default 1
#' * 'min_junction_reads', default `100`
#' * 'share_y_axis', default `FALSE`
#' * 'gene_coords_default', default NULL
#' * 'label_junctions', default `TRUE`
#' * 'show_gene_model', default `TRUE`
#' * 'show_tx_model', default `TRUE`
#' * 'show_detected_tx', default `TRUE`
#' * 'font_sizing', default "Default", note it must match available
#' choices in the UI: "-4 smaller" to "+4 larger".
#' * 'exon_font_sizing', default "Default", see comment for
#' 'font-sizing', the same applies.
#' * 'junction_alpha', default `0.7` for 70% opacity
#' * 'junction_arc_factor', default "Default", note it must match
#' available choices in the UI: "-2 flat" to "+3 higher".
#' * 'junction_arc_minimum', default `500`
#' * 'aboutExtra', default is html tag data suitable to render
#' in an HTML section of the R-shiny web page. This entry
#' is useful to describe the data or experiment underlying
#' the Splicejam visualizations shown.
#' 
#' The environment is prepared by calling `sashimiAppConstants()`
#' which also defines `aboutExtra` and other UI elements.
#' It calls `sashimiDataConstants()` which will verify the
#' environment is usable.
#'
#' ## Running within shiny-server
#'
#' Create a file `'app.R'` and include lines which create
#' the environment, for example:
#' 
#' ```R
#' # Prepare data
#' sjenv <- sashimiDataConstants(gtf=gtf, filesDF=filesDF)
#' 
#' # Add any customizations
#' sjenv$default_gene <- "Gria1";
#' sjenv$panel_height <- 100;
#' sjenv$gene_panel_height <- 200;
#' 
#' # Launch the app
#' launchSashimiApp()
#' ```
#' 
#' This file should be suitable to run the shiny server.
#' 
#' Note that the default will load Farris et al data into the
#' environment, controlled by `empty_uses_farrisdata=TRUE`.
#' 
#' ## Other data notes
#' 
#' The `filesDF` object should be a `data.frame` with
#' at least three colnames:
#'
#' * `"sample_id"`
#' * `"type"` (with values either `"bw"` or `"junction"`)
#' * `"url"` (a URL or file path to each file.)
#' 
#' In addition, it recognizes `"scale_factor"` with `numeric`
#' values used to scale (by multiplication) the corresponding
#' scores.
#' 
#' If coverage or junctions are available in multiple files for
#' the same `'sample_id'`, for example when using 
#' replicates, use the same `sample_id` for each file, and
#' the coverage and junctions will be combined by taking the sum.
#' The optional `"scale_factor"` is applied to each file first,
#' which permits adjusting the files upfront.
#' 
#' The reason for using the sum is that the Splicejam visualization
#' is intended to show the evidence supporting particular exons
#' and junctions, and this evidence is built up across replicates.
#'
#' For more direct control over the data preparation, including
#' `tx2geneDF`, `detectedTx`, `exonsByGene`, and `flatExonsByGene`,
#' see `sashimiAppConstants()` which calls `sashimiDataConstants()`,
#' both of these functions return an environment that contains the
#' required data.
#'
#' When the R-shiny app is created, the `ui` and `server` components
#' have their environments set to `envir` - so their context will
#' include the variables defined in that environment.
#'
#' ## Troubleshooting
#'
#' 1. Error `"SSL peer certificate or SSH remote key was not OK"` or
#' `"SSL certificate problem: unable to get local issuer certificate"`
#' seen in the console output of the R shiny app.
#'
#'    * This error occurs when `launchSashimiApp()` references
#'    coverage data on an external web server, for example when
#'    using `farrisdata` example data. The remote web server
#'    certificates used with `https` web address are using
#'    a certificate whose issuer is not recognized.
#'    * It means the authority that signed the certificate
#'    is not recognized as an approved authority.
#'    * There may be two workarounds to this issue:
#'
#'    1. If the certificate authority should be approved, sometimes
#'       on linux hosts it is enough to add certificate extensions,
#'       for example on Ubuntu, or Ubuntu Docker images, one may run
#'       ```
#'       apt-get ca-certificates
#'       # or
#'       sudo apt-get ca-certificates
#'       ```
#'       Other linux systems may require a different installation.
#'    2. The certificate verification can be skipped temporarily.
#'       Define this option within the same R session:
#'       ```
#'       httr::set_config(httr::config(ssl_verifypeer=FALSE))
#'       ```
#'       Of course, this option should only be used in trusted
#'       circumstances, for example when accessing a local and trusted
#'       host.
#'
#' 2. Error `"covNames must be in colnames(GenomicRanges::values(gr))"`,
#' or other error related to `names(covNames)` or `names(covNamesL)`.
#'
#'    * These errors ultimately mean there was no available coverage data,
#'    which is usually caused by inability to access the coverage
#'    file itself. The web server may be offline, or may be denying
#'    connection. The path to the file may be incorrect.
#'    * Typically the coverage data is obtained from memoise cached data,
#'    and only when the cache is not available, or somehow incorrect,
#'    the data is retrieved from the file or remote server.
#'
#' @family splicejam R-shiny functions
#'
#' @returns `shiny::shiny.appobj` which is a Shiny App object,
#'    suitable to run an R-shiny app by printing to console.
#'
#' @return output from `shiny::shinyApp()` which is an object of
#'    class "shiny.appobj", whose default print method is to run
#'    the app.
#'
#' @param envir `environment` default `parent.frame()` uses the environment
#'    which called this function.
#'    The environment contains data needed for sashimi plots.
#'    If `envir=NULL` by default it will use `parent.frame()`.
#'    Otherwise  call `sashimiDataConstants()` or `sashimiAppConstants()`
#'    which returns an `environment` with the necessary data objects,
#'    and that can be passed to this function.
#' @param options `list` of R-shiny app options, for example two
#'    common options are: `host` to indicate the host or IP address
#'    the R-shiny app will bind to respond to requests; and
#'    `port` for the port number. For example:
#'    `launchSashimiApp(options=list(host="0.0.0.0", port=8080))`
#'    where `port="0.0.0.0"` will listen to any request sent
#'    to the current machine, whether by host name, or any valid
#'    IP address; and `port=8080` will only listen to port 8080.
#' @param ... additional arguments are passed to `sashimiAppConstants()`.
#'
#' @examples
#' # Note: disabled for web page examples
#' # launchSashimiApp();
#'
#' @export
launchSashimiApp <- function
(...,
 envir=parent.frame(),
 options=list(width=1200),
 verbose=FALSE)
{

   # retrieve an environment that contains the required data
   envir <- sashimiAppConstants(envir=envir,
      verbose=verbose,
      ...);

   # Define these functions specifically so we can set the environment.
   # - Now the environment will contain the data obtained above.
   ui <- sashimiAppUI;
   server <- sashimiAppServer;
   environment(ui) <- envir;
   environment(server) <- envir;

   if (length(options) == 0) {
      options <- list();
   } else if (!is.list(options)) {
      options <- as.list(options);
   }
   if (!"width" %in% names(options)) {
      options$width <- 1200;
   }
   ##
   shiny::shinyApp(ui=ui,
      server=server,
      options=options
   );
}
