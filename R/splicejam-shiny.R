
#' Launch Sashimi R-shiny application
#'
#' Launch Sashimi R-shiny application
#'
#' This function launches the Sashimi visualization
#' R-shiny app.
#'
#' The R objects required to prepare sashimi plots are
#' defined by the function `sashimiAppConstants()`, which
#' documents each R object required and how it is used.
#'
#' In general, if you have any of `txdb`, `tx2geneDF`, `filesDF`,
#' `gtf` (optional) defined, this function will populate the
#' other objects for you as needed. It looks into the `parent.frame()`
#' which is the environment calling `launchSashimiApp()`,
#' unless you provide specific environment `envir`.
#'
#' It will also populate the new objects inside this environment,
#' unless you pass `assign_global=FALSE`.
#'
#' The `sashimiAppConstants()` function returns an `environment`
#' in which the required sashimi plot data is stored and
#' used by the R-shiny app. The default environment is `globalenv()`
#' however any custom environment can be used, for example
#' with `myenv <- new.env()`.
#'
#' The most straightforward way to run a new Sashimi R-shiny
#' app is to define: `filesDF` and `gtf`;
#' or: `filesDF`, `txdb`, `tx2geneDF`.
#' You may define them in the current environment,
#' or in a custom environment
#' and pass the environment using argument `envir`.
#'
#' When providing `gtf` it will create `txdb` and `tx2geneDF`
#' as needed. If providing `txdb` and `tx2geneDF` it will use
#' that data asis. It caches the necessary intermediate files
#' when `use_memoise=TRUE` (default) so the time spent is
#' usually first-time cost only.
#'
#' If no data are defined, the default `empty_uses_farrisdata=TRUE`
#' will cause the app to use Farris et al 2019 data.
#'
#' ## Running in shiny-server
#'
#' Create a file `app.R` and include lines which define the
#' data above (gtf, filesDF; or txdb, tx2geneDF, filesDF)
#' then call `launchSashimiApp()` at the end.
#'
#' If no data are defined, the default `empty_uses_farrisdata=TRUE`
#' will cause the app to use Farris et al 2019 data
#' together with the `farrisdata` R package, and will run the
#' Shiny App using that data.
#'
#' ## Data Preparation
#'
#' The data derived from the GTF file is listed below. Any data object
#' that already exists in the `environment` is used in subsequent steps:
#'
#' * `txdb` - TranscriptDb from which `exonsByGene` and `exonsByTx`
#'    are derived.
#' * `tx2geneDF` - `data.frame` with transcript-to-gene relationship,
#'    with colnames `"gene_name"` and `"transcript_id"`. See
#'    `makeTx2geneFromGtf()` for details.
#' * `detectedTx` - `character` vector of detected transcripts, used
#'    to match `tx2geneDF$transcript_id`. When `detectedTx` is NULL,
#'    all entries in `tx2geneDF` are used. Note that we found it is
#'    *much better to use only a subset of detected transcripts*,
#'    mainly because many GTF sources include a large number of potential
#'    alternative isoforms, many of which have no supported evidence in
#'    any one given cell type. See `defineDetectedTx()` for one method
#'    to define detected transcripts.
#' * `detectedGenes` - `character` vector of genes that match
#'    `tx2geneDF$gene_name`. When `detectedGenes` is NULL, it is
#'    inferred using `detectedTx` and `tx2geneDF$transcript_id`.
#' * `exonsByTx`, `cdsByTx` - derived from `txdb` and annotated to include
#'    values from `tx2geneDF$gene_name`.
#' * `flatExonsByGene`, `flatExonsByTx` - `GRangesList` objects derived
#'    from `exonsByGene` and `exonsByTx`, using `detectedTx`.
#'    They also use `cdsByTx` to indicate coding regions (CDS) of exons.
#'
#' Note that if `detectedTx` is not defined, it will use all transcripts
#' at this stage, which can be substantially slower than using only
#' the subset of "observed/detected" transcripts.
#'
#' An alternative is to supply one `detectedGenes` gene value, which will prepare
#' only one gene for `flatExonsByGene` in the R-shiny app. However, the
#' R-shiny app has the option to search all non-detected genes, which
#' are prepared one by one inside the R-shiny app. This process is slightly
#' slower when using the app by a few seconds, and will use all transcripts
#' for `detectedTx`.
#'
#' The `filesDF` object should be a `data.frame` with at least three colnames:
#'
#' * `"sample_id"`
#' * `"type"` (with values either `"bw"` or `"junction"`)
#' * `"url"` (a URL or file path to each file.)
#'
#' If coverage or junctions are available in multiple files, for example
#' sequencing replicates, use the same `sample_id` for each file, and
#' the coverage and junctions will be combined using the sum, after
#' multiplying an optional `"scale_factor"` to each file.
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
