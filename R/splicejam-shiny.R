
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
#' The most straightforward way to run a new Sashimi R-shiny
#' app is to define `filesDF` and `gtf` in the global environment.
#' The `gtf` is a path or URL to a GTF file (which can be gzipped).
#' This GTF file will be used to derive all related annotation
#' data.
#'
#' * `txdb` -- TranscriptDb from which other objects are derived
#' * `tx2geneDF` -- `data.frame` with transcript-to-gene relationship
#' * `detectedTx` -- if not already defined, all transcripts are
#' used. *Much better to use only a subset of detected transcripts.*
#' * `detectedGenes` -- inferred from `detectedTx`, using `tx2geneDF`.
#' * `flatExonsByGene`, `flatExonsByTx` -- these objects will combine
#' CDS exons and non-CDS exons to represent CDS and UTR regions.
#'
#' Note that if `detectedTx` is not defined, it will use all transcripts
#' at this stage, which can be substantially slower than using only
#' the subset of "observed/detected" transcripts.
#'
#' The first time running `launchSashimiApp()`
#' will populate several R objects in the global environment,
#' and these objects will be re-used during subsequent calls to
#' this function. To make changes in the content, these objects
#' can be edited or deleted so the object is created again.
#' For example, if `detectedTx` is edited, the object
#' `detectedGenes` should be removed so `detectedGenes`
#' will be created again during the next call to `launchSashimiApp()`.
#'
#' The `filesDF` object should be a `data.frame` with colnames
#' `"sample_id"`, `"type"` (with values either `"bw"` or `"junction"`),
#' and `"url"` (a URL or file path to each file.) If coverage
#' or junctions are available in separate files, use the same
#' `sample_id` value for each file. Files with the same `sample_id`
#' value are combined using the sum, after multiplying each file
#' by a value in the optional `"scale_factor"` column.
#'
#' This function calls `sashimiAppConstants()` which does the
#' heavy work of defining or deriving all necessary data objects,
#' then assigns the result to the relevant environment. The default
#' environment is `globalenv()` (also known as `.GlobalEnv`).
#'
#' @family splicejam R-shiny functions
#'
#' @param ... additional arguments are passed to `shiny::shinyApp()`.
#'
#' @examples
#' # Note: disabled for web page examples
#' # launchSashimiApp();
#'
#' @export
launchSashimiApp <- function
(...,
 options=list(width=1200))
{
   ## Temporarily disabled the environment methods below
   if (1 == 2) {
      ## attach dots for local environment
      dots <- list(...);
      if (length(dots) > 0) {
         dots_name <- gsub(" ", "_", paste0("dots ", Sys.time()));
         printDebug("launchSashimiApp(): ",
            "Attaching dots to environment:",
            dots_name);
         attach(dots,
            name=dots_name,
            warn.conflicts=FALSE);
         printDebug("launchSashimiApp(): ",
            "Attached objects:",
            ls(dots_name));
         on.exit(detach(dots_name,
            character.only=TRUE));
      }
   }
   ##
   shiny::shinyApp(ui=sashimiAppUI,
      server=sashimiAppServer,
      onStart=sashimiAppConstants,
      options=options
   );
}
