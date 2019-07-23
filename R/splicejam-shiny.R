
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
#' Typically, the first time running `launchSashimiApp()`
#' will populate several R objects in the global environment,
#' and they will be re-used during subsequent calls to
#' this function. These R objects are intended to be
#' available to review, for possible troubleshooting
#' operations, or to update the data as needed.
#'
#' For example, the `data.frame` object `filesDF` contains
#' a colname `"scale_factor"` used to help normalize the
#' visible height on the y-axis of individual entries.
#' The `"scale_factor"` values can be edited, then the
#' `launchSashimiApp()` function will use the existing
#' data with the new `"scale_factor"` values.
#'
#' More about using R environments:
#'
#' One comment about the environment. The `sashimiAppConstants()`
#' function uses the `exists()` function to check for existing
#' objects, and the `<<-` operator, which updates
#' R objects, both these methods search up the parent
#' environment chain to find a matching object name.
#'
#' For most cases, the global environment is used, which
#' can be convenient for creating and updating R objects
#' for use outside the R-shiny sashimi app. However, to
#' avoid populating objects in the global environment,
#' variables can be passed through the `dots` argument `...`
#' when calling `launchSashimiApp()`.
#'
#' For example `launchSashimiApp(filesDF=farris_sashimi_files_df)`
#' will define a local variable `filesDF` in the context
#' of `launchSashimiApp()`. The `exists()` function will recognize
#' this object and use it as-is. Any updates to `filesDF` using
#' `<<-` will update the object in the local function environment.
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
      options=list(width=1200)
   );
}
