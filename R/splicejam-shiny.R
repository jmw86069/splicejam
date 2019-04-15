
#' Launch Sashimi R-shiny application
#'
#' Launch Sashimi R-shiny application
#'
#' This function launches the Sashimi visualization
#' R-shiny app.
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
(...)
{
   #
   shiny::shinyApp(ui=sashimiAppUI,
      server=sashimiAppServer,
      onStart=sashimiAppConstants,
      options=list(width=1200),
      ...);
}
