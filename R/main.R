#' Start web application
#'
#' Launch Shiny application in default web browser
#'
#' @export
#'
#' @examples \dontrun{
#' start()
#' }
start <- function() {
  shiny::runApp(system.file("shiny", package = "WebCytoMetry"),
                launch.browser = T)
}
