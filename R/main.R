#' Start web application
#'
#' Launch Shiny application in default web browser
#'
#' @export
#'
#' @examples \dontrun{
#' webstart()
#' }
webstart <- function() {
  shiny::runApp(system.file("shiny", package = "WebCytoMetry"), launch.browser = T)
}
