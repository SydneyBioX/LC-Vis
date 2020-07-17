#' @title shiny app LC-N2G
#' @import shiny
#' @importFrom  shinyjs useShinyjs
#' @import shinythemes
#' @import dynamicTreeCut
#' @import reshape2
#' @import ggplot2
#' @import fields
#' @import visNetwork
#' @import grid
#' @import tidyverse
#' @importFrom  DT dataTableOutput
#' @import directlabels
#' @import GA
#' @import gridExtra
#' @import ggdendro

#' @export run_App

run_App <- function() {
  appDir <- system.file("LC-N2G", package = "LCN2G")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}


