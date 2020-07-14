#' @title shiny app LC-N2G
#' @import shiny
#' @import shinyjs
#' @import shinythemes
#' @importFrom WGCNA labels2colors plotDendroAndColors
#' @import dynamicTreeCut
#' @import reshape2
#' @import ggplot2
#' @import plotly
#' @import fields
#' @import visNetwork
#' @import grid
#' @import tidyverse
#' @import DT
#' @import directlabels
#' @import psych
#' @import GA
#' @import mclust

#' @export run_App

run_App <- function() {
  appDir <- system.file("LC-N2G", package = "LCN2G")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

