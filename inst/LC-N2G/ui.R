#options(jave.parameters = "-Xmx1024m")
#pkgs<- c("shiny","shinyjs","shinythemes","WGCNA","dynamicTreeCut","reshape2",
#         "ggplot2","plotly","fields","visNetwork","grid","tidyverse",
#         "DT","directlabels","psych","GA","mclust")
#pkg_type <- c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C")

#if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
#for (i in 1:length(pkgs)){
#    pkg <- pkgs[i]
#    test <- require(pkg,character.only = T)
#    if(!test){
#        if(pkg_type[i] == "C"){
#            install.packages(pkg,ask = F,update = F)
#            require(pkg,character.only = T)
#        }else{
#            BiocManager::install(pkg,ask = F,update = F)
#            library(pkg,character.only = T)
#        }
#    }
#}



source("utils.R")
getwd()
shinyUI(fluidPage(theme = shinythemes::shinytheme("cerulean"),
                  shinyjs::useShinyjs(),
        navbarPage("LC-N2G",
                   id = "LCN2G",
                   #tabPanel(icon("home"),
                    #       fluidRow(tags$img(src="./Fig1-4.jpg",width = "800px",height = "1200px")
                     #          )),

                   tabPanel("Data Preprocessing",
                            sidebarLayout(
                              sidebarPanel(
                                    fileInput("GE","Browse  your gene expression input file (.csv, row sample)", accept = ".csv"),
                                    fileInput("N","Browse  your nutrition input file (.csv, row sample)", accept = ".csv"),
                                    numericInput("mcpm","Filter out gene expression (cpm) below:",5),
                                    numericInput("maxcpm","Filter out gene expression (cpm) above:",500),
                                    numericInput("mvar","Filter out gene expression (cpm) sd below :",0.1),
                                    radioButtons("mnorm","Transformation method",choices = c("cpm" = "cpm","log2 cpm" = "lcpm","none" = "nonorm"),selected = "nonorm"),
                                    sliderInput("mabscor","Filter out gene expression (cpm) maxium absolute correlation with nutrition variables below:",min = 0,max = 1,step = 0.05,value = 0.6),
                                    actionButton("preprocess", "Analysis")
                              ),
                              mainPanel(
                                plotOutput("Preprocess",height = "700px"),
                                textOutput("ngene")
                              )
                            )
                            ),

                   tabPanel("Gene Cluster",
                            sidebarLayout(
                              sidebarPanel(
                                wellPanel(
                                  helpText("Paramerters for clustering using gene-gene relationship:"),
                                  radioButtons("mclust1","Cluster method:",
                                               c("Cluster by gene expression" = "m0",
                                                 "Cluster by gene-gene correlation (WGCAN) " = "m1"),selected = "m1")
                                ),
                                wellPanel(
                                  helpText("Parameters for clustering using gene-nutrition relationship:"),
                                  radioButtons("mclust2","Cluster method:",
                                               c("Cluster by gene-nutrition covariance " = "m0",
                                                 "Cluster by gene-nutrition correlation" = "m1"),selected = "m1")
                                ),
                                wellPanel(
                                  helpText("Other Parameters for clustering"),
                                  sliderInput("alpha","Mix coefficient for two kind of cluster method (0 for gene-gene relationship,1 for gene-nutrition relationship)",
                                              min = 0,max = 1,step = 0.05,value = 0.5),
                                  radioButtons("scaledist","Scaling method for distance matrix:",
                                               c("No scaling" = "m0",
                                                 "Linear transformation scaling" = "m1",
                                                 "Exponential transformation scaling" = "m2"),selected = "m2"),
                                  radioButtons("mcut","Tree cutting method:",
                                               c("Dynamic tree cut" = "c1",
                                                 "By threshold" = "c2")),
                                  textInput("k","Desired number of groups",value = 10),
                                  sliderInput("thresh","Threshold of correlation for network visualization for clustering result:",
                                              min = 0,max = 1,step = 0.05,value = 0.5),
                                  actionButton("cluster", "Analysis")
                                )
                              ),
                              mainPanel(
                                tabsetPanel(
                                  id = "Clust_output",
                                  tabPanel("Dendro",plotOutput("denro"),DT::dataTableOutput("clusttab")),
                                  tabPanel("Table",DT::dataTableOutput("exprtab")),
                                  tabPanel("Summary",visNetwork::visNetworkOutput("vnet",height = "700px"))
                                )
                                )
                            )),

                   tabPanel("Gene Nutrition Visualization",
                            sidebarLayout(
                              sidebarPanel(
                                wellPanel(
                                  helpText("Parameters for z-axis"),
                                  radioButtons("zinput","z-axis chosen by:",
                                               c("Cluster" = "z1","Gene name" = "z2")),
                                  uiOutput("zaxis10"),
                                  textInput("zaxis2","z-axis (gene):",value = "Ucp2"),
                                  checkboxInput("scalez","zscore for z-axis:",value = T)
                                ),

                                wellPanel(
                                  helpText("Parameters for minimum local consistency nutrition variables selection"),
                                  selectInput("Nutrivar","Nutrition variables considered",
                                                      choices = choices_xy(),selected = choices_xy0(),multiple = T),
                                  uiOutput("m0"),
                                  actionButton("visualization0", "Analysis")
                                          ),

                                wellPanel(
                                  helpText("Parameters for 2D visualization of Nutrigenetics"),
                                  uiOutput("xaxis0"),
                                  uiOutput("yaxis0"),
                                  checkboxInput("scalexy","zscore for xy-axis:",value = T),
                                  actionButton("visualization", "Analysis")
                                )

                              ),
                              mainPanel(
                                textOutput("opt"),
                                plotOutput("nutrivisual",height = "700px"))
                            ))

        )))
