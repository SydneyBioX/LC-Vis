source("utils.R")
source("preprocessing.R")
source("Cluster.R")
source("Visualization.R")

shinyServer(function(input, output,session) {
  mcpm = reactive(as.numeric(input$mcpm))
  mvar = reactive(as.numeric(input$mvar))
  mabscor = reactive(as.numeric(input$mabscor))
  maxcpm = reactive(as.numeric(input$maxcpm))

#  showModal(modalDialog(
#    footer = "./readme.html"
#    easyClose = TRUE
#  ))

  #preprocessing

  output$DemoN<- downloadHandler(filename = function() {
    paste("Nutrition_Demo", ".csv", sep="")
  },
  content = function(file) {
    dat1 = N_d()
    write.csv(dat1, file)
  }
  )

  output$DemoG<- downloadHandler(filename =function() {
    paste("Gene Expression_demo", ".csv", sep="")
  },
  content = function(file) {
    dat1 = GE_d()
    write.csv(dat1, file)
  }
  )

  observeEvent(input$preprocess,
               {output$Preprocess <-
                 renderPlot(
                   {
                     input$preprocess
                     isolate(preprocessing(input$GE$datapath, input$N$datapath, mcpm(),maxcpm(),mvar(),mabscor(),input$mnorm))}
                   )
               output$ngene <- renderText({
                 paste0("Number of gene after preprocessing is ",numgene(),".")
               })})

  #js for cutree
  observe({shinyjs::toggleState("k", input$mcut == "c2")})

  #cluster
  observeEvent(input$cluster,
               {output$denro <-
                 renderPlot(
                   {
                     input$cluster
                     isolate(Clust_dendro(data,input$mclust1,input$mclust2,input$alpha,input$mcut,input$scaledist,input$k))}
                 )
               output$clusttab <- DT::renderDataTable(
                 {
                   input$cluster
                   isolate(DT::datatable(Clust_tab(data,input$mclust,input$mcut,input$k)))}
               )
               output$exprtab <- DT::renderDataTable(
                 {
                   input$cluster
                   isolate(DT::datatable(Expr_tab(data,input$mclust,input$mcut,input$k)))}
               )

               output$vnet <- visNetwork::renderVisNetwork(
                 {
                   input$cluster
                   isolate(Clust_summary(data,input$mclust1,input$mclust2,input$alpha,input$mcut,input$scaledist,input$thresh,input$k))}
               )
              })

  #Visualization

  #js for z-axis
  observe({shinyjs::toggleState("zaxis1", input$zinput == "z1")})
  observe({shinyjs::toggleState("zaxis2", input$zinput == "z2")})

  output$xaxis0 <- renderUI({
    selectInput("xaxis","x-axis:",choices = c(input$Nutrivar))
  })


  output$yaxis0 <- renderUI({
    selectInput("yaxis","y-axis:",choices = input$Nutrivar[input$Nutrivar != input$xaxis])
  })

  output$m0 <- renderUI({
    sliderInput("m","Number of Nutrition variables for visualization",min = 1,max = length(input$Nutrivar),step = 1,value = 2)
  })

  output$zaxis10 <- renderUI({
    selectInput("zaxis1","z-axis(cluster):",
                choice = choices_z(),multiple = T)
  })
  #lc-opt
  observeEvent(input$visualization0,
               {output$opt <- renderText(
                 {
                   input$visualization0
                   isolate(LC_opt(input$zinput,input$zaxis1,input$zaxis2,input$mclust,input$scalez,input$m,input$Nutrivar))}
                   #LC_opt_fixedz(input$zaxis2,input$scalez,input$m,input$Nutrivar)}
               )}
  )

#  output$opt <- renderText({
#    input$visualization0
#    isolate(LC_opt_fixedz(input$zaxis2,input$scalez,input$m,input$Nutrivar))
#  })

  #lc-p

  observeEvent(input$visualization,
               {output$nutrivisual <- renderPlot(
                 {
                   input$visualization
                   isolate(Visual(input$xaxis,input$yaxis,input$zinput,input$zaxis1,input$zaxis2,input$mclust,input$scalez,input$scalexy))
                   }
               )}
               )
  })
