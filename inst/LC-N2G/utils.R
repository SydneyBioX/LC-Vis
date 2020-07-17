labels2colors <- function (labels, zeroIsGrey = TRUE, colorSeq = NULL, naColor = "grey",
          commonColorCode = TRUE)
{
  if (is.null(colorSeq))
    colorSeq = standardColors()
  if (is.numeric(labels)) {
    if (zeroIsGrey)
      minLabel = 0
    else minLabel = 1
    if (any(labels < 0, na.rm = TRUE))
      minLabel = min(c(labels), na.rm = TRUE)
    nLabels = labels
  }
  else {
    if (commonColorCode) {
      factors = factor(c(as.matrix(as.data.frame(labels))))
      nLabels = as.numeric(factors)
      dim(nLabels) = dim(labels)
    }
    else {
      labels = as.matrix(as.data.frame(labels))
      factors = list()
      for (c in 1:ncol(labels)) factors[[c]] = factor(labels[,
                                                             c])
      nLabels = sapply(factors, as.numeric)
    }
  }
  if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
    nRepeats = as.integer((max(labels) - 1)/length(colorSeq)) +
      1
    warning(paste("labels2colors: Number of labels exceeds number of avilable colors.",
                  "Some colors will be repeated", nRepeats, "times."))
    extColorSeq = colorSeq
    for (rep in 1:nRepeats) extColorSeq = c(extColorSeq,
                                            paste(colorSeq, ".", rep, sep = ""))
  }
  else {
    nRepeats = 1
    extColorSeq = colorSeq
  }
  colors = rep("grey", length(nLabels))
  fin = !is.na(nLabels)
  colors[!fin] = naColor
  finLabels = nLabels[fin]
  colors[fin][finLabels != 0] = extColorSeq[finLabels[finLabels !=
                                                        0]]
  if (!is.null(dim(labels)))
    dim(colors) = dim(labels)
  colors
}

plotDendroAndColors <- function (dendro, colors, groupLabels = NULL, rowText = NULL,
                                 rowTextAlignment = c("left", "center", "right"),
                                 rowTextIgnore = NULL, textPositions = NULL, setLayout = TRUE,
                                 autoColorHeight = TRUE, colorHeight = 0.2, colorHeightBase = 0.2,
                                 colorHeightMax = 0.6, rowWidths = NULL, dendroLabels = NULL,
                                 addGuide = FALSE, guideAll = FALSE, guideCount = 50, guideHang = 0.2,
                                 addTextGuide = FALSE, cex.colorLabels = 0.8, cex.dendroLabels = 0.9,
                                 cex.rowText = 0.8, marAll = c(1, 5, 3, 1), saveMar = TRUE,
                                 abHeight = NULL, abCol = "red", ...)
{
  oldMar = par("mar")
  if (!is.null(dim(colors))) {
    nRows = dim(colors)[2]
  }
  else nRows = 1
  if (!is.null(rowText))
    nRows = nRows + if (is.null(textPositions))
      nRows
  else length(textPositions)
  if (autoColorHeight)
    colorHeight = colorHeightBase + (colorHeightMax - colorHeightBase) *
      (1 - exp(-(nRows - 1)/6))
  if (setLayout)
    layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight,
                                             colorHeight))
  par(mar = c(0, marAll[2], marAll[3], marAll[4]))
  plot(dendro, labels = dendroLabels, cex = cex.dendroLabels,
       ...)
  if (addGuide)
    addGuideLines(dendro, count = if (guideAll)
      length(dendro$height) + 1
      else guideCount, hang = guideHang)
  if (!is.null(abHeight))
    abline(h = abHeight, col = abCol)
  par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
  plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels,
                     rowText = rowText, rowTextAlignment = rowTextAlignment,
                     rowTextIgnore = rowTextIgnore, textPositions = textPositions,
                     cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
  if (saveMar)
    par(mar = oldMar)
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = T)
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

choices_xy0 <- function(data){
  load("Preprocess.RData")
  MacroNutrition = data$MacroNutrition
  choice = colnames(MacroNutrition)
  names(choice) = choice
#  choice = c(choice,"PC1" = "PC1","PC2" = "PC2","PC3" = "PC3","PC4" = "PC4")
  return(choice)
}

l_choices_xy <- function(data){
  return(length(choices_xy()))
}

choices_xy <- function(data){
  choice = choices_xy0()
  choice = c(choice,"PC1" = "PC1","PC2" = "PC2","PC3" = "PC3","PC4" = "PC4")
  return(choice)
}

choices_z <-function(data){
  load("Clust_data.RData")
  clust_color = clust_data$clust_color
  choice = names(table(clust_color))
  names(choice) = choice
  return(choice)
}

Tps_out = function(x,y,z)
{
  sumframe<-structure(list(xvalue = x, yvalue = y, zvalue = z), .Names = c("xvalue", "yvalue", "zvalue"), class = "data.frame")
  surf<-fields::Tps(cbind(sumframe$xvalue, sumframe$yvalue), sumframe$zvalue, lambda=0.01)

  surf.out=fields::predictSurface(surf)

  return(surf.out)
}

lc = function(D,M,Z,lmd = NULL){

  eps = 0.001

  if (is.vector(Z)){
    Z = matrix(Z,ncol = 1)
  }else{
    Z = as.matrix(Z)
  }

  D = scale(D)
#  Z = scale(Z)

  D = as.matrix(D)


  if(is.null(lmd)){
    adp_weight = 1/sqrt(2)*(apply(Z,2,var) + eps)
    lmd = adp_weight
  }

#  print(Z)
  D_z  = abs(outer(Z[,1],Z[,1],"-"))


  D_M = exp(as.matrix(-dist(D %*% M)^2))
 # print(D_M)
 # print(D_z)

  R = sum(D_M * abs(D_z))
  return(R)
}


lc_test = function(D,M,Z,n = 100){

  D = scale(D)
  Z = scale(Z)
  lc_0 = lc(D,M,Z)
  lc = c()
  for(i in 1:n){
    Z1 = sample(Z,replace = F)
    lc[i] = lc(D,M,Z1)
  }
  p = sum(lc<lc_0)/n
  return(p)
}

lc_opt = function(D,Z,k = 2,miter = 50){
  # k is number of 1 in metric
  D = scale(D)
#  Z = scale(Z)
  flag = choose(ncol(D),k) <= 10000
  if(flag){
    res = lc_exhaust(D,Z,k,type = 2)

    return(res)
  }else{
  lc_GA = function(M,D,Z,k){
    M = diag(M)
    res = -lc(D,M,Z) - 10000*abs(tr(M) - k)
    return(res)
  }

  lc_GA = partial(lc_GA,D = D,Z = Z, k = k)
  GA = GA::ga(type = "binary",fitness = lc_GA,nBits = ncol(D),maxiter = miter)

  return(GA@solution)
  }
}

lc_exhaust <- function(D,Z,k,type = 2,isorder = T){
  D = scale(D)
  searchspace = combn(colnames(D),k)
  res = c()
  for(i in 1:ncol(searchspace)){
    D1 = D[,searchspace[,i]]
    res[i] = lc(D1,diag(k),Z)
  }
  if (type == 1){
    # return all result, for histogram
    if (isorder){
      t = order(res,decreasing = F)
      res = res[t]
      searchspace = searchspace[,t]

      LC = list(searchspace = searchspace,lc = res)
    }else
    {
      LC = list(searchspace = searchspace,lc = res)
    }


    return(LC)
  }else{

    t = which.min(res)
    label = searchspace[,t]
    out = rep(0,ncol(D))
    names(out) = colnames(D)

    out[label] = 1
    return(out)
    }
}




bl <- function(n,k,q = 1000){
  res = c()
  for(i in 1:q){
    D = matrix(0,nrow = n,ncol = k+1)
    for (j in 1:(k+1)){
      D[,j] = rnorm(n)
    }
    D = scale(D)
    res[i] = lc(D[,1:k],diag(k),D[,(k+1)])
  }
  return(res)
}

