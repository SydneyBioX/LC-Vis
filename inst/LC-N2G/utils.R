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

