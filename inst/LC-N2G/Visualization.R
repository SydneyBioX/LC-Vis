Visual_contour <- function(x,y,z){
  surf_out = Tps_out(x,y,z)
  surf_x = surf_out$x
  surf_y = surf_out$y
  surf_z = surf_out$z
  rownames(surf_z) = surf_x
  colnames(surf_z) = surf_y
  surf_melt = melt(surf_z)
  colnames(surf_melt) = c("x","y","z")
  surf_melt = surf_melt[!is.na(surf_melt[,3]),]
  p = ggplot(surf_melt) + geom_tile(aes(x = x,y = y,fill = z)) + 
    stat_contour(aes(x = x,y = y,z = z,color = ..level..))
  p = direct.label(p,"bottom.pieces")
  return(p)
}

Visual_fix_z <- function(xaxis,yaxis,zaxis2,scalez,scalexy){
  load("Preprocess.RData")
  if (scalez){
    gene_cpm = t(scale(t(data$gene_cpm)))
  }else{
    gene_cpm = data$gene_cpm
  }
  
  if (scalexy){
    MacroNutrition = scale(data$MacroNutrition)
  }else{
    MacroNutrition = data$MacroNutrition
  }
  
  MacroNutrition_pca = prcomp(MacroNutrition,scale. = T)
  MacroNutrition_pc = MacroNutrition_pca$x
  MacroNutrition = cbind(MacroNutrition,MacroNutrition_pc[,1:4])
  x = MacroNutrition[,xaxis]
  y = MacroNutrition[,yaxis]
  z = gene_cpm[zaxis2,]
  D = data.frame(x = x,y = y)
  M = matrix(c(1,0,0,1),nrow = 2)
  Z = z
  p_val = lc_test(D,M,Z)
  
  if(all(x>0)&all(y>0)){
    p = Visual_contour(x,y,z) + xlim(c(0,max(x) + 10)) + ylim(c(0,max(y) + 10))  + 
      scale_fill_gradientn(colours = jet.colors(256)) + theme_classic()
  }else{
    p = Visual_contour(x,y,z) + 
      scale_fill_gradientn(colours = jet.colors(256))+ theme_classic()
  }
  
  p = p + labs(x = xaxis,y = yaxis,fill = zaxis2,title = paste0("gene ",zaxis2," expression vs. ",xaxis," and ",yaxis,"(p-value: ",p_val,").")) +
    geom_point(data = data.frame(x = x,y = y),aes(x = x,y = y),shape = 4)+ theme_classic()
  return(p)
}

Visula_clust_z <- function(xaxis,yaxis,zaxis1,mclust,scalez,scalexy){
  load("Preprocess.RData")
  
  if (scalez){
    gene_cpm = t(scale(t(data$gene_cpm)))
  }else{
    gene_cpm = data$gene_cpm
  }
  
  if (scalexy){
    MacroNutrition = scale(data$MacroNutrition)
  }else{
    MacroNutrition = data$MacroNutrition
  }
  MacroNutrition_pca = prcomp(MacroNutrition,scale. = T)
  MacroNutrition_pc = MacroNutrition_pca$x
  MacroNutrition = cbind(MacroNutrition,MacroNutrition_pc[,1:4])
  
  x = MacroNutrition[,xaxis]
  y = MacroNutrition[,yaxis]
  D = data.frame(x = x,y = y)
  M = matrix(c(1,0,0,1),nrow = 2)
  
  load("Clust_data.RData")
  clust_res = clust_data$clust_res
  clust_color = clust_data$clust_color

  plotres = list()
  for(i in 1:length(zaxis1)){
    gene_group = gene_cpm[clust_color == zaxis1[i],]
    gene_mean = apply(gene_group,2,mean)
    gene_var = apply(gene_group,2,var)
    gene_eigen = prcomp(t(gene_group),scale. = T)$x[,1]
    gene_sel_name = rownames(gene_group)[which.max(abs(cor(t(gene_group),gene_eigen)))]
    gene_sel = gene_group[gene_sel_name,]
    
    p1 = lc_test(D,M,gene_mean)
    plotres[[4*i - 3]] = Visual_contour(x,y,gene_mean) + 
      labs(x = xaxis,y = yaxis,fill = zaxis1[i],title = paste0("mean of cluster ",zaxis1[i],"(p-value: ",p1,").")) + 
      scale_fill_gradientn(colours = jet.colors(256)) +
      geom_point(data = data.frame(x = x,y = y),aes(x = x,y = y),shape = 4) + theme_classic()
      
    p2 = lc_test(D,M,gene_var)
    plotres[[4*i - 2]] = Visual_contour(x,y,gene_var) +
      labs(x = xaxis,y = yaxis,fill = zaxis1[i],title = paste0("variance of cluster ",zaxis1[i],"(p-value: ",p2,").")) + 
      scale_fill_gradientn(colours = jet.colors(256)) +
      geom_point(data = data.frame(x = x,y = y),aes(x = x,y = y),shape = 4) + theme_classic()
    
    p3 = lc_test(D,M,gene_eigen)
    plotres[[4*i - 1]] = Visual_contour(x,y,gene_eigen) +
      labs(x = xaxis,y = yaxis,fill = zaxis1[i],title = paste0("eigen-gene of cluster ",zaxis1[i],"(p-value: ",p3,").")) + 
      scale_fill_gradientn(colours = jet.colors(256)) +
      geom_point(data = data.frame(x = x,y = y),aes(x = x,y = y),shape = 4) + theme_classic()
    
    p4 = lc_test(D,M,gene_sel)
    plotres[[4*i]] = Visual_contour(x,y,gene_sel) +
      labs(x = xaxis,y = yaxis,fill = zaxis1[i],title = paste0(gene_sel_name,"(cluster:",zaxis1[i]," p-value:",p4,").")) + 
      scale_fill_gradientn(colours = jet.colors(256)) +
      geom_point(data = data.frame(x = x,y = y),aes(x = x,y = y),shape = 4) + theme_classic()
  }

  return(plotres)
}

Visual <- function(xaxis,yaxis,zinput,zaxis1,zaxis2,mclust,scalez,scalexy){
  
  if (zinput == "z2"){
    p = Visual_fix_z(xaxis,yaxis,zaxis2,scalez,scalexy)
    return(p)
  } else{
    return(multiplot(plotlist = Visula_clust_z(xaxis,yaxis,zaxis1,mclust,scalez,scalexy),cols = 2))
  }
}

LC_opt_fixedz <- function(zaxis2,scalez,m,Nutrivar){
  load("Preprocess.RData")
  
  # Macronutrition should be scaled to make the distance comparable
  MacroNutrition = scale(data$MacroNutrition)
  D = MacroNutrition[,Nutrivar]
  if (scalez){
    gene_cpm = t(scale(t(data$gene_cpm)))
  }else{
    gene_cpm = data$gene_cpm
  }
  
  Z = gene_cpm[zaxis2,]
  res = as.vector(lc_opt(D,Z,m))
  
  res = colnames(D)[which(res == 1)]
  res = paste(res,collapse = " ,")
  info = paste0("The minimum local consistency nutrition variables in ",m,"-dimenstion for gene ", zaxis2, " is: ",res)
  return(info)
}

LC_opt_clust_z = function(zaxis1,mclust,scalez,m,Nutrivar){
  load("Preprocess.RData")
  load("Clust_data.RData")
  
  MacroNutrition = scale(data$MacroNutrition)
  D = MacroNutrition[,Nutrivar]
  
  if (scalez){
    gene_cpm = t(scale(t(data$gene_cpm)))
  }else{
    gene_cpm = data$gene_cpm
  }
  
  clust_res = clust_data$clust_res
  clust_color = clust_data$clust_color
  info = c()
  for(i in 1:length(zaxis1)){
    gene_group = gene_cpm[clust_color == zaxis1[i],]
    
    Eig_group = prcomp(t(gene_group))$x[,1]
    Z1 = apply(gene_group,2,mean)
    Z2 = apply(gene_group,2,sd)
    Z3 = Eig_group
    
    Z = cbind(Z1,Z2,Z3)
    
    res = as.vector(lc_opt(D,Z,m))
    res = colnames(D)[which(res == 1)]
    res = paste(res,collapse = " ,")
    info[i] = paste0("The minimum local consistency nutrition variables in ",m,"-dimenstion for cluster ", zaxis1[i] ," is: ",res)
  }
  info = paste(info,collapse = ". \n ")
  return(info)
}

LC_opt = function(zinput,zaxis1,zaxis2,mclust,scalez,m,Nutrivar){
  
  if (zinput == "z2"){
    info = LC_opt_fixedz(zaxis2,scalez,m,Nutrivar)
    return(info)
  } else{
    return(LC_opt_clust_z(zaxis1,mclust,scalez,m,Nutrivar))
  }
  
}