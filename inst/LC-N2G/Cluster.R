distscaling <- function(clust_dist,scaledist,sgm = NULL){
  clust_dist = as.matrix(clust_dist)
  if (scaledist == "m0"){
    clust_dist = clust_dist
  }else if(scaledist == "m1"){
    clust_dist = clust_dist/max(clust_dist)
  }else{
    if (is.null(sgm)){
      clust_dist1 = as.matrix(clust_dist)
      sgm = mean(clust_dist1[upper.tri(clust_dist1)])
    }
    clust_dist = 1 - exp(-(clust_dist/(sqrt(2)*sgm))^2)
  }
  return(clust_dist)
}

corClust1 <- function(data,mcut,scaledist){
  load("Preprocess.RData")
  gene_cor = data$gene_cor
  dist_cor = stats::dist(t(gene_cor))
  dist_cor = distscaling(dist_cor,scaledist)
  return(dist_cor)
}

covClust1 <- function(data,mcut,scaledist){
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  MacroNutrition = data$MacroNutrition
  MacroNutrition_scale = scale(MacroNutrition)
  gene_cov = cov(MacroNutrition_scale,t(gene_cpm))
  dist_cov = stats::dist(t(gene_cov))
  dist_cov = distscaling(dist_cov,scaledist)
  return(dist_cov)
}

expClust0 <- function(data,mcut,scaledist){
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  dist_gene = stats::dist(gene_cpm)
  dist_gene = distscaling(dist_gene,scaledist)
  return(dist_gene)
}

corClust0 <- function(data,mcut,scaledist){
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  dist_cor = 1 - abs(cor(t(gene_cpm)))
  dist_cor = distscaling(dist_cor,scaledist)
  return(dist_cor)
}

hybClust <- function(data,mclust1,mclust2,alpha,scaledist){
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  gene_cor = data$gene_cor
  MacroNutrition = data$MacroNutrition
  gene_entrz = data$gene_entrz
  MacroNutrition_scale = scale(MacroNutrition)
  if(mclust1 == "m0"){
    clust_dist1 = expClust0(data,mcut,scaledist)
  }else{
    clust_dist1 = corClust0(data,mcut,scaledist)
  }
  if(mclust2 == "m0"){
    clust_dist2 = covClust1(data,mcut,scaledist)
  }else{
    clust_dist2 = corClust1(data,mcut,scaledist)
  }
  clust_dist = alpha*clust_dist2 + (1 - alpha)*clust_dist1
  return(clust_dist)
}

Clust_dendro <- function(data,mclust1,mclust2,alpha,mcut,scaledist,k = NULL){
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  gene_cor = data$gene_cor
  MacroNutrition = data$MacroNutrition
  gene_entrz = data$gene_entrz
  
  clust_dist = as.dist(hybClust(data,mclust1,mclust2,alpha,scaledist))
  
  clust = hclust(clust_dist)
  
  if (mcut == "c2")
  {
    clust_res = cutree(clust,k = k)
  }else{
    clust_res = cutreeDynamicTree(clust)
  }
  names(clust_res) = colnames(gene_cor)
  clust_color = labels2colors(clust_res)
  
  clust_data = list(clust = clust,clust_res = clust_res,clust_color = clust_color,clust_dist = clust_dist)
  save(clust_data, file = "Clust_data.RData")
  p1 = plotDendroAndColors(clust, clust_color,
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  return(p1)
}

Clust_tab <- function(data,mclust1,mclust2,alpha,mcut,scaledist,k = NULL){
  load("Clust_data.RData")
  clust = clust_data$clust
  clust_res = clust_data$clust_res
  clust_color = clust_data$clust_color
  color_tab = table(clust_color)
  res = t(matrix(color_tab))
  colnames(res) = names(color_tab)
  rownames(res) = "Cluster Size"
  res = as.data.frame(res)
  return(res)
}
  
Expr_tab <- function(data,mclust1,mclust2,alpha,mcut,scaledist,k = NULL){
  load("Preprocess.RData")
  load("Clust_data.RData")
  gene_cpm = data$gene_cpm
  gene_cor = data$gene_cor
  MacroNutrition = data$MacroNutrition
  gene_entrz = data$gene_entrz
  clust = clust_data$clust
  clust_res = clust_data$clust_res
  clust_color = clust_data$clust_color
  gene_mean = apply(gene_cpm,1,mean)
  gene_med = apply(gene_cpm,1,median)
  gene_sd = apply(gene_cpm,1,sd)
  res = data.frame(Mean = gene_mean,Median = gene_med,SD = gene_sd,Cluster = clust_color)
  res = as.matrix(res)
  rownames(res) = rownames(gene_cpm)
  res = as.data.frame(res)
  return(res)
}

Clust_summary <- function(data,mclust1,mclust2,alpha,mcut,scaledist,thresh,k = NULL){
  
  load("Preprocess.RData")
  load("Clust_data.RData")
  
  gene_cpm = data$gene_cpm
  gene_cor = data$gene_cor
  MacroNutrition = data$MacroNutrition
  gene_entrz = data$gene_entrz
  clust = clust_data$clust
  clust_res = clust_data$clust_res
  clust_color = clust_data$clust_color
  clust_dist = clust_data$clust_dist
  cluster_vis = clust_vis(3,thresh,gene_cpm,clust_color,MacroNutrition)
  p1 = visNetwork(cluster_vis$nodes,cluster_vis$edges)%>% visHierarchicalLayout(levelSeparation = 500)
  
  return(p1)
}

visparameter <- function(topn,thresh,gene_cpm,clust_color,MacroNutrition){
  temp = table(clust_color)
  temp = temp[names(temp)!="grey"]
  gene_top = list()
  connection = list()
  for(i in 1:length(temp)){
    gene_group = t(gene_cpm)[,clust_color == labels2colors(i)]
    topn = min(topn,ncol(gene_group))
    group_pca = prcomp(gene_group,scale. = F)
    eigengene_group = group_pca$x[,1]

    eigen_cor1 = cor(gene_group,eigengene_group)
    gene_top[[i]] = eigen_cor1[order(abs(eigen_cor1[,1]),decreasing = T),][c(1:topn)]
    
    eigen_cor2 = cor(MacroNutrition,eigengene_group)
    connection[[i]] = as.vector(eigen_cor2)
    names(connection[[i]]) = colnames(MacroNutrition)
  }
  return(list(gene_top = gene_top,connection = connection))
}

clust_vis <- function(topn,thresh,gene_cpm,clust_color,MacroNutrition){
  temp = table(clust_color)
  temp = temp[names(temp)!="grey"]
  vis_parameter = visparameter(topn,thresh,gene_cpm,clust_color,MacroNutrition)
  gene_top = vis_parameter$gene_top
  connection = vis_parameter$connection
  nodes1 = data.frame(id = 1:(ncol(MacroNutrition) - 1),
                      level = 1,
                      label = colnames(MacroNutrition[,-1]),
                      title = colnames(MacroNutrition[,-1]),
                      color = "lightblue")
  nodes2 = data.frame(id = ncol(MacroNutrition):(ncol(MacroNutrition)+length(temp)-1),
                      level = 2,
                      label = "",
                      title = paste0("Module:",labels2colors(1:length(temp))),
                      color = labels2colors(1:length(temp))
  )
  edge12 = data.frame(from = 0, to = 0, color = "", title = 0)
  for (i in 1:length(connection))
  {
    for (j in 2:length(connection[[i]]))
      if (abs(connection[[i]][j])>thresh)
      {
        edge12 = rbind(edge12,
                       data.frame(from = j - 1,
                                  to =  ncol(MacroNutrition) + i - 1,
                                  color = ifelse(connection[[i]][j]>0,"red","blue"),
                                  title = connection[[i]][j])
        )
      }
  }
  edge12 = edge12[-1,]
  nodes3 = data.frame(id = 0,level = 3, label = "",title = "",color = "")
  edge23 = data.frame(from = 0, to = 0, color = "", title = 0)
  for (i in 1:length(gene_top))
  {
    for (j in 1:length(gene_top[[i]]))
    {
      id = ncol(MacroNutrition)+length(temp) + (i - 1)*topn + j - 1
      nodes3 = rbind(nodes3,data.frame(id = id,
                                       level = 3,
                                       label = names(gene_top[[i]][j]),
                                       title = names(gene_top[[i]][j]),
                                       color = labels2colors(i))
      )
      edge23 = rbind(edge23,data.frame(from = ncol(MacroNutrition) + i - 1,
                                       to = id,
                                       color = ifelse(gene_top[[i]][j]>0,"red","blue"),
                                       title = gene_top[[i]][j]))
    }
  }
  nodes3 = nodes3[-1,]
  edge23 = edge23[-1,]
  nodes = rbind(nodes1,nodes2,nodes3)
  edges = rbind(edge12,edge23)
  return(list(nodes = nodes,edges = edges))
}
