preprocessing <- function(GE,N,mcpm,maxcpm,mvar,mabscor,mnorm){
  source("utils.R")
  
  if (is.null(GE)|is.null(N)){
    load("data.RData")
    gene_cpm = data$gene_cpm
    MacroNutrition = data$MacroNutrition
  }
  tryCatch(
    {
      gene_cpm <- read.csv(file = GE)
      MacroNutrition <- read.csv(file = N)
    },
    error = function(e) {
      load("data.RData")
      gene_cpm = data$gene_cpm
      MacroNutrition = data$MacroNutrition
    }
  )
  

  gene_cpm = na.omit(gene_cpm)
  MacroNutrition = na.omit(MacroNutrition)
  A = rowSums(gene_cpm)
  gene_cpm = gene_cpm[A>0,]
  
  if (mnorm == "cpm"){
    gene_cpm = gene_cpm/rowSums(gene_cpm) * 1000000
  }else if(mnorm == "lcpm"){
    gene_cpm = gene_cpm/rowSums(gene_cpm) * 1000000
    gene_cpm = log2(gene_cpm + 1/2)
  }else if(mnorm == "none"){
    gene_cpm = gene_cpm
  }
  

  A = A[A>0]
  B = apply(gene_cpm,1,sd)
  C = melt(log2(A),variable.names = "Sample")
  D = melt(log2(B),variable.names = "Sample")
  
  p1 = ggplot(C, aes(x = value)) +
    geom_density(fill = "lightblue",alpha = 0.5, size = 1)  + 
    labs(x = expression(log[2](CPM)),title = paste0("density plot of gene expression(log2(cpm))")) + 
    geom_vline(xintercept = log2(mcpm),color = "red",size = 1.25) + 
    geom_vline(xintercept = log2(maxcpm),color = "red",size = 1.25) + theme_classic()
  
  p2 = ggplot(D, aes(x = value)) + 
    geom_density(fill = "green",alpha = 0.5, size = 1)  + 
    labs(x = expression(log[2](SD_CPM)),title = paste0("density plot of gene expression standard deviation(log2(sd_cpm))")) + 
    geom_vline(xintercept = log2(mvar),color = "red",size = 1.25) + theme_classic()
  
  gene_cpm = gene_cpm[A>mcpm&B>mvar&A<maxcpm,]

  
  
  gene_cor = cor(MacroNutrition,t(gene_cpm),method = "spearman")
  gene_cor_max = apply(abs(gene_cor),2,max)
  E = melt(gene_cor_max,variable.names = "Sample")
  
  p3 = ggplot(E, aes(x = value)) + 
    geom_density(fill = "yellow",alpha = 0.5, size = 1.25) + 
    labs(x = "max(abs(corr))",title = paste0("density plot of gene expression(cpm) maxium absolute correlation with diet")) + 
    geom_vline(xintercept = mabscor,color = "red",size = 1.25) + theme_classic()
  
  gene_cpm = gene_cpm[gene_cor_max>mabscor,]
  
  gene_cor = gene_cor[,gene_cor_max>mabscor]
  
  data = list(gene_cpm = gene_cpm,gene_cor = gene_cor,MacroNutrition = MacroNutrition)
  save(data,file = "Preprocess.RData")
  return(multiplot(p1,p2,p3,cols = 2))
}

numgene <- function()
{
  load("Preprocess.RData")
  gene_cpm = data$gene_cpm
  gene_cor = data$gene_cor
  MacroNutrition = data$MacroNutrition
  gene_entrz = data$gene_entrz
  n = nrow(gene_cpm)
  return(n)
}