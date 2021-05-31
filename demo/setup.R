packages.cran <- c("madsim","factoextra","FactoMineR","UpSetR","igraph","reshape2","devtools","data.table","NanoStringNorm","pheatmap","nanostringr","ggrepel")
packages.bioC <- c("limma","GEOquery","DESeq2","DESeq",
                   "vsn","WGCNA","RankProd","edgeR",
                   "DEFormats","gcrma","hgu133plus2.db",
                   "biomaRt", "affy", "affyio", "impute",
                   "GO.db", "AnnotationDbi", "Biobase"
                   )
# Indisponible sur R 4.0 :
#packages = c("DESeq2","DESeq","RankProd","edgeR", "DEFormats", "WGCNA")
#pkgs = c("DEFormats", "madsim", "edgeR", "DESeq", "DESeq2", "RankProd", "limma", "nanostringr", "NanoStringNorm", "reshape2", "WGCNA", "factoextra", "FactoMineR", "UpSetR", "data.table", "dplyr")
#problème : vsn doit être installé avant NanostringNorm
installed <- rownames(installed.packages())
packages.cran <- packages.cran[!packages.cran%in%installed]
packages.bioC <- packages.bioC[!packages.bioC%in%installed]

if (length(packages.cran)>0){
  install.packages(packages.cran,dependencies=TRUE)
}

if (!requireNamespace("nanoR", quietly = TRUE)){
  install.packages("~/GIT/GENEXPRESSO/nanoR", repos=NULL,type="source")
  
  
  #renv::install("~/GIT/GENEXPRESSO/nanoR", repos=NULL,type="source")
  #renv::install("C:/Users/Antonin\ COLAJANNI/Documents/GIT/GIT/GENEXPRESSO/nanoR", repos=NULL,type="source")
}
if (!requireNamespace("NAPPA", quietly = TRUE)){
  install.packages("~/GIT/GENEXPRESSO/NAPPA_2.0.tar.gz", repos = NULL, type = "source")
  #renv::install("~/GIT/GENEXPRESSO/NAPPA_2.0.tar.gz", repos = NULL, type = "source")
}

if (length(packages.bioC)>0){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (package in packages.bioC){
    BiocManager::install(package, dependancies = TRUE)
  }
}

if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  library(devtools)
  install_github("vqv/ggbiplot")
}







#if (!requireNamespace("UpSetR", quietly = TRUE)){
#  install.packages("UpSetR")
#}

#if (!requireNamespace("reshape2", quietly = TRUE)){
#  install.packages("reshape2")
#}

#if (!requireNamespace("WGCNA", quietly = TRUE)){
#  BiocManager::install("WGCNA") 
#}

if (!requireNamespace("RankProd", quietly = TRUE)){
  BiocManager::install("RankProd")
}
#if (!requireNamespace("igraph", quietly = TRUE)){
#  BiocManager::install("igraph")
#}

if (!requireNamespace("factoextra", quietly = TRUE)){
  BiocManager::install("factoextra")
}
if (!requireNamespace("FactoMineR", quietly = TRUE)){
  BiocManager::install("FactoMineR")
}
if (!requireNamespace("GEOquery", quietly = TRUE)){
  BiocManager::install("GEOquery")
}
if (!requireNamespace("edgeR", quietly = TRUE)){
  BiocManager::install("edgeR")
}
if (!requireNamespace("DEFormats", quietly = TRUE)){
  BiocManager::install("DEFormats")
}
if (!requireNamespace("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
}
if (!requireNamespace("DESeq", quietly = TRUE)){
  BiocManager::install("DESeq")
}
if (!requireNamespace("madsim", quietly = TRUE)){
  BiocManager::install("madsim")
}
