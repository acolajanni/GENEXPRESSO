packages.cran <- c("devtools","data.table","NanoStringNorm","pheatmap","nanostringr","ggrepel")
packages.bioC <- c("limma","GEOquery","DESeq2","vsn")

installed <- rownames(installed.packages())
packages.cran <- packages.cran[!packages.cran%in%installed]
packages.bioC <- packages.bioC[!packages.cran%in%installed]

if (length(packages.cran)>0){
  install.packages(packages.cran,dependencies=TRUE)
}

if (!requireNamespace("nanoR", quietly = TRUE)){
  install.packages("~/nanoR",repos=NULL,type="source")
}

if (!requireNamespace("NAPPA", quietly = TRUE)){
  install.packages("~/GIT/CPRD/NAPPA_2.0.tar.gz",repos=NULL,type="source")
}

if (length(packages.cran)>0){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(packages.bioC)
}

if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  library(devtools)
  install_github("vqv/ggbiplot")
}

if (!requireNamespace("UpSetR", quietly = TRUE)){
  install.packages("UpSetR")
}

if (!requireNamespace("reshape2", quietly = TRUE)){
  install.packages("reshape2")
}

if (!requireNamespace("WGCNA", quietly = TRUE)){
  BiocManager::install("WGCNA") 
}

if (!requireNamespace("RankProd", quietly = TRUE)){
  BiocManager::install("RankProd")
}
