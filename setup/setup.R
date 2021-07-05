# Define the package directory (ends with "/GENEXPRESSO")
setwd("/Path/to/GENEXPRESSO/")


# List of necessary packages
# From CRAN :
packages.cran <- c("madsim","factoextra","FactoMineR","UpSetR","igraph",
                   "reshape2","devtools","data.table","NanoStringNorm",
                   "pheatmap","nanostringr","ggrepel", "viridis")

# From Bioconductor
packages.bioC <- c("limma","GEOquery","DESeq2","DESeq",
                   "vsn","WGCNA","RankProd","edgeR",
                   "DEFormats","gcrma","hgu133plus2.db",
                   "biomaRt", "affy", "affyio", "impute",
                   "GO.db", "AnnotationDbi", "Biobase",
                   "hgu133plus2.db"
                   )

# List of needed package for GENEXPRESSO :  
installed <- rownames(installed.packages())
packages.cran <- packages.cran[!packages.cran%in%installed]
packages.bioC <- packages.bioC[!packages.bioC%in%installed]


# Two packages are needed (NAPPA and nanoR), and available directly in github : acolajanni/GENEXPRESSO 
if (!requireNamespace("nanoR", quietly = TRUE)){
  install.packages("./nanoR", repos=NULL,type="source")
}
if (!requireNamespace("NAPPA", quietly = TRUE)){
  install.packages("./NAPPA_2.0.tar.gz", repos = NULL, type = "source")
}



# Installation of Bioconductor packages :
if (length(packages.bioC)>0){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  for (package in packages.bioC){
    BiocManager::install(package, dependancies = TRUE)
  }
}

# Installation of cran packages :
if (length(packages.cran)>0){
  install.packages(packages.cran,dependencies=TRUE)
}

# ggbiplot from github is necessary too
if (!requireNamespace("ggbiplot", quietly = TRUE)) {
  library(devtools)
  install_github("vqv/ggbiplot")
}

# Once all the previous steps are comple, GENEXPRESSO should install properly
if (!requireNamespace("GENEXPRESSO", quietly = TRUE)) {
  library(devtools)
  install_github("acolajanni/GENEXPRESSO")
}