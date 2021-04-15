# Nom               : make_comparison_table
# Type              : Programme
# Objet             : 
# Input             : Fichier de données nanostrings
# Output            : data.to.comp.csv
# Auteur            : EDarbo
# R version         : 3.6
# Date de creation  : 04.03.2021
#______________________________________________________________________________

# go to the directory you downloaded
setwd("~/GIT/CPRD")
# define path to usefull directories
data.dir <- "./DATA/NANOSTRING"
RCC.dir <- file.path(data.dir,"GSE146204_RAW")



# install packages
source(file.path("./SCRIPT","setup.R"))

# source functions
source(file.path("./SCRIPT","functions.R"))

# use GEOQuery to upload annotation data
GOIDs <- c("GSE146204")
GEOdata <- getGEO(GOIDs,GSEMatrix=T, AnnotGPL=F)
sample.annots <- get.annotations(GEOdata)
sample.annots <- sample.annots[,c("title","tp53 status")]

# get the names of the RCC files containing the raw data
rcc.files <- system(paste("ls",RCC.dir),intern=T)
# map file names with sample IDs
samp.GSM <- unlist(sapply(strsplit(rcc.files,"_"),"[[",1))
m <- match(samp.GSM,row.names(sample.annots))
sample.annots$fileName <- sub(".RCC","",rcc.files)
colnames(sample.annots) <- sub(" ",".",colnames(sample.annots))
sample.annots$tp53.status <- sub("Wild type/NA/ND","WildType",sample.annots$tp53.status)
# read RCC files from your computer
raw.data <- import.raw.rcc(RCC.dir,sample.annots)
samples.IDs <- raw.data$samples.IDs

# options for tools :
# nappa.NS, nappa.param1, nappa.param2,nanostringnorm,desq2,nanoR.top100,nanostringR,nanoR.total

# compute all comparisons: each tool will normalize the data and compute differentation expression
# The function returns the data.frame with 2 columns : the log Fold Change and the SYMBOL of the gene
tools <- c("nappa.NS","nappa.param1", "nappa.param2","nappa.param3","nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2","desq2","nanoR.top100","nanoR.total","nanostringR")

data.to.comp <- tools.inspect(raw.data,tool="nappa.default",nanoR=F)
data.to.comp

for (tool in tools){
  print(tool)
  nanoR=F
  raw <- raw.data
  if (tool%in%c("nanoR.top100","nanoR.total")){
    nanoR <- T
    raw <- RCC.dir
  }
  tmp <- tools.inspect(raw,tool,nanoR)
  data.to.comp <- merge(data.to.comp,tmp,by="SYMBOL",all=T)
}

row.names(data.to.comp) <- data.to.comp$SYMBOL
data.to.comp <- data.to.comp[,-1]

# if you need to remove genes not present in all analyses: NA (missing) data
data.to.comp <- na.omit(data.to.comp)

# transpose the matrix to use the genes as descriptive values to compare the tools
data.to.comp <- as.data.frame(t(data.to.comp))
data.to.comp
#Export as csv
write.csv(data.to.comp,"~/GIT/CPRD/OUTPUT/TABLE/data.to.comp.csv", row.names = TRUE)

################################################################################
################################################################################
# compute PCA
library("factoextra")
library("FactoMineR")
res.pca <- PCA(data.to.comp)
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 90))
var <- get_pca_var(res.pca)
fviz_pca_var(res.pca, col.var = "cos2")
fviz_pca_ind (res.pca, col.ind = "cos2",
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

################################################################################
################################################################################
# compute intersection with upsetR

Upset <- copy(data.to.comp)
for (i in names(Upset)){
  for (u in 1:nrow(Upset)){
    if(data.to.comp[[i]][u] > 0.05){
      Upset[[i]][u] = 0
    }else{
      Upset[[i]][u] = 1
    }
  }
}
Upset = as.data.frame(t(Upset))

# Avec les non-intersections
upset(Upset, sets = names(Upset), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on" )
# Sans les non-intersections
upset(Upset, sets = names(Upset), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = NULL )

