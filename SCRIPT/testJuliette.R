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

nappaNS <- tools.inspect(raw.data, "nappa.NS", nanoR = F)
nappaNS <- as.data.frame(nappaNS)

nanostringR <- tools.inspect(raw.data, "nanostringR", nanoR = F)

nappa.default <- tools.inspect(raw.data, "nappa.default", nanoR = F)
nappa.default <- as.data.frame(nappa.default)

test <- read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)
boucle1 <- wilcoxDEG(test, 42,22)
bucle3 <- wilcoxDEGoption2(test,42,22)


wilcox <- data.frame()
for (i in 1:nrow(test)){
  x = test[i,1:42]
  y = test[i,43:64]
  test1 <- wilcox.test(as.numeric(x),as.numeric(y))
  test2 <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "less")
  test3 <- wilcox.test(as.numeric(x),as.numeric(y), alternative = "greater")
  pval1 <- test1$p.value
  pval2 <- test2$p.value
  pval3 <- test3$p.value
  pval <- c(pval1,pval2,pval3)
  pval <- as.data.frame(pval)
  pvat(pval)
  wilcox <- rbind(wilcox,pval)
}
row.names(wilcox) = row.names(data)
return(wilcox)
