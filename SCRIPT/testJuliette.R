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

#nappaNS <- tools.inspect(raw.data, "nappa.NS", nanoR = F)
#nappaNS <- as.data.frame(nappaNS)

tools <- c("nanostringR", "nappa.NS", "nappa.default", "nappa.param1","nappa.param3", "nappa.param2", "nanostringnorm.default", "nanostringnorm.param1", "nanostringnorm.param2", "nanoR.top100", "nanoR.total")
for (i in tools){
  
}
