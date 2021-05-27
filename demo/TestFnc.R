GSE13507 = tools.norm.Microarray(GEOiD = "GSE13507", FetchOnGEOdb = TRUE, tools = "none")
GSE32894  = tools.norm.Microarray(GEOiD = "GSE32894", FetchOnGEOdb = TRUE, tools = "none")
GSE32548= tools.norm.Microarray(GEOiD = "GSE32548", FetchOnGEOdb = TRUE, tools = "none")
GSE31684= tools.norm.Microarray(GEOiD = "GSE31684", FetchOnGEOdb = TRUE, tools = "none")


GSE_Test = tools.norm.Microarray(GEOiD = "./data/GSE31684/", FetchOnGEOdb = F, tools = "rma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")

#renv::install("C:/Users/Antonin\ COLAJANNI/Documents/GIT/GIT/org.Hs.eg.db_3.10.0.tar.gz", repos=NULL,type="source")
#renv::install("C:/Users/Antonin\ COLAJANNI/Documents/GIT/GIT/hgu133plus2.db", repos=NULL,type="source")
#renv::install("C:/Users/Antonin\ COLAJANNI/Documents/GIT/GIT/hgu133plus2cdf_2.18.0.tar.gz", repos=NULL,type="source")



BiocManager::install("hgu133plus2cdf",dependencies = TRUE,force=TRUE)

hgu133plus2cdf 
hgu133plus2cdf 