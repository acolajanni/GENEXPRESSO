###########################
#####  Normalization ######
###########################

library(nanostringr)
library(NanoStringNorm)
library(NAPPA)
library(nanoR)

#' Normalize a Nanostring dataset
#'
#' Normalization through different methods in several packages
#'
#' @param raw.data rcc type file. List of 4 elements: 
#' Samples.IDs is a dataframe with four columns: sample name, status (WildType, Mutated), filename. 
#' rcc.df contains expression levels with samples in column and genes in rows. 
#' annots.df is a dataframe with 3 columns: geneclass (endogenous, housekeeping, negative, positive), gene name, and the accession number. 
#' FOV is a dataframe that contains Field of views information.
#' 
#' @param tool Normalization tool. "nappa.NS", "nappa.param1","nappa.param2","nappa.param3" are different parameters used with the NAPPA() function from the NAPPA package.
#'  "nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2" use the NanoStringNorm() normalization function from the package NanoStringNorm
#'  "nanostringR" uses the HKnorm() function from the package nanostringr.
#'  "nanoR.top100","nanoR.total" uses the nsNormalize() function from the nanoR package.
#'  For the nanoR package, it is needed to give the file path to rcc files
#'
#' @param nanoR Logical value. TRUE is necessary if the tool used is "nanoR.top100" or "nanoR.total". By default, the value is set on FALSE
#' @param dir directory of rcc files. 
#' This parameter is only necessary if the nanoR normalizations are wanted.
#'
#' @return dataframe of normalized expression values
#' 
#' @import "nanostringr" "NanoStringNorm" "NAPPA" "nanoR"
#' @export
#'
#' @examples
#' # Retrieve Nanostring Data
#' # Data = Simul.data(type = "Nanostring")
#' # Normalize data using one method :
#' # Norm.data = tools.norm.Nanostring(raw.data = Data, tool = "nappa.NS",nanoR=F)
#' #with nanoR : Give the rcc files location
#' #RCC.dir <- file.path("./DATA/NANOSTRING","GSE146204_RAW")
#' #Norm.data = tools.norm.Nanostring(raw.data = Data, tool = "nanoR.top100",dir = RCC.dir,nanoR=T)
tools.norm.Nanostring <- function(raw.data,tool,nanoR=F,dir = NULL){
  # if the method is "nanoR.top100" or "nanoR.total", the function needs other informations contained in other files 
  if(!nanoR){
    rcc.samples <- raw.data$rcc.df
    annots.df <- raw.data$annots.df
    samples.IDs <- raw.data$samples.IDs
    FOV <- raw.data$FOV
    dir = NULL
  }
  
  if (missing(dir) & (nanoR = T)) {
    data.dir <- "./DATA/NANOSTRING"
    RCC.dir <- file.path(data.dir,"GSE146204_RAW")
    samples.IDs <- raw.data$samples.IDs
  }
  else if(nanoR){
    samples.IDs <- raw.data$samples.IDs
    raw.data = dir
  }
  
  tools.fnc <- switch(tool,
                      nanostringR={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        colnames(rcc.samples)[1:3] <- c("Code.Class", "Name", "Accession")
                        res.norm <- HKnorm(rcc.samples)
                        row.names(res.norm) <- res.norm$Name
                        res.norm <- res.norm[,-(1:3)]
                        colnames(res.norm) <- samples.IDs$ID
                        res.norm
                      },
                      nappa.NS={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        hk <- rcc.samples[grep("Housekeeping",rcc.samples$CodeClass),]
                        hk$CodeClass <- "Endogenous"
                        rcc.samples <- rbind(rcc.samples,hk)
                        res.norm <- NAPPA(rcc.samples,tissueType = "tumour",nposcontrols = 6,scaleFOV = F,output="All",raise.low.counts=1,poscontrol.method = "geometric.mean",background.method = "none",hk.method = "subtract")$GeneExpression
                        colnames(res.norm) <- samples.IDs$ID
                        res.norm
                      },  
                      nappa.default={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples <- rbind(FOV,rcc.samples)
                        hk <- rcc.samples[grep("Housekeeping",rcc.samples$CodeClass),]
                        hk$CodeClass <- "Endogenous"
                        rcc.samples <- rbind(rcc.samples,hk)
                        res.norm <- NAPPA(rcc.samples,tissueType = "tumour")
                        colnames(res.norm) <- samples.IDs$ID
                        res.norm
                        
                      },  
                      
                      nappa.param1={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples <- rbind(FOV,rcc.samples)
                        hk <- rcc.samples[grep("Housekeeping",rcc.samples$CodeClass),]
                        hk$CodeClass <- "Endogenous"
                        rcc.samples <- rbind(rcc.samples,hk)
                        res.norm <- NAPPA(rcc.samples,tissueType = "tumour",nposcontrols = 6,scaleFOV = TRUE,output="All",raise.low.counts=25)$GeneExpression
                        colnames(res.norm) <- samples.IDs$ID
                        res.norm
                      },
                      nappa.param3={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples <- rbind(FOV,rcc.samples)
                        hk <- rcc.samples[grep("Housekeeping",rcc.samples$CodeClass),]
                        hk$CodeClass <- "Endogenous"
                        rcc.samples <- rbind(rcc.samples,hk)
                        res.norm <- NAPPA(rcc.samples,tissueType = "tumour",nposcontrols = 6,scaleFOV = TRUE,output="All",raise.low.counts=1)$GeneExpression
                        colnames(res.norm) <- samples.IDs$ID
                        res.norm
                      },
                      nappa.param2={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples <- rbind(FOV,rcc.samples)
                        hk <- rcc.samples[grep("Housekeeping",rcc.samples$CodeClass),]
                        hk$CodeClass <- "Endogenous"
                        rcc.samples <- rbind(rcc.samples,hk)
                        res.norm <- NAPPA(rcc.samples,tissueType = "tumour",sampleNumber=10,hk.method="shrunken.subtract",nposcontrols = 6,scaleFOV = TRUE,output="All",raise.low.counts=25)$GeneExpression
                        colnames(res.norm) <- samples.IDs$ID
                      },
                      nanostringnorm.default={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples$CodeClass[rcc.samples$CodeClass=="Housekeeping"] <- "Endogenous"
                        res.norm <- NanoStringNorm(rcc.samples,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                      },
                      
                      nanostringnorm.param1={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples$CodeClass[rcc.samples$CodeClass=="Housekeeping"] <- "Endogenous"
                        res.norm <- NanoStringNorm(rcc.samples,CodeCount = 'geo.mean',Background = 'mean.2sd',SampleContent = 'low.cv.geo.mean',round.values = TRUE,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                      },
                      nanostringnorm.param2={
                        design <- model.matrix(~0+samples.IDs$tp53.status)
                        colnames(design) <- c("Mutated","WildType")
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples$CodeClass[rcc.samples$CodeClass=="Housekeeping"] <- "Endogenous"
                        res.norm <- NanoStringNorm(rcc.samples,Background = 'mean.2sd',SampleContent = 'low.cv.geo.mean',CodeCount="sum",traits=design,round.values = TRUE,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                        
                      },
                      nanoR.top100={
                        nano <- parseRCC(dir = raw.data)
                        nano <- nsBackgroundCorrect(nano)
                        nano <- nsPositiveControlNormalization(nano)
                        res.norm <- nsNormalize(nano, method="top100")$bg.corr.counts
                        res.norm <- log2(res.norm[!grepl("Neg|Pos",res.norm$CodeClass),-c(1:3)]+1)
                        colnames(res.norm) <- samples.IDs$title
                        res.norm
                        
                      },
                      nanoR.total={
                        nano <- parseRCC(dir = raw.data)
                        nano <- nsBackgroundCorrect(nano)
                        nano <- nsPositiveControlNormalization(nano)
                        res.norm <- nsNormalize(nano, method="total")$bg.corr.counts
                        res.norm <- log2(res.norm[!grepl("Neg|Pos",res.norm$CodeClass),-c(1:3)]+1)
                        colnames(res.norm) <- samples.IDs$title
                        res.norm
                      },
                      stop("Enter something that switches me!")          
  )
  return(res.norm)
}

#' Normalize a Nanostring dataset
#'
#' @param GEOiD GEO accession number or directory. 
#' If you give a GEO accession number, then set FetchOnGEOdb = TRUE.
#' Otherwise, if it is a directory, ignore the FetchOnGEOdb parameter. Note that if the data are stored localy, it should be inside a .tar archive.
#' @param FetchOnGEOdb logical value. if TRUE, GEOiD parameter should be a GEO accession number.
#' @param tools.normalize 
#' @param tools.bgcorrect 
#' @param tools.pmcorrect 
#' @param tools.express.summary.stat 
#' @param tools 
#' @import "affy" "GEOquery" "gcrma"
#' @return
#' @export
#'
#' @examples
tools.norm.Microarray <-function(GEOiD , FetchOnGEOdb = FALSE , tools, tools.normalize, tools.bgcorrect, tools.pmcorrect, tools.express.summary.stat){
  
  #########################
  ###  Work in progress ###
  #########################
  
  
  
  if (!FetchOnGEOdb){
    # If we already have an "AffyBatch", we can skip thoses steps
    if (class(GEOiD) != "AffyBatch") {
      # Getting the directory of .Cel files
      celpath = celpath = file.path(GEOiD)
      files <- list.files(path = celpath, pattern = "CEL.gz", full.names = TRUE)
      abatch <- ReadAffy(filenames = files)
    }  
    else {
      message("Going straight to normalization step")
    }
  }
  else {
    # Download the archive
    getGEOSuppFiles(GEOiD, fetch_files = TRUE, baseDir = "./data")  
    # get the directory
    celpath = paste0("./data/",GEOiD,"/")
    # extracting it
    tarfile = paste0(celpath, GEOiD,"_RAW.tar")
    # Extract the .cel files from the archive
    untar(tarfile = tarfile, exdir = celpath)
    # Retrieve the filenames
    files <- list.files(path = celpath, pattern = "CEL.gz", full.names = TRUE)
    # Construct the AffyBatch object
    abatch <- ReadAffy(filenames = files)
  }
  message("Switch")
  # Depending on the applied methods, different functions will be called
  tools.fnc = switch(tools,
                     # rma, gcrma and mas5 are pretty straigthforward functions
                     rma = {
                       eset = rma(abatch)
                       return(eset)
                     },
                     
                     gcrma = {
                       eset = gcrma(abatch)
                     },
                     
                     mas5 = {
                       eset = mas5(abatch)
                     }, 
                     
                     none = {
                       return(abatch)
                     },
                     
                     # Custom normalization through different method : 
                     # Normalization, background correction, pm correction, summary stat expression
                     custom = {
                       
                       # Error management for each tool method given in argument
                       bgcorrect.methods = bgcorrect.methods()
                       bgcorrect.methods = bgcorrect.methods[bgcorrect.methods == "rma" | bgcorrect.methods == "mas" ]
                       if (!tools.bgcorrect%in%bgcorrect.methods){
                         stop("Enter a background correction method that switches me!")
                       }
                       
                       
                       norm.methods = normalize.AffyBatch.methods()
                       norm.methods = norm.methods [norm.methods != "methods" & norm.methods != "vsn" ]
                       if (!tools.normalize%in%norm.methods){
                         stop("Enter a normalization method that switches me!")
                       }
                       
                       pmcorrect.methods = pmcorrect.methods()
                       pmcorrect.methods = pmcorrect.methods[pmcorrect.methods != "methods"]
                       if (!tools.pmcorrect%in%pmcorrect.methods){
                         stop("Enter a pm correction method that switches me!")
                       }
                       
                       sum.stat = express.summary.stat.methods()
                       if (!tools.express.summary.stat%in% sum.stat){
                         stop("Enter a summary stat method that switches me!")
                       }
                       
                       # If all tools.normalize, tools.bgcorrect, .... are correct
                       # We can call the expresso() function with the four tools parameters : 
                       eset <- expresso(abatch, 
                                        bgcorrect.method= tools.bgcorrect,
                                        normalize.method= tools.normalize,
                                        pmcorrect.method= tools.pmcorrect,
                                        summary.method= tools.express.summary.stat)
                       
                       }
                     )
  
  exprSet = exprs(eset)
  return(exprSet)
}


#' Mapping Probes ID of Affymetrix microarray chip in expression set dataframe
#'
#' @param exprSet Dataframe with samples in columns, probes ID in rows
#' @param annotation function like : hgu133plus2SYMBOL if the chip is hgu 133 plus 2.0 and the wanted mapping is with the gene SYmbols
#' 
#' @import "affy" "dplyr"
#' 
#' @return Dataframe with about 20.000 rows. Gene symbols will replace the probes ID
#' @export
#'
#' @examples
mapping.affymetrix.probe <- function(exprSet){
  
  #########################
  ###  Work in progress ###
  #########################
  
  
  # Remove control probes from the expression set
  ControlProbes <- grep("AFFX",row.names(exprSet)) 
  
  if (length(ControlProbes) != 0){
    expr = exprSet[ - ControlProbes,]
  }
  else {
    expr = exprSet
  }
  
  # Retrieve probe names
  probes.ALL=row.names(expr)
  # We want the annotation through the gene Symbol : 
  symbol.ALL = unlist(mget(probes.ALL, hgu133plus2SYMBOL))
  
  # Recreating the dataframe with the matching probes
  table.ALL=cbind(SYMBOL = symbol.ALL,  expr)
  table.ALL$PROBES = row.names(expr)
  table.ALL = relocate(table.ALL, PROBES, SYMBOL)
  row.names(table.ALL) = NULL

  expr = table.ALL
  # Retrieving sample names 
  samples = colnames(expr)
  samples = samples[!samples %in% c("PROBES","SYMBOL")]
  
  # grouping expr by the Gene symbols
  tmp = expr %>%
    group_by(SYMBOL) 
  
  # Summarised the group with the median 
  Mapped = tmp %>%
    summarise(across(all_of(samples), ~ median(.x)  ))
  
  Mapped = as.data.frame( na.omit(Mapped) )
  row.names(Mapped) = Mapped$SYMBOL
  Mapped$SYMBOL = NULL
  
  
  return(Mapped)
}
