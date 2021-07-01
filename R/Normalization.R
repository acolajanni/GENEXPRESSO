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



#' Compute normalization / size factors for RNAseq count matrix
#' 
#' This function calls various function to compute normalization, size factors, or normalized dataset 
#' to call another function that analyses differentially expressed genes
#' 
#' @param count.matrix 
#' Dataframe of count with samples in columns and genes SYMBOL in rows.
#' 
#' @param tool 
#' Character string among "TMM","TMMwsp", "RLE", "Upperquartile", "voom", "vst", "vst2".
#' 
#' "TMM","TMMwsp", "RLE", "Upperquartile" calls the \link{edgeR}{calcNormFactors} function.
#' "voom" calls the \link{limma}{voom} function.
#' "vst" calls the \link{DESeq}{estimateSizeFactors} function on a CountDataSet.
#' "vst2" does the same but also calls the \link{DESeq2}{varianceStabilizingTransformation} function.
#' 
#' @param design 
#' Vector of 1 and 2 of the same length of colnames(count.matrix).
#' 1 for the first group and 2 for the second.
#'
#' @import "DEFormats" "edgeR" "DESeq2" "DESeq" "limma"
#'
#' @return 
#' "TMM","TMMwsp", "RLE", "Upperquartile" and "vst" returns a vector of the same size as colnames(count.matrix)
#' "voom" returns an "Elist" class object.
#' "vst2" returns the normalized count.matrix with the variance stabilizing transformation
#' 
#' @export
#'
#' @examples
#' # load a count matrix (example with a random dataset)
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute design vector
#' design = c(rep(1,10), rep(2,10)) # 10 from group 1, 10 from group 2
#' 
#' Norm = tools.norm.RNAseq(Data, "TMM", design)
#' 
tools.norm.RNAseq <- function(count.matrix, tool, design){
  # EdgeR normalisation need a DGElist :
  if (tool%in%c("TMM", "TMMwsp", "RLE", "Upperquartile", "voom")){
    edgeR.dgelist = DGEList(counts = count.matrix, group = factor(design))
  }
  storage.mode(count.matrix) = "integer"
  
  
  tools_norm_RNAseq.fnc <- switch(tool,
                                  ##################################### edgeR
                                  TMM = {
                                    nf = calcNormFactors(edgeR.dgelist, method = "TMM")
                                  },
                                  
                                  TMMwsp = {
                                    nf = calcNormFactors(edgeR.dgelist, method = "TMMwsp")
                                  },
                                  
                                  RLE = {
                                    nf = calcNormFactors(edgeR.dgelist, method = "RLE")
                                  },
                                  
                                  Upperquartile = {
                                    nf = calcNormFactors(edgeR.dgelist, method = "upperquartile")
                                  },
                                  #####################################
                                  
                                  ##################################### DESeq
                                  vst = {
                                    # CountDataSet object
                                    DESeq.cds = newCountDataSet(countData = count.matrix,
                                                                conditions = factor(design))
                                    
                                    # Computing size factors (normalization factors)
                                    DESeq.cds = estimateSizeFactors(DESeq.cds)
                                    nf = sizeFactors(DESeq.cds)
                                    
                                  },
                                  ##################################### DESeq2
                                  vst2 = {
                                    design = data.frame(design,row.names=colnames(count.matrix))
                                    design$design = as.factor(design$design)
                                    # DESeqDataSet (deseq2) object
                                    dds<-DESeqDataSetFromMatrix(count.matrix,
                                                                colData = design,
                                                                design= ~design)
                                    
                                    # Size factors + dispersions
                                    dds = estimateSizeFactors(dds)
                                    dds = estimateDispersions(dds)
                                    
                                    # Returning normalized data
                                    DESeq.vst = getVarianceStabilizedData(dds)
                                    # Returns a matrix
                                    return(DESeq.vst)
                                    
                                  },
                                  ##################################### limma
                                  voom = {
                                    #limma - voom needs a DGElist to normalize
                                    nf = calcNormFactors(count.matrix, method = "TMM")
                                    voom.data = voom(count.matrix, 
                                                     design = model.matrix(~factor(design)),
                                                     lib.size = colSums(count.matrix) * nf)
                                    
                                    # Returns an "Elist" (voom format)
                                    return(voom.data)
                                  },
                                  stop("Enter a normalization method that switches me !")     )
  
  
  if (tool%in%c("TMM", "TMMwsp", "RLE", "Upperquartile")){
    # adding sample names to normalization factors 
    nf = nf[["samples"]][["norm.factors"]]
    names(nf) = colnames(count.matrix)
  } 
  return(nf)
  
}


#' Normalize a Nanostring dataset
#' 
#' Import a dataset in GEOdatabase and normalize it with several methods
#' 
#' @param GEOiD GEO accession number or directory. 
#' If you give a GEO accession number, then set FetchOnGEOdb = TRUE.
#' Otherwise, if it is a directory, ignore the FetchOnGEOdb parameter. 
#' Note that if the data are stored localy, it should be inside a .tar archive.
#' 
#' @param FetchOnGEOdb 
#' logical value. if TRUE, GEOiD parameter should be a GEO accession number.
#' 
#' @param tools 
#' Character string among : 
#' "rma", "gcrma", "mas5","none", "liwong", "RMA.inv.mas", "mas.mas.inv.mas", 
#' "mas.mas.const.liwong", "mas.mas.inv.med", "mas.mas.inv.liwong"
#' "rma" calls \link[affy]{rma}.
#' "mas5" calls \link[affy]{mas5}.
#' "gcrma" calls \link[gcrma]{gcrma}.
#' The other methods calls \link[affy]{expresso}.
#' "liwong" has the following parameters : (bgcorrect.method= "none", 
#' normalize.method= "invariantset", pmcorrect.method= "pmonly", summary.method= "liwong")
#' Methods that start with "RMA." has the parameters : 
#' bgcorrect.method = "rma" , pmcorrect.method = "pmonly" 
#' Those that starts with "mas.mas" have the following parameters : 
#' bgcorrect.method = "rma" , pmcorrect.method = "pmonly" 
#' ".inv" is for  normalize.method = "invariantset", 
#' ".const" is for normalize.method = "constant" 
#' ".liwong" is for summary.method= "liwong"
#' ".med" is for summary.method= "medianpolish"
#' ".mas" is for summary.method= "mas"
#' 
#' 
#' @import "affy" "GEOquery" "gcrma"
#' @return
#' "none" returns an "affybatch" object
#' Other methods returns a matrix of normalized data : Probe ID in rows sample in columns.
#' 
#' @export
#'
#' @examples
#' # giving a GSE ID
#' GEO = "GSE31684"
#' 
#' # Need to fetch those data on GEO: (returns an affybatch)
#' # Abatch = tools.norm.Microarray(GEOiD = GEO, FetchOnGEOdb = TRUE, tools = "none")
#' 
#' # Normalizing it : 
#' # norm = tools.norm.Microarray(GEOiD = Abatch, FetchOnGEOdb = FALSE, tools = "rma" ) 
tools.norm.Microarray <-function(GEOiD , FetchOnGEOdb , tools){

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
    message("Downloading data")
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
                     
                     liwong = {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "none",
                                        normalize.method= "invariantset",
                                        pmcorrect.method= "pmonly",
                                        summary.method= "liwong")
                     },

                     RMA.inv.mas = {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "rma",
                                        normalize.method= "invariantset",
                                        pmcorrect.method= "pmonly",
                                        summary.method= "mas")
                       
                       
                     },
                     mas.mas.inv.mas= {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "mas",
                                        normalize.method= "invariantset",
                                        pmcorrect.method= "mas",
                                        summary.method= "mas")
                     },

                     mas.mas.const.liwong = {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "mas",
                                        normalize.method= "constant",
                                        pmcorrect.method= "mas",
                                        summary.method= "liwong")
                     },
                     

                     mas.mas.inv.med = {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "mas",
                                        normalize.method= "invariantset",
                                        pmcorrect.method= "mas",
                                        summary.method= "medianpolish")
                     },

                     mas.mas.inv.liwong = {
                       eset <- expresso(abatch, 
                                        bgcorrect.method= "mas",
                                        normalize.method= "invariantset",
                                        pmcorrect.method= "mas",
                                        summary.method= "liwong")
                     },

                     stop("Enter something that switches me!")
                     
                     )
  # Returning an "expression set" = Matrix of expression values
  exprSet = exprs(eset)
  return(exprSet)
}


#' Mapping Probes ID of Affymetrix microarray chip in expression set dataframe
#'
#' @param exprSet Dataframe with samples in columns, probes ID in rows
#' @param annotation function like : hgu133plus2SYMBOL if the chip is hgu 133 plus 2.0 and the wanted mapping is with the gene SYmbols
#' 
#' @import "affy" "dplyr" "hgu133plus2.db"
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
  symbol.ALL = unlist(mget(probes.ALL, hgu133plus2.db::hgu133plus2SYMBOL))
  # Recreating the dataframe with the matching probes
  SYMBOL = unname(symbol.ALL)
  expr = cbind(SYMBOL,expr)
  
  expr <- as.data.frame(apply(expr, 2, as.numeric))
  expr$SYMBOL = SYMBOL
  
  
  
  # Retrieving sample names 
  samples = colnames(expr)
  samples = samples[!samples %in% c("PROBES","SYMBOL")]
  

  # grouping expr by the Gene symbols
  tmp = expr %>%
    group_by(SYMBOL)
    

  # Summarised the group with the median 
  Mapped = tmp %>%
    summarise(across(all_of(samples), ~ median(.x , na.rm = TRUE)  ))
  

  Mapped = as.data.frame( na.omit(Mapped) )
  row.names(Mapped) = Mapped$SYMBOL
  Mapped$SYMBOL = NULL
  
  return(Mapped)
}
