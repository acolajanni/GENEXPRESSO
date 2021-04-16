suppressWarnings(library(data.table))
suppressWarnings(library(NAPPA))
suppressWarnings(library(NanoStringNorm))
suppressWarnings(library(nanostringr))
suppressWarnings(library(nanoR))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggbiplot))
suppressWarnings(library(limma))
suppressWarnings(library(GEOquery))
suppressWarnings(library(ggbiplot))
library(ggrepel)
library(DESeq2)
source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))
library("factoextra")
library("FactoMineR")
library(edgeR)
library(DEFormats)
library(DESeq)
library(UpSetR)

get.annotations <- function(x){
  annots <- phenoData(x[[1]])@data
  characteristics_ch1 <- annots[,grepl("characteristics_ch1",colnames(annots))]
  if (class(characteristics_ch1)=="factor"){
    print("Analysing Factor")
    characteristics_ch1 <- as.vector(characteristics_ch1)
    characteristics_ch1 <- strsplit(characteristics_ch1,";")
    OS <- unlist(sapply(characteristics_ch1,"[[",2))
    STATUS <- unlist(sapply(characteristics_ch1,"[[",3))
    OS <- as.numeric(as.vector(unlist(sapply(strsplit(OS,": "),"[[",2))))
    OS_IND <- as.numeric(as.vector(unlist(sapply(strsplit(STATUS,": "),"[[",2))))
    characteristics <- data.frame(OS,OS_IND)
  }
  else {
    print("Analysing Matrix")
    characteristics <- annots[,grepl(":ch1",colnames(annots))]
    colnames(characteristics) <- gsub(":ch1","",colnames(characteristics))
    print(head(characteristics))
    characteristics <- as.data.frame(characteristics)
  }
  annots <- annots[,grepl("data_processing|platform_id|description|scan_protocol|title|contact_name",colnames(annots))]
  annots$GSMID <- row.names(annots)
  annots <- cbind(annots,characteristics)
}

import.raw.rcc <- function(data.dir,annots){
  rcc <- read.markup.RCC(rcc.path = data.dir,rcc.pattern = "*.RCC|*.rcc",include = NULL,nprobes = -1)
  rcc.df <- rcc[[1]]
  gsm <- unlist(sapply(strsplit(colnames(rcc.df)[-c(1:3)],"_"),"[[",1))
  m <- match(row.names(annots),gsm)
  colnames(rcc.df)[-c(1:3)] <- as.vector(annots$title[m])
  annots.df <- rcc.df[,1:3]
  rcc.df <- rcc.df[,-c(1:3)]
  row.names(rcc.df) <- annots.df$Name
  temp <- rcc[[2]]
  temp <- data.frame(CodeClass=row.names(temp),Name="",Accession="",temp)
  temp <- rbind(temp,temp[c("FovCount","FovCounted"),])
  temp$CodeClass <- as.vector(temp$CodeClass)
  temp[c(nrow(temp)-1,nrow(temp)),"CodeClass"] <- c("FOV Count","FOV Counted")
  colnames(temp)[-c(1:3)] <- colnames(rcc.df)
  
  #rcc.df <- rbind(temp,rcc.df)
  return(list(samples.IDs=annots,rcc.df=rcc.df,annots.df=annots.df,FOV=temp))
}

tools.inspect <- function(raw.data,tool,nanoR=F){
  if (!nanoR){
    rcc.samples <- raw.data$rcc.df
    annots.df <- raw.data$annots.df
    samples.IDs <- raw.data$samples.IDs
    FOV <- raw.data$FOV
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
                        #rcc.samples$CodeClass[rcc.samples$Name%in%test.samp] <- "Housekeeping"
                        res.norm <- NanoStringNorm(rcc.samples,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                      },
                      
                      nanostringnorm.param1={
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples$CodeClass[rcc.samples$CodeClass=="Housekeeping"] <- "Endogenous"
                        #rcc.samples$CodeClass[rcc.samples$Name%in%test.samp] <- "Housekeeping"
                        res.norm <- NanoStringNorm(rcc.samples,CodeCount = 'geo.mean',Background = 'mean.2sd',SampleContent = 'low.cv.geo.mean',round.values = TRUE,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                      },
                      nanostringnorm.param2={
                        design <- model.matrix(~0+samples.IDs$tp53.status)
                        colnames(design) <- c("Mutated","WildType")
                        rcc.samples <- cbind(annots.df,rcc.samples)
                        rcc.samples$CodeClass[rcc.samples$CodeClass=="Housekeeping"] <- "Endogenous"
                        #rcc.samples$CodeClass[rcc.samples$Name%in%test.samp] <- "Housekeeping"
                        res.norm <- NanoStringNorm(rcc.samples,Background = 'mean.2sd',SampleContent = 'low.cv.geo.mean',CodeCount="sum",traits=design,round.values = TRUE,take.log = TRUE,return.matrix.of.endogenous.probes = TRUE)
                      },
                      desq2={
                        rcc.samples <- rcc.samples[!grepl("^Pos|^Neg",annots.df$CodeClass),]
                        samples.IDs <- samples.IDs[match(colnames(rcc.samples),samples.IDs$title),]
                        dds <- DESeqDataSetFromMatrix(countData = rcc.samples,
                                                      colData = samples.IDs,
                                                      design= ~ tp53.status)
                        dds <- DESeq(dds)
                        res.diff <- results(dds)
                        res.diff <- data.frame(deseq2=-log10(res.diff$padj),SYMBOL=row.names(res.diff))
                      },
                      nanostringDiff={
                        endogenous <- as.matrix(rcc.samples[grepl("Endog",annots.df$CodeClass),])
                        positive <- as.matrix(rcc.samples[grepl("Pos",annots.df$CodeClass),])
                        negative <- as.matrix(rcc.samples[grepl("Neg",annots.df$CodeClass),])
                        housekeeping <- as.matrix(rcc.samples[grepl("House",annots.df$CodeClass),])
                        NanoStringData <- createNanoStringSet(endogenous,positive,negative,housekeeping,samples.IDs)
                        pheno <- pData(NanoStringData)
                        group <- pheno$tp53.status
                        design.full=model.matrix(~0+group)
                        contrast <- c(1,-1)
                        NanoStringData <- estNormalizationFactors(NanoStringData)
                        result <- glm.LRT(NanoStringData,design.full,contrast=contrast)
                        res.diff <- data.frame(NSDiff=result$table,SYMBOL=row.names(result$table))
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
  if (!tool%in%c("desq2","nanostringDiff")){
    design <- model.matrix(~0+samples.IDs$tp53.status)
    colnames(design) <- c("Mutated","WildType")
    fit <- lmFit(res.norm,design)
    cm <- makeContrasts(diff=Mutated-WildType ,levels=design)
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    res.diff <- topTable(fit2, coef="diff",genelist=row.names(res.norm), number=Inf)
    res.diff <- data.frame(nanoR.total=-log10(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
    colnames(res.diff) <- c(tool,"SYMBOL")
  }
  return(res.diff)
}

make_designMatrix <- function(dataset,cond1 = "A", cond2 = "B",ncond1=(ncol(dataset)/2),ncond2=(ncol(dataset)/2)){
  status.control = rep(cond1,ncond1)
  status.test = rep(cond2,ncond2)
  status = c(status.control,status.test)
  design = model.matrix(~0+status)
  colnames(design) <- c(cond1,cond2)
  return(design)
}  

DEG_limma <- function(dataset,design){
  cm <- makeContrasts(diff = A-B, levels=design)
  fit <- lmFit(Micro,design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res.diff <- topTable(fit2, coef="diff",genelist=row.names(dataset), number=Inf)
  res.diff_limma <- data.frame(PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
  colnames(res.diff_limma) <- c("Limma","Gene.ID")
  return(res.diff_limma)
}

DEG_GEOlimma <- function(dataset,design){
  cont.matrix <- makeContrasts(constrast = A-B, levels=design)
  fit <- lmFit(Micro,design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  de <- topTable(fit22, number = nrow(Micro))
  res.diff_geolimma <- data.frame(PValue=(de$adj.P.Val),genes=row.names(de))
  colnames(res.diff_geolimma) <- c("GEOlimma","Gene.ID")
  return(res.diff_geolimma)
}



### Visualisation
PCA_tools <- function(data.to.comp){
  res.pca <- PCA(data.to.comp,graph=F)
  fviz_pca_ind (res.pca, col.ind = "cos2",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE)
}

UpsetPlot <- function(data.to.comp, threshold, log = FALSE){
  if(missing(threshold)){
    threshold = 0.05
  }
  Upset <- copy(data.to.comp)
  for (i in names(Upset)){
    for (u in 1:nrow(Upset)){
      if (log == F){
        if(data.to.comp[[i]][u] > threshold){
          Upset[[i]][u] = 0}
        else{
          Upset[[i]][u] = 1}
      }
      else{
        if(data.to.comp[[i]][u] < threshold){
          Upset[[i]][u] = 0}
        else{
          Upset[[i]][u] = 1}
      }
    }
  }

  Upset = as.data.frame(t(Upset))
  upset(Upset, sets = names(Upset), sets.bar.color = "#56B4E9",
          order.by = "freq", empty.intersections = "on" )

}


