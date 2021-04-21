# Nom               : Nanostring_functions.R
# Type              : Programme
# Objet             : Normalisation + analyse DEG en 2 fonctions séparées
# Input             : Dataset Nanostring
# Output            : data.to.comp
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 18.04.2021
#______________________________________________________________________________

source(file.path("./SCRIPT","functionwilcoxtest.R"))
source(file.path("./SCRIPT","functions.R"))

tools_norm.inspect <- function(raw.data,tool,nanoR=F){
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

tools_DEG.inspect <- function(raw.data, tool, data, tool_norm) {

  rcc.samples <- raw.data$rcc.df
  annots.df <- raw.data$annots.df
  samples.IDs <- raw.data$samples.IDs
  FOV <- raw.data$FOV

  method = paste0(tool_norm,"_",tool)


  tools.fnc <- switch(tool,
                      desq2={
                        rcc.samples <- rcc.samples[!grepl("^Pos|^Neg",annots.df$CodeClass),]
                        samples.IDs <- samples.IDs[match(colnames(rcc.samples),samples.IDs$title),]
                        dds <- DESeqDataSetFromMatrix(countData = rcc.samples,
                                                      colData = samples.IDs,
                                                      design= ~ tp53.status)
                        dds <- DESeq(dds)
                        res.diff <- results(dds)
                        #res.diff <- data.frame(deseq2=-log10(res.diff$padj),SYMBOL=row.names(res.diff))
                        res.diff <- data.frame(deseq2=(res.diff$padj),SYMBOL=row.names(res.diff))
                        colnames(res.diff) = c("DESeq2", "SYMBOL")
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
                      
                      limma = {
                        group = table(samples.IDs$tp53.status)
                        n1 = as.integer(group[1])
                        n2 = as.integer(group[2])
                        
                        design <- model.matrix(~0+samples.IDs$tp53.status)
                        colnames(design) <- c("A","B")
                        
                        res.diff = DEG_limma(data,design)

                        #fit <- lmFit(data,design)
                        #fit2 <- contrasts.fit(fit, cm)
                        #fit2 <- eBayes(fit2)
                        #res.diff <- topTable(fit2, coef="diff",genelist=row.names(data), number=Inf)
                        #res.diff <- data.frame(nanoR.total=-log10(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
                        res.diff <- data.frame(limma=(res.diff$PValue),SYMBOL=res.diff$SYMBOL)
                        colnames(res.diff) = c(method, "SYMBOL")

                      },
                      
                      Wilcox = {
                        group = table(samples.IDs$tp53.status)
                        n1 = as.integer(group[1])
                        n2 = as.integer(group[2])
                        res.diff = wilcoxDEG(data, n1, n2)
                        res.diff = res.diff[-(-1)]
                        res.diff$SYMBOL = row.names(res.diff)
                        colnames(res.diff) = c(method, "SYMBOL")
                      },
                      
                      stop("Enter something that switches me!") 
          
  )
  
  return(res.diff)
}


#design <- model.matrix(~0+samples.IDs$tp53.status)
#colnames(design) <- c("Mutated","WildType")
##



# Appel de fonction
# Load data
raw.data = readRDS(file = "./DATA/NANOSTRING/Nanostring_Data.rds" )
##
data.dir <- "./DATA/NANOSTRING"
RCC.dir <- file.path(data.dir,"GSE146204_RAW")
raw <- RCC.dir
samples.IDs <- raw.data$samples.IDs
##

tools_norm <- c("nappa.NS","nappa.param1", "nappa.param2","nappa.param3","nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2","nanoR.top100","nanoR.total","nanostringR")
tools_diff = c("limma", "Wilcox" )
data.to.comp <- tools_DEG.inspect(raw.data = raw.data,data =  raw.data, tool = "desq2", tool_norm = NULL)

for (tool_norm in tools_norm){
  nanoR=F
  raw <- raw.data
  if (tool_norm%in%c("nanoR.top100","nanoR.total")){
    nanoR <- T
    RCC.dir <- file.path(data.dir,"GSE146204_RAW")
    raw <- RCC.dir
  }
  tmp <- tools_norm.inspect(raw,tool_norm,nanoR)

  for (tool_diff in tools_diff){
    print(paste(c(tool_norm, tool_diff)))
    tmp2 = tools_DEG.inspect(raw.data = raw.data, data =  tmp, tool = tool_diff, tool_norm = tool_norm)
    data.to.comp <- merge(data.to.comp,tmp2,by="SYMBOL",all=T) 
  }
}

row.names(data.to.comp) <- data.to.comp$SYMBOL
data.to.comp <- data.to.comp[,-1]

# if you need to remove genes not present in all analyses: NA (missing) data
data.to.comp <- na.omit(data.to.comp)
data.to.comp <- as.data.frame(t(data.to.comp))

################################################################################
################################################################################
# compute PCA
PCA_tools(data.to.comp)

################################################################################
################################################################################
# compute intersection with upsetR
#UpsetPlot(data.to.comp, threshold = 0.05, log = F)

### ATTENTION, si vous lancez comme ça, vous risquez de planter
A = UpsetPlot(data.to.comp,threshold = 0.05)
names = colnames(A)

upset(A, sets = names[1:11], sets.bar.color = "#56B4E9",
             order.by = "freq", 
      empty.intersections = "on" )

upset(A, sets = names[12:21], sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on" )

