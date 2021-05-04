# Nom               : Nanostring_functions.R
# Type              : Programme
# Objet             : Normalisation + analyse DEG en 2 fonctions séparées
# Input             : Dataset Nanostring
# Output            : data.to.comp
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 18.04.2021
#______________________________________________________________________________

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
  
  if (!tool%in%c("desq2")){
    group = table(samples.IDs$tp53.status)
    n1 = as.integer(group[1])
    n2 = as.integer(group[2])
  }
  if (tool%in%c("RankProduct", "RankSum")){
    design = rep(c(0,1),c(n1,n2))
  }


  tools.fnc <- switch(tool,
                      desq2={
                        rcc.samples <- rcc.samples[!grepl("^Pos|^Neg",annots.df$CodeClass),]
                        samples.IDs <- samples.IDs[match(colnames(rcc.samples),samples.IDs$title),]
                        dds <- DESeqDataSetFromMatrix(countData = rcc.samples,
                                                      colData = samples.IDs,
                                                      design= ~ tp53.status)
                        
                        dds <- DESeq(dds)

                        resG <- results(dds, altHypothesis="greater", format = "DataFrame")
                        resL <- results(dds, altHypothesis="less", format = "DataFrame")
                        head(resG)
                        res.diff = data.frame("DESeq2_Up"=(resG$padj),"DESeq2_Down"=(resL$padj), "SYMBOL" = row.names(resG))
                        
                        #res.diff <- data.frame(deseq2=-log10(res.diff$padj),SYMBOL=row.names(res.diff))
                      },

                      limma = {

                        design <- model.matrix(~0+samples.IDs$tp53.status)
                        colnames(design) <- c("A","B")
                        
                        res.diff = DEG_limma(data,design)
                        res.diff = DEG_alternative(res.diff)
                        
                        method_Up = paste0(tool_norm,"_",tool,"_Up")
                        method_Down = paste0(tool_norm,"_",tool,"_Down")

                        #fit <- lmFit(data,design)
                        #fit2 <- contrasts.fit(fit, cm)
                        #fit2 <- eBayes(fit2)
                        #res.diff <- topTable(fit2, coef="diff",genelist=row.names(data), number=Inf)
                        #res.diff <- data.frame(nanoR.total=-log10(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
                        res.diff <- data.frame(A=(res.diff$PValue_Up), B=(res.diff$Pvalue_Down) ,SYMBOL=res.diff$SYMBOL)
                        colnames(res.diff) = c(method_Up,method_Down, "SYMBOL")

                      },
                      
                      Wilcox = {
                        res.diff = wilcoxDEG(data, n1, n2)
                        res.diff = res.diff[(-1)]
                        
                        res.diff$SYMBOL = row.names(res.diff)
                        method_less = paste0(tool_norm,"_",tool,"_less")
                        method_greater = paste0(tool_norm,"_",tool,"_greater")
                        colnames(res.diff) = c(method_less, method_greater, "SYMBOL")
                      },
                      
                      RankProduct = {
                        res.diff = RankProducts(data, design, rand = 123, logged = TRUE , na.rm = TRUE , calculateProduct = TRUE)
                        res.diff = res.diff[["pval"]]
                        
                      },
                      
                      RankSum = {
                        res.diff = RankProducts(data, design, rand = 123,logged = TRUE ,na.rm = TRUE ,calculateProduct = FALSE)
                        res.diff = res.diff[["pval"]]
                      },
                      
                      stop("Enter something that switches me!") 
          
  )
  if (tool%in%c("RankProduct","RankSum")){
    res.diff = as.data.frame(res.diff)
    res.diff$SYMBOL = row.names(data)
    method_Up = paste0(method,"_Up")
    method_Down = paste0(method,"_Down")
    colnames(res.diff) = c(method_Up,method_Down,"SYMBOL")
  }
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
#tool_norm = "nappa.NS" 
#tool = "RankSum"
#data = tools_norm.inspect(raw.data, "nappa.NS")

tools_norm <- c("nappa.NS","nappa.param1", "nappa.param2","nappa.param3","nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2","nanoR.top100","nanoR.total","nanostringR")
tools_norm <- c("nappa.NS","nanostringnorm.default","nanoR.top100","nanoR.total","nanostringR")

tools_diff = c("limma","Wilcox","RankProduct","RankSum" )
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
Upset = UpsetPlot(data.to.comp,threshold = 0.05)
methods = colnames(Upset)

Upreg = methods[grepl("Up|less",methods)]
Downreg = methods[grepl("Down|greater",methods)]

dev.new()
upset(Upset, sets = Upreg, sets.bar.color = "#56B4E9",
             order.by = "freq", 
      empty.intersections = NULL )

upset(Upset, sets = Downreg, sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = NULL )
