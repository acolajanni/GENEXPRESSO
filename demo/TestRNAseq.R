library(edgeR)
library(DESeq)
library(DESeq2)
library(limma)
library(voom)

#########################
#########################
#########################
#########################
# EdgeR avec la normalisation de DESeq2 : 
data = Simul.data("RNAseq",n.cond1 = 15, n.cond2 = 15,nb.genes = 1000)

condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,type,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

# Calcul des size factors de DESeq2
SizeF = estimateSizeFactors(dds)
SizeF = sizeFactors(SizeF)
names(SizeF) = NULL

# Analyse de EdgeR
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2),
              norm.factors = SizeF) # functions for model with n1 = n2 


Disp = estimateCommonDisp(DGE)
Disp = estimateTagwiseDisp(Disp)

# Application de size factors de DESeq2
Disp$samples$norm.factors = SizeF
# Vérification
Disp$samples #ça fonctionne

pvalue = exactTest(Disp)
Norm.DESeq.ExactTest = topTags(pvalue)
Norm.DESeq.ExactTest = Norm.DESeq.ExactTest$table

#########################
#########################
#########################
#########################
# DESeq2 avec la normalisation de TMM de edgeR 

data = Simul.data("RNAseq",n.cond1 = 15, n.cond2 = 15,nb.genes = 1000)

# Dataset into DEGlist (edgeR class)
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 
# Normalization 
res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")
NormF = res.norm$samples$norm.factors

condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,type,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

dds = estimateSizeFactors(dds)
##
avant = sizeFactors(dds)
sizeFactors(dds) = NormF
##
apres = sizeFactors(dds)
# On vérifie les changements : 
avant
apres # ça fonctionne


DEG <- DESeq(dds, test = "Wald")
DEG <- results(DEG)

DEG= data.frame(DEG$log2FoldChange, DEG$padj, row.names = row.names(DEG))
Norm.TMM.Wald = DEG[order(DEG$DEG.padj),]
Norm.TMM.Wald = Norm.TMM.Wald[1:10,]



######## Méthodes standard : 
#### EdgeR TMM

# Dataset into DEGlist (edgeR class)
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 
# Normalization 
res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")

Disp = estimateCommonDisp(res.norm)
Disp = estimateTagwiseDisp(Disp)
# Retrieving pvalues
pvalue = exactTest(Disp)
TMM = topTags(pvalue)$table



#### DESeq2 Wald

#data into DESeqDataSet (DESeq2 class) :
condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )

type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

# Retrieve DEG analysis (test Wald)
DEG <- DESeq(dds, test = "Wald")
DEG <- results(DEG)

DEG= data.frame(DEG$log2FoldChange, DEG$padj, row.names = row.names(DEG))
Wald = DEG[order(DEG$DEG.padj),]
Wald = Wald[1:10,]



#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 

count.matrix = table.count.subset

design = factor(tab.subset$`Sample Type`)
class = ifelse(design == "Primary Tumor", yes = 1, no = 2)
class = design

A = tools.norm.RNAseq(table.count.subset, tool = "TMM", design = class)
B = tools.DEG.RNAseq(table.count.subset, nf = A,
                     tool.norm = "TMM", 
                     tool.DEG = "GLM", 
                     design = class)


tools.norm.RNAseq <- function(data, tool, design){
  
  
  if (tool%in%c("TMM", "TMMwsp", "RLE", "Upperquartile", "voom")){
    edgeR.dgelist = DGEList(counts = count.matrix, group = factor(class))
  }
  
  
  tools_norm_RNAseq.fnc <- switch(tool,
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
                                  
                                  vst = {
                                    DESeq.cds = newCountDataSet(countData = count.matrix,
                                                                conditions = factor(class))
                                    DESeq.cds = estimateSizeFactors(DESeq.cds)
                                    nf = sizeFactors(DESeq.cds)
                                    
                                  },
                                  
                                  voom = {
                                    nf = calcNormFactors(count.matrix, method = "TMM")
                                    voom.data = voom(count.matrix, 
                                                     design = model.matrix(~factor(class)),
                                                     lib.size = colSums(count.matrix) * nf)
                                    return(voom.data)
                                  },
                                  stop("Enter a normalization method that switches me !")
              
                                  
  )
  if (tool%in%c("TMM", "TMMwsp", "RLE", "Upperquartile")){
    nf = nf[["samples"]][["norm.factors"]]
    names(nf) = colnames(count.matrix)
  }
  return(nf)
}

tools.DEG.RNAseq <- function(count.matrix.raw, nf, tool.norm, tool.DEG, design){
  # limma only apply for voom normalization
  if (class(nf) == "EList"){
    tool.norm = "voom"
    tool.DEG = "limma"
    voom.data = nf 
    
    method = "limma.voom"
  }
  else{
    method = paste0(tool.nom,"+",tool.DEG)
  }
  
  # Specific class is needed for edgeR analysis
  if (tool.DEG %in% c("ExactTest","GLM")){
    edgeR.dgelist = DGEList(counts = count.matrix.raw, group = factor(design))
    
    # Estimating dispersion, retrieve and apply normalization factors
    edgeR.dgelist[["samples"]][["norm.factors"]] = nf
    edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
    edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist,
                                        trend = "movingave")
    
  }
  # Same for DESeq2
  else if (tool.DEG%in%c("nbinom.Wald","nbinom.LRT")){
    design = data.frame(design,row.names=colnames(count.matrix.raw))
    design$design = as.factor(design$design)
    dds<-DESeqDataSetFromMatrix(count.matrix.raw,
                                colData = design,
                                design= ~design)
    sizeFactors(dds) = nf
    
    
    dds = estimateDispersions(dds, sharingMode = "maximum", 
                              method = "pooled",
                              fitType = "local")
    
  }
  # Same again for DESeq
  else if (tool.DEG == "nbinom"){
    DESeq.cds = newCountDataSet(countData = count.matrix.raw,
                                conditions = factor(design))
    # Initializing size factors
    sizeFactors(DESeq.cds) = nf
  }
  
  tools_norm_RNAseq.fnc <- switch(tool.DEG,
                                  # EdgeR exact test
                                  ExactTest = {
                                    
                                    DEG = exactTest(edgeR.dgelist)
                                    DEG.pval = DEG$table$PValue
                                    
                                    
                                  },
                                  # EdgeR General linear model
                                  GLM = {
                                    design = model.matrix(~0+group, data = edgeR.dgelist$samples)
                                    colnames(design) <- levels(edgeR.dgelist$samples$group)
                                    
                                    # Fitting the GLM model
                                    fit <- glmQLFit(edgeR.dgelist, design)
                                    # Testing outliers
                                    qlf <- glmQLFTest(fit, contrast=c(-1,1))
                                    # Retrieving pvalues
                                    DEG.pval = qlf[["table"]]$PValue

                                  },
                                  # DESeq analysis
                                  nbinom = {
                                    # Dispersions
                                    DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum",     
                                                                    method = "pooled", 
                                                                    fitType = "local")
                                    
                                    # Searching for DE genes
                                    DESeq.test = nbinomTest(DESeq.cds, "1", "2")
                                    DEG.pval = DESeq.test$pval
                                    
                                  },
                                  # DESeq2 Wald Test
                                  nbinom.Wald = {
                                    DEG <- DESeq(dds, test = "Wald")
                                    
                                  },
                                  # DESeq2 Likelihood ratio test for GLMs
                                  nbinom.LRT = {
                                    DEG <- DESeq(dds, test = "LRT")
                                    
                                  },
                                  limma = {
                                    voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(design)))
                                    voom.fitbayes = eBayes(voom.fitlimma)
                                    DEG.pval = voom.fitbayes$p.value[, 2]
                                    #voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")
                                    #return(voom.fitbayes)
                                  },
                                  stop("Enter something that switches me !")
  )
  
  if (tool.DEG%in% c("nbinom.Wald","nbinom.LRT")){
    DEG.pval = results(DEG)$pvalue
  }
  
  DEG.padj = p.adjust(DEG.pval, method = "BH")
  
  DEG = data.frame(SYMBOL = row.names(count.matrix.raw), DEG.padj)
  colnames(DEG) = c("SYMBOL", method)
  
  return(DEG)
}








