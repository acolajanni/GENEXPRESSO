library(edgeR)
library(DESeq)
library(DESeq2)

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

# Analyse de EdgeR
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 


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
colData = data.frame(condition,type,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

# Retrieve DEG analysis (test Wald)
DEG <- DESeq(dds, test = "Wald")
DEG <- results(DEG)

DEG= data.frame(DEG$log2FoldChange, DEG$padj, row.names = row.names(DEG))
Wald = DEG[order(DEG$DEG.padj),]
Wald = Wald[1:10,]
