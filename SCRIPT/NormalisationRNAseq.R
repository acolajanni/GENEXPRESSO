setwd("~/GIT/CPRD")



library("DESeq2")
library("Biobase")
library("cqn")

conditions <-factor(c("condition1","condition1","condition1","condition1","condition1","condition1","condition2", "condition2","condition2","condition2","condition2","condition2"))
SimulRNASEQ = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)

#il faut écrire single ou paired pour les reads, ça dépend du type de méthodologie utilisée, 12 car 12 échantillons dans ce jeu
type = factor(rep("single-read",12))

colData = data.frame(conditions,type,row.names=colnames(SimulRNASEQ))

#Pour utiliser les fonctions de normaliser, de DE avec DESeq, il faut créer un objet où on doit mettre notre dataframe dedans
dds<-DESeqDataSetFromMatrix(SimulRNASEQ,colData,design=~conditions)
class(dds)
colData(dds)
design(dds)

#petite visu sur les données si c'est toujours ok
dim(counts(dds))
head(counts(dds))
summary(counts(dds))
#-----------------------------------------------------------------------------
#Normalisation (median of ratio) + DEG
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)

DEG_Wald <- DESeq(dds, test = "Wald")
results(DEG_Wald)

DEG_LRT <- DESeq(dds, test = "LRT",reduced = ~1)
results(DEG_LRT)
# Likelihood ratio test 

# Lignes de codes issues de functions.R
res.diff <- results(DEG_Wald)
res.diff <- data.frame(deseq2=-log10(res.diff$padj),SYMBOL=row.names(res.diff))
res.diff

res.diff2 <- results(DEG_LRT)
res.diff2 <- data.frame(deseq2=-log10(res.diff2$padj),SYMBOL=row.names(res.diff2))
summary(res.diff2)

colData(dds)
#
sizeF=sizeFactors(dds)

#J'ai chois de stocker les données dans une dataframe, car avec 100 gènes on dépasse la capacité d'affichage de la console de Rstudio
#et c'est beaucoup plus clair
StockVisuNorm = data.frame(counts(dds,normalized = TRUE))

#simple stat descriptive
summary( counts ( dds ,  normalized = TRUE))
counts ( dds ,  normalized = TRUE)


#comptage moyen gène et son log2 ratio
par(mfrow=c(1 ,1))
DESeq2::plotMA(dds)

#------------------------------------------------------------------------------
# Méthode CQN

GC= round(runif(100, min=0.01, max=0.99), digits=4)
Lenght = round(runif(100, min=1000, max=9000))
help(cqn)
cqn = cqn(SimulRNASEQ, x = GC,  lengths = Lenght, sizeFactors = sizeF, verbose = TRUE)

cqnOffset <- cqn$glm.offset
cqnNormFactors <- exp(cqnOffset)
cqnNormFactors

#normFactors_sameScale <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))

#fusion RPKM + CQN
RPKM.cqn = cqn$offset
RPKM.cqn
RPKM.cqn = cqn$y + cqn$offset
RPKM.cqn = as.data.frame(RPKM.cqn) # = jeu de données normalisé ???

########################" Essai sans réussite

conditions <-factor(c("condition1","condition1","condition1","condition1","condition1","condition1","condition2", "condition2","condition2","condition2","condition2","condition2"))
colData = data.frame(conditions,type,row.names=colnames(RPKM.cqn))
print(class(RPKM.cqn))

RPKM.cqn<-DESeqDataSetFromMatrix(RPKM.cqn,colData, ~conditions)
RPKM.cqn <- DESeqDataSetFromMatrix(round(RPKM.cqn), 
                                  colData=colData, 
                                  design=~conditions)

# DESeq2 a besoin de renormaliser derrière donc inutile
# essayons avec EdgeR
# MAis 
d.mont <- DGEList(counts = montgomery.subset, lib.size = sizeFactors.subset,                   
                    group = rep(1:2,each=ncol(data)/2))

#---------------------------------------------------------------------------------
# RPKM ? génération nbR
# reads per kilobases per millions
