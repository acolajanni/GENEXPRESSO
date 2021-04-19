setwd("~/GIT/CPRD")

BiocManager::install(c("Biobase","DESeq2", "cqn"))

source(file.path("./SCRIPT","SimulateRNAseq.R"))

library("DESeq2")
library("Biobase")
library("cqn")

SimulRNASEQ = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
conditions = factor(c("condition 1","condition 1","condition 1","condition 1","condition 1","condition 1","condition 2", "condition 2","condition 2","condition 2","condition 2","condition 2"))

#il faut écrire single ou paired pour les reads, ça dépend du type de méthodologie utilisée, 12 car 12 échantillons dans ce jeu
type = factor(rep("single-read",12))

colData = data.frame(conditions,type,row.names=colnames(SimulRNASEQ))

#Pour utiliser les fonctions de normaliser, de DE avec DESeq, il faut créer un objet où on doit mettre notre dataframe dedans
dds=DESeqDataSetFromMatrix(SimulRNASEQ,colData,design=~conditions)
class(dds)
colData(dds)
design(dds)

#petite visu sur les données si c'est toujours ok
dim(counts(dds))
head(counts(dds))
summary(counts(dds))

#Normalisation (median of ratio?)
dds <- DESeq(dds)
colData(dds)
# "effective library size" (sans connaitre la profondeur de base ?)
sizeF=sizeFactors(dds)

#J'ai chois de stocker les données dans une dataframe, car avec 100 gènes on dépasse la capacité d'affichage de la console de Rstudio
#et c'est beaucoup plus clair
StockVisuNorm = data.frame(counts(dds,normalized = TRUE))

#simple stat descriptive
summary( counts ( dds ,  normalized = TRUE))

#comptage moyen gène et son log2 ratio ?
par(mfrow=c(1 ,1))
DESeq2::plotMA(dds)

#faire un boxplot: extraire donnée avec "rld"
#------------------------------------------------------------------------------
# Méthode CQN  (inutile) même échantillon ; entre gène

#on génère du Chargaff et des longueurs aléatoires de gènes pour la normalisation
GC= round(runif(100, min=0.01, max=0.99), digits=4)
Lenght = round(runif(100, min=1000, max=9000))

#object avec les facteurs, sizeF réutilisé 
cqn = cqn(SimulRNASEQ, x = GC,  lengths = Lenght, sizeFactors = sizeF, verbose = TRUE)
#normalisation
cqnOffset <- cqn$glm.offset
cqnNormFactors <- exp(cqnOffset)8

#normGeom <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors))) moy géométrique

#--------------------------------------------------------------------------------

#fusion RPKM + CQN    #idem

#RPKM.cqn = cqn$y + cqn$offset (marche mais pas cohérent, besoin nb read normlt)

# RPKM ? génération nbR

