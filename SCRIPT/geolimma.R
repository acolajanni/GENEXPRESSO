# Nom               : geolimma.R
# Type              : Programme
# Objet             : test l'analyse DEG de geolimma
# Input             : Dataset normalisé de Nanostring + dataset microarray
# Output            : analyse DEG 
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 14.04.2021
#______________________________________________________________________________
library(limma)

# Problème 1 : GEOlimma ??? Gene Expression Omnibus - Articles récents mais pas de package dispo
Micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv", header = TRUE,row.names = 1)

Nano = read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)

#________________________________________________ 
#Limma : 
# Nano : (extrait de functions.R)
design <- model.matrix(~0+samples.IDs$tp53.status)
colnames(design) <- c("Mutated","WildType")
fit <- lmFit(res.norm,design)
cm <- makeContrasts(diff=Mutated-WildType ,levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
res.diff <- topTable(fit2, coef="diff",genelist=row.names(res.norm), number=Inf)
res.diff <- data.frame(nanoR.total=-log10(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
colnames(res.diff) <- c(tool,"SYMBOL")
#________________________________________________ 
#Limma : 
# Microarray : (extrait de functions.R)

# construction de la matrice de design
status.control = rep("Control",ncol(Micro)/2)
status.test = rep("Test",ncol(Micro)/2)
status = c(status.control,status.test)
design = model.matrix(~0+status)
colnames(design) <- c("Control","Test")

#Application du modèle linéaire:
fit <- lmFit(Micro,design)
cm <- makeContrasts(diff=Control-Test ,levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
res.diff <- topTable(fit2, coef="diff",genelist=row.names(Micro), number=Inf)
res.diff <- data.frame(PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
colnames(res.diff) <- c("Limma","Gene.ID")
res.diff
