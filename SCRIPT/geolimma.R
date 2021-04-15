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
source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))

# Problème 1 : GEOlimma ??? Gene Expression Omnibus - Articles récents mais pas de package dispo
Micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv", header = TRUE,row.names = 1)

Nano = read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)

make_designMatrix <- function(dataset,cond1 = "Control", cond2 = "Test"){
  status.control = rep(cond1,ncol(dataset)/2)
  status.test = rep(cond2,ncol(dataset)/2)
  status = c(status.control,status.test)
  design = model.matrix(~0+status)
  colnames(design) <- c(cond1,cond2)
  return(design)
}

DEG_limma <- function(dataset,design){
  fit <- lmFit(Micro,design)
  cm <- makeContrasts(diff = Control-Test,levels=design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res.diff <- topTable(fit2, coef="diff",genelist=row.names(dataset), number=Inf)
  res.diff_limma <- data.frame(PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
  colnames(res.diff_limma) <- c("Limma","Gene.ID")
  return(res.diff_limma)
}

#________________________________________________ 
#Limma : 
# Nano : Déja implémenté

#________________________________________________ 
#Limma : 
# Microarray : (extrait de functions.R)

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

#Application du modèle linéaire:
fit <- lmFit(Micro,design)
cm <- makeContrasts(diff=Control-Test ,levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
res.diff <- topTable(fit2, coef="diff",genelist=row.names(Micro), number=Inf)
res.diff_limma <- data.frame(PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
colnames(res.diff_limma) <- c("Limma","Gene.ID")
head(res.diff_limma)
tail(res.diff_limma)

#________________________________________________
# GEOlimma sur Micro array

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

# Lignes de codes issues de DE_source.R : Application du modèle linéaire
cont.matrix <- makeContrasts(compare="Control-Test",levels=design)
fit2  <- contrasts.fit(fit, cont.matrix)
load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
de <- topTable(fit22, number = nrow(Micro))
res.diff_geolimma <- data.frame(PValue=(de$adj.P.Val),genes=row.names(de))
colnames(res.diff_geolimma) <- c("GEOlimma","Gene.ID")
head(res.diff_geolimma)
tail(res.diff_geolimma)



  
