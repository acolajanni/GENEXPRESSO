# Nom               : geolimma.R
# Type              : Programme
# Objet             : test l'analyse DEG de geolimma
# Input             : Dataset normalisé de Nanostring + dataset microarray
# Output            : analyse DEG 
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 14.04.2021
#______________________________________________________________________________

source(file.path("./SCRIPT","functions.R"))


# Problème 1 : GEOlimma ??? Gene Expression Omnibus - Articles récents mais pas de package dispo
Micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)

Micro = read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)

#________________________________________________ 
#Limma : 
# Nano : Déja implémenté

#________________________________________________ 
#Limma : 
# Microarray : (extrait de functions.R)

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

#Récupération des pvalues
res.diff_limma = DEG_limma(Micro,design)
head(res.diff_limma)
tail(res.diff_limma)
#________________________________________________
# GEOlimma sur Micro array

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

# récupération des pvalues
res.diff_GEOlimma= DEG_GEOlimma(Micro,design)
head(res.diff_GEOlimma)
tail(res.diff_GEOlimma)




# Test nanostring - limma vs geolimma
  Micro = read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)
  
  design <- model.matrix(~0+samples.IDs$tp53.status)
  colnames(design) <- c("Mutated","WildType")
  fit <- lmFit(Micro,design)
  cont.matrix <- makeContrasts(diff=Mutated-WildType ,levels=design)

  cont.matrix <- makeContrasts(constrast = A-B, levels=design)
  fit <- lmFit(Micro,design)
  fit2  <- contrasts.fit(fit, cont.matrix)
  load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  de <- topTable(fit22, number = nrow(Micro))
  res.diff_geolimma <- data.frame(PValue=(de$adj.P.Val),genes=row.names(de))
  colnames(res.diff_geolimma) <- c("GEOlimma","Gene.ID")
############# geolimma
  design <- model.matrix(~0+samples.IDs$tp53.status)
  colnames(design) <- c("A","B")
  fit <- lmFit(Micro,design)
  cont.matrix <- makeContrasts(diff=A-B,levels=design)
  res.diff_limma = DEG_limma(Micro,design)
  colnames(res.diff_limma) <- c("imma","Gene.ID")
  data_to_comp = merge(res.diff_geolimma,res.diff_limma,by = "Gene.ID",all=T)
  
# Conclusion : résultats exactement identiques

