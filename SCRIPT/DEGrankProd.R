#Votre choose directory
setwd("~/GIT/CPRD")

#Accès à la doc (introuvable depuis mon R)
browseVignettes("RankProd")

#Packages nécessaires à l'installation de RankProd
#Conseillé de faire avant DL : Dans le terminal Linux : 
#sudo apt-get update
#sudo apt-get install libgmp3-dev 
#sudo apt-get install libmpfr-dev

#Installation de RankProd

#Appel des packages pour le script
library("gmp")
library("Rmpfr")
library("RankProd")

#MICROARRAY-DNA---------------------------------------------------------------------------------------------------------
#source(file.path("./SCRIPT","Simulatemicroarraysseq.R"))

# paramétrage : le nom de la matrice, label contrôle/étudié (pour 2 classes), resultat reproductible,
# matrice log2 ou pas, devenir des NA, si nom de gènes dans la sortie, si rank product ou sum rank. Si sum rank,
# on pourra utiliser d'autres paramètres le réglage de la pvalue avec tail.time

#J'assigne les labels aux classes contrôle (g1) et étudié (g2) dans l'ordre de la structuration des colonnes
g1 = 6
g2 = 6
cl = rep(c(0,1),c(g1,g2))
cl # vérif'
# FALSE pour NA
outRPmic = RankProducts(Simulmicro, cl, rand = 123, logged = F , na.rm = T , calculateProduct = TRUE)
outRPmic2 = RankProducts(Simulmicro, cl, rand = 123, logged = T , na.rm = T , calculateProduct = TRUE)
outRPmic3 = RankProducts(Simulmicro, cl, rand = 123, logged = T , na.rm = F , calculateProduct = TRUE)
outRPmic4 = RankProducts(Simulmicro, cl, rand = 123, logged = F , na.rm = F , calculateProduct = TRUE)

valeur = outRPmic[["pfp"]]
value = outRPmic[["pval"]]
value

A = topGene(outRPmic,num.gene = 100)
head(A$Table1)
head(A$Table2)

B = topGene(outRPmic2,num.gene = 100)
head(B$Table1)
head(B$Table2)

C = topGene(outRPmic3,num.gene = 100)
head(C$Table1)
head(C$Table2)

D = topGene(outRPmic4,num.gene = 100)
head(D$Table1)
head(D$Table2)

#Tentative fonction (fonctionne a priori)------------------------------------------

ParaRP = function(data,g1,g2,ran,ProdRank,logged) {
  cl = rep(c(0,1),c(g1,g2))
  
  if (ProdRank == FALSE) {
      if (logged == TRUE) {
      outRPmic = RankProducts(data, cl, rand = ran, logged = TRUE , na.rm = TRUE , calculateProduct = FALSE)
    
      } else {
      outRPmic = RankProducts(data, cl, rand = ran, logged = FALSE , na.rm = TRUE , calculateProduct = FALSE) }

# 2nd fonction avec juste Rank Sum plutôt que ça avec ajout de paramètre libre en argument de fonction pour les 2
  } else {
      if (logged == TRUE) {
      outRPmic = RankProducts(data, cl, rand = ran, logged = TRUE , na.rm = TRUE) 
    
      } else {
      outRPmic = RankProducts(data, cl, rand = ran, logged = FALSE , na.rm = TRUE) }

  pvalue = outRPmic[["pval"]]
  return(pvalue)
  }
}

A = ParaRP(Simulmicro,6,6,123,FALSE,TRUE)

#MICROARRAY-AFFY.----------------------------------------------------------------------------------------------

#One origin dans notre cas 



#NANOSRING (NAPPA DEFAULT)-faisable ???????-----------------------------------------------------------------------------------------------------------------
library(readr)
NanoNormNappa <- read_csv("DATA/NANOSTRING/NanoNormNappa.csv")
View(NanoNormNappa)

