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

# Génération donnée microarray (afflymetrix non loggé ici)-----------------------------------------------------

source(file.path("./SCRIPT","Simulatemicroarraysseq.R"))

#MICROARRAY-DNA, pas pertinent ici car données affly------------------------------------------------------------------------------------------------

# paramétrage : le nom de la matrice, label contrôle/étudié (pour 2 classes), resultat reproductible,
# matrice log2 ou pas, devenir des NA, si nom de gènes dans la sortie, si rank product ou sum rank. Si sum rank,
# on pourra utiliser d'autres paramètres le réglage de la pvalue avec tail.time

#J'assigne les labels aux classes contrôle (g1) et étudié (g2) dans l'ordre de la structuration des colonnes
g1 = 6
g2 = 6
cl = rep(c(0,1),c(g1,g2))
cl
# execution du calcul
outRPmic = RankProducts(Simulmicro, cl, rand = 123, logged = FALSE , na.rm = FALSE , calculateProduct = TRUE)
#recup des p-value
value = outRPmic[["pval"]]

#MICROARRAY-AFFY.-----------------------------------------

#Une origin dans notre cas + 2 classes
# Cl :
g1 = 6
g2 = 6
cl = rep(c(0,1),c(g1,g2))
#Origine
origin = rep(1,g1+g2)
origin
#calcul DEG pour une origine, deux classes. Ces étapes sont inutiles si on a une seule origine, importante si plusieurs.
Simulmicro.sub = Simulmicro[,which(origin==1)]
Simulmicro.cl.sub = cl[which(origin==1)]
Simulmicro.origin.sub = origin[which(origin==1)]
#RP
RP.out = RankProducts(Simulmicro.sub,Simulmicro.cl.sub, logged=FALSE, na.rm=FALSE,plot=FALSE,  rand=123)
PvalueRPaffy = RP.out[["pval"]]
colnames(PvalueRPaffy) = c("RP p-val surexpression (c1 < c2)","RP p-val sous-expression (c1 > c2)")
#RS
RS.out = RankProducts(Simulmicro.sub,Simulmicro.cl.sub, logged=FALSE, na.rm=FALSE,plot=FALSE,  rand=123, calculateProduct = FALSE)
PvalueRSaffy = RS.out[["pval"]]
colnames(PvalueRSaffy) = c("RS p-val surexpression (c1 < c2)","RS p-val sous-expression (c1 > c2)")

#comparaison
PValueCompare = cbind(PvalueRPaffy,PvalueRSaffy)
PValueCompare = PValueCompare[,c(1,3,2,4)]

#Résultat variable, ce n'est pas exactement la même chose.

#NANOSRING (NAPPA DEFAULT)-------------------------------------------------------------------------------------

library(readr)
NanoNormNappa <- read_csv("DATA/NANOSTRING/NanoNormNappa.csv")
View(NanoNormNappa)

g1 = length(NanoNorNappabis)/2 #attention, pas toujours la même taille d'échantillonnage de part et d'autres !
g2 = length(NanoNorNappabis)/2
cl = rep(c(0,1),c(g1,g2))

names.gene = NanoNormNappa[,1]
NanoNorNappabis = NanoNormNappa[,-1]

#RP
RP.outNano = RankProducts(as.matrix(NanoNorNappabis), cl, logged=TRUE, na.rm = FALSE ,plot=FALSE, gene.names = names.gene, rand=123)
#RS
RS.outNano = RankProducts(as.matrix(NanoNorNappabis), cl, logged=TRUE, na.rm = FALSE ,plot=FALSE, gene.names = names.gene, rand=123)

# summary(NanoNorNappabis) (V19 et V24 contiennent des valeurs négatives d'après les valeurs min)

# write.csv(NanoNorNappabis,"~/Bureau/M1/ProjetRD/NanoNappaBis.csv", row.names = TRUE)

# Ne marche pas ( pb données + matrice) 
