# Nom               : testPaul.R
# Type              : Programme Test, Brouillon
# Objet             : Analyse coexpression Kendall + Spearman
# Input             : Dataset (microarray, RNA-seq, Nanostring)
# Output            : iplot
# Auteur            : Paul SOURBE
# R version         : 3.6
# Date de creation  : 20.04.2020
#______________________________________________________________________________
library("reshape2")
#install.packages('igraph')
library('igraph')
library(data.table)
source(file.path("./SCRIPT","functions.R"))

data = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)

#Transposition de la matrice car cor() agit sur les colonnes et non les lignes
data = as.data.frame(t(data))

################# KENDALL ######################
testKendall <- cor(data, method = "kendall")

################# SPEARMAN #####################

testSpearman <- cor(data, method = "spearman")

################################################################################

################# WGCNA ####################
##### corAndPvalue #####
library("WGCNA")
testWGCNA <- corAndPvalue(data)

##### TOM #####
# Mesure de chevauchement topologique (TOM) :
# est une mesure de similarité par paire entre les nœuds du réseau (gènes).
# TOM(i,j ) est élevé si les gènes i,j ont de nombreux voisins communs (le chevauchement de leurs voisins de réseau est important).
# TOM(i,j ) élevé implique que les gènes ont des profils d'expression similaires.
# TOM, en tant que mesure de similarité, peut être transformé en une mesure de dissimilarité distTOM = 1-TOM.
TOMsimilarityFromExpr()
TOMsimilarity()
hierarchicalConsensusTOM()
TOMplot()





