# Nom               : Test_norm_RNAseq
# Type              : Programme
# Objet             : Tester la normalisation (les fonctions iront dans functions.R)
# Input             : dataset RNA-seq
# Output            : ?
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 15.03.2021
#______________________________________________________________________________

library(edgeR)
library(DEFormats)

#help("read.csv")

# importation de notre jeu de données
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
#Besoin de stocker le jeu de données sous forme de matrice
data.matrix(data)

#help("readDGE")
#help('DGEList')
#help('simulateRnaSeqData')

# On créer une DGElist pour pouvoir faire d'autres analyses
data = DGEList(count = data, 
        group = rep(1:2,each=ncol(data)/2)) #2 groupes (control / test) ==> 2 rep(1:2) / chaque grp a autant d'indiv chacun (ncol(data)/2)

#help("calcNormFactors.DGEList")

# On peut maintenant calculer des facteurs de normalisation : (on les calcule juste, on les applique pas donc voir comment faire ?)

# Méthode par défaut : TMM
TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")
TMM_norm

# Methode utiles si le jeu de données est rempli de 0
TMMwsp_norm <- calcNormFactors.DGEList(data, method = "TMMwsp") 
TMMwsp_norm

# Relative log expression ==> moyenne geometric
RLE_norm <- calcNormFactors.DGEList(data, method = "RLE")
RLE_norm

# méthode normalisation par rapport au 75% quantile
upperquartile_norm <- calcNormFactors.DGEList(data, method = "upperquartile")
upperquartile_norm

#Maintenant, besoin d'estimer la dispersion puis tagwise dispersion
# Dispersion : permet d'établir un critère de dispersion des données de reads mappé sur un gène 
# Sans ce paramètre on peut pas établir si la différence de read mappé est due à une différence de longueur de gène
# Ensuite : Tagwise : établit la dispersion des données pour chacune des valeurs du dataset
Disp = estimateCommonDisp(TMM_norm)
Disp = estimateTagwiseDisp(Disp)


help("estimateCommonDisp")
help('estimateTagwiseDisp')

help(exactTest)
help(topTags)

# Test des DEG avec une méthode proche du Test exact de fisher
DEG = exactTest(Disp)
DEG
# on affiche par classement, les gènes les plus différentiellement exprimés
topTags(DEG, sort.by = 'PValue')

# page 22/122 sur la doc de EdgeR

# Problème : Exporter la normalisation :
# Sur edgeR, il calcule un "facteur de normalisation" qu'il applique dans les calculs (Dispersion, etc.)
# Mais ce facteur ne transforme pas le jeu de données. Donc comment faire pour exporter un jeu de donnée transformé ?¨

#Question : est ce qu'on s'embête vraiment à croiser les méthodes de normalisation avec les méthodes DEG pour la RNAseq ? (du moins pas maintenant je pense ?)
