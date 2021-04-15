# Nom               : Wilcox.R
# Type              : Programme
# Objet             : Test de Wilcoxon sur données microarrays et nanostring
# Input             : donnees microarrays et nanostring
# Output            : 
# Auteur            : Juliette CASEMAJOR
# R version         : 3.6
# Date de creation  : 15.04.2021
#______________________________________________________________________________

#Importation des packages

########## MICROARRAYS

#Importation des donnees microarrays a partir du fichier de simulation
micro <- read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv",row.names = 1,header = T, sep = ',')
micro

#Les donnees n'ont pas besoin d'etre normalisees car la simulation genere des
#donnees normalisees automatiquement
#On peut passer à test de wilcoxon (Test U) : comparaison de deux echantillons independants
#On a 100 genes et deux groupes de patients (2x6)

#Creation d'un dataframe vide
wilcoxtsmicro <- data.frame()
wilcoxlmicro <- data.frame()
wilcoxgmicro <- data.frame()

#Test de wilcoxon effectué sur chaque ligne du dataframe
#Chaque valeur de pvalue est extraite puis rajoutée à un tableau de resultats
for (i in 1:nrow(micro)) {
  x = micro[i,1:6]
  y = micro[i,7:12]
  test <- wilcox.test(as.numeric(x),as.numeric(y))
  pval <- test$p.value
  wilcoxtsmicro <- rbind(wilcoxtsmicro,pval)
}

for (i in 1:nrow(micro)) {
  x = micro[i,1:6]
  y = micro[i,7:12]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="less")
  pval <- test$p.value
  wilcoxlmicro <- rbind(wilcoxlmicro,pval)
}

for (i in 1:nrow(micro)) {
  x = micro[i,1:6]
  y = micro[i,7:12]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="greater")
  pval <- test$p.value
  wilcoxgmicro <- rbind(wilcoxgmicro,pval)
}

wilcoxmicro <- cbind(wilcoxtsmicro,wilcoxlmicro,wilcoxgmicro)
#Changement du nom de la colonne en p-values
colnames(wilcoxmicro)[1] <- "wilcox two-sided"
colnames(wilcoxmicro)[2] <- "wilcox less"
colnames(wilcoxmicro)[3] <- "wilcox greater"


########## NANOSTRING

