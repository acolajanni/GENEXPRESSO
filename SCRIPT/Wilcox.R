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

#Creation des dataframes vides
wilcoxtsmicro <- data.frame()
wilcoxlmicro <- data.frame()
wilcoxgmicro <- data.frame()

#Test de wilcoxon two-sided effectué sur chaque ligne du dataframe
#Chaque valeur de pvalue est extraite puis rajoutée à un tableau de resultats
for (i in 1:nrow(micro)) {
  #Selection du 1er groupe de patients
  x = micro[i,1:6]
  #Selection du deuxieme groupe de patients
  y = micro[i,7:12]
  #Test
  test <- wilcox.test(as.numeric(x),as.numeric(y))
  #On recupere la pvalue
  pval <- test$p.value
  #On l'ajoute au tableau
  wilcoxtsmicro <- rbind(wilcoxtsmicro,pval)
}

#Test de wilcoxon less effectue sur chaque ligne du dataframe :
#On fait l'hypothese que l'expression du groupe 1 est inferieure a celle du groupe 2
for (i in 1:nrow(micro)) {
  x = micro[i,1:6]
  y = micro[i,7:12]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="less")
  pval <- test$p.value
  wilcoxlmicro <- rbind(wilcoxlmicro,pval)
}

#Test de wilcoxon greater effectue sur chaque ligne du dataframe :
#On fait l'hypothese que l'expression du groupe 1 est superieure a celle du groupe 2
for (i in 1:nrow(micro)) {
  x = micro[i,1:6]
  y = micro[i,7:12]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="greater")
  pval <- test$p.value
  wilcoxgmicro <- rbind(wilcoxgmicro,pval)
}

#On fusionne les tableaux par les colonnes
wilcoxmicro <- cbind(wilcoxtsmicro,wilcoxlmicro,wilcoxgmicro)

#Changement du nom des colonnes selon le type de test
colnames(wilcoxmicro)[1] <- "wilcox two-sided"
colnames(wilcoxmicro)[2] <- "wilcox less"
colnames(wilcoxmicro)[3] <- "wilcox greater"

########## NANOSTRING
#Fichier importe pour normalisation de donnees Nanostring avec la methode Nappa
nano <- read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv",row.names = 1,header = T, sep = ',')

#On a 64 patients : phenotype 1 pour patients 1 à 42
#phenotype 2 pour patients 43 a 64
#Et 770 genes

#Creation des dataframes vides
wilcoxtsnano<- data.frame()
wilcoxlnano <- data.frame()
wilcoxgnano <- data.frame()

#Test de wilcoxon two-sided effectué sur chaque ligne du dataframe
#Chaque valeur de pvalue est extraite puis rajoutée à un tableau de resultats
for (i in 1:nrow(nano)) {
  #Selection du 1er groupe de patients
  x = nano[i,1:42]
  #Selection du deuxieme groupe de patients
  y = nano[i,43:64]
  #Test
  test <- wilcox.test(as.numeric(x),as.numeric(y))
  #On recupere la pvalue
  pval <- test$p.value
  #On l'ajoute au tableau
  wilcoxtsnano<- rbind(wilcoxtsnano,pval)
}

#Test de wilcoxon less effectue sur chaque ligne du dataframe :
#On fait l'hypothese que l'expression du groupe 1 est inferieure a celle du groupe 2
for (i in 1:nrow(nano)) {
  x = nano[i,1:42]
  y = nano[i,43:64]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "less")
  pval <- test$p.value
  wilcoxlnano<- rbind(wilcoxlnano,pval)
}

#Test de wilcoxon greater effectue sur chaque ligne du dataframe :
#On fait l'hypothese que l'expression du groupe 1 est superieure a celle du groupe 2
for (i in 1:nrow(nano)) {
  x = nano[i,1:42]
  y = nano[i,43:64]
  test <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "greater")
  pval <- test$p.value
  wilcoxgnano<- rbind(wilcoxgnano,pval)
}

#On fusionne les tableaux par les colonnes
wilcoxnanonappa <- cbind(wilcoxtsnano,wilcoxlnano,wilcoxgnano)

#Changement du nom des colonnes selon le type de test
colnames(wilcoxnanonappa)[1] <- "Nappa.wilcox two-sided"
colnames(wilcoxnanonappa)[2] <- "Nappa.wilcox less"
colnames(wilcoxnanonappa)[3] <- "Nappa.wilcox greater"

#Faire une fonction qui fait la même chose pour chaque tools utilise ??
