# Nom               : functionwilcoxtest.R
# Type              : Programme
# Objet             : Fonction pour Wilcoxon sur données microarrays et nanostring
# Input             : donnees microarrays et nanostring
# Output            : 
# Auteur            : Juliette CASEMAJOR
# R version         : 3.6
# Date de creation  : 16.04.2021
#______________________________________________________________________________

wilcoxDEG <- function(data, n1, n2, type){
  if (type == "nanostring") {
    for (i in nrow(data)){
      x = data[i,1:n1]
      y = data[i,n1+1:n1+n2]
    }
  }
}

