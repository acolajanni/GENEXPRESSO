# Nom               : functionwilcoxtest.R
# Type              : Programme
# Objet             : Fonction pour Wilcoxon sur données microarrays et nanostring
# Input             : donnees microarrays et nanostring
# Output            : 
# Auteur            : Juliette CASEMAJOR
# R version         : 3.6
# Date de creation  : 16.04.2021
#______________________________________________________________________________

########### OPTION 1 
########### La fonction traite les deux types de donnees separement et renvoie
########### du two-sided pour nanostring ou du less/great pour microarrays en
########### fonction du type donne par l'utilisateur

wilcoxDEG <- function(data, n1, n2, type){
  if (type == "nanostring") {
    wilcoxnano <- data.frame()
    for (i in 1:nrow(data)){
      x = data[i,1:n1]
      y = data[i,(n1+1):(n1+n2)]
      test <- wilcox.test(as.numeric(x),as.numeric(y))
      pval <- test$p.value
      wilcoxnano<- rbind(wilcoxnano,pval)
    }
    colnames(wilcoxnano) <- "Wilcox two-sided"
    row.names(wilcoxnano) = row.names(data)
    return(wilcoxnano)
  }
  else {
    wilcoxlessmicro <- data.frame()
    for (i in 1:nrow(data)){
      x = data[i,1:n1]
      y = data[i,(n1+1):(n1+n2)]
      test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="less")
      pval <- test$p.value
      wilcoxlessmicro<- rbind(wilcoxlessmicro,pval)
    }
    wilcoxgreatermicro <- data.frame()
    for (i in 1:nrow(data)){
      x = data[i,1:n1]
      y = data[i,(n1+1):(n1+n2)]
      test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="greater")
      pval <- test$p.value
      wilcoxgreatermicro<- rbind(wilcoxgreatermicro,pval)
    }
    wilcoxmicro <- cbind(wilcoxlessmicro,wilcoxgreatermicro)
    colnames(wilcoxmicro)[1] <- "wilcox less"
    colnames(wilcoxmicro)[2] <- "wilcox greater"
    row.names(wilcoxmicro) = row.names(data)
    return(wilcoxmicro)
  }
}

########### OPTION 2
########### La fonction traite les deux types de donnees de la meme maniere et
########### donne un dataframe avec les 3 alternatives pour le test de Wilcox

wilcoxDEG4 <- function(data, n1, n2){
  wilcoxts <- data.frame()
  wilcoxless <- data.frame()
  wilcoxgreater <- data.frame()
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test <- wilcox.test(as.numeric(x),as.numeric(y))
    pval <- test$p.value
    wilcoxts<- rbind(wilcoxts,pval)
  }
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="less")
    pval <- test$p.value
    wilcoxless<- rbind(wilcoxless,pval)
  }
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test <- wilcox.test(as.numeric(x),as.numeric(y),alternative="greater")
    pval <- test$p.value
    wilcoxgreater<- rbind(wilcoxgreater,pval)
  }
  wilcox <- cbind(wilcoxts, wilcoxless, wilcoxgreater)
  colnames(wilcox)[1] <- "wilcox two-sided"
  colnames(wilcox)[2] <- "wilcox less"
  colnames(wilcox)[3] <- "wilcox greater"
  row.names(wilcox) = row.names(data)
  return(wilcox)
}

