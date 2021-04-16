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

