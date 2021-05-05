###########################
########    DEG   #########
###########################


library(limma)
source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))


#' Creates a Model matrix
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row
#' @param cond1 Name of the first experimental condition
#' @param cond2 Name of the second experimental condition
#' @param ncond1 Number of sample in the first experimental condition
#' @param ncond2 Number of sample in the the second experimental condition
#'
#' @return Model matrix
#' @export
#'
#' @examples
#' 
#' ajouter un exemple
#' 
make_designMatrix <- function(dataset,cond1 = "A", cond2 = "B",ncond1=(ncol(dataset)/2),ncond2=(ncol(dataset)/2)){
  status.control = rep(cond1,ncond1)
  status.test = rep(cond2,ncond2)
  status = c(status.control,status.test)
  design = model.matrix(~0+status)
  colnames(design) <- c(cond1,cond2)
  return(design)
}  

#' Limma Function to analyse Differentially expressed genes
#'
#' Uses the limma function from the limma package to standardize the returned object
#'
#' @param dataset ataframe of expression values with samples in columns and genes in row
#' @param cond1 Name of the first experimental condition
#' @param cond2 Name of the second experimental condition
#' @param ncond1 Number of sample in the first experimental condition
#' @param ncond2 Number of sample in the the second experimental condition
#'
#' @return 
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL)
#' @export
#'
#' @examples
#' 
#' 
#' 
DEG_limma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  fit <- lmFit(dataset,design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  res.diff <- topTable(fit2, coef="diff",genelist=row.names(dataset), number=Inf)
  res.diff_limma <- data.frame(FoldChange = (res.diff$logFC) ,PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
  
  colnames(res.diff_limma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_limma)
}

#' Title
#'
#' @param dataset 
#' @param cond1 
#' @param cond2 
#' @param ncond1 
#' @param ncond2 
#'
#' @return
#' @export
#'
#' @examples
DEG_GEOlimma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  
  fit <- lmFit(dataset,design)
  fit2  <- contrasts.fit(fit, cm)
  load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  de <- topTable(fit22, number = nrow(dataset))
  res.diff_geolimma <- data.frame(FoldChange = (de$logFC), PValue=(de$adj.P.Val),genes=row.names(de))
  
  colnames(res.diff_geolimma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_geolimma)
}

#' Title
#'
#' @param dataset 
#' @param cond1 
#' @param cond2 
#' @param ncond1 
#' @param ncond2 
#'
#' @return
#' @export
#'
#' @examples
DEG_alternative <- function(res.diff_limma){
  
  res.up = copy(res.diff_limma)
  res.down = copy(res.diff_limma)
  for (row in 1:nrow(res.diff_limma)){
    if (res.diff_limma[['logFC']][row] > 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.down[['PValue']][row] = 1 
      
    }else if (res.diff_limma[['logFC']][row] < 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.up[['PValue']][row] = 1
    }
  }
  res = data.frame(PValue_Up = res.up$PValue, Pvalue_Down = res.down$PValue, SYMBOL = res.diff_limma$SYMBOL)
  return(res)
}

#' Title
#'
#' @param dataset 
#' @param cond1 
#' @param cond2 
#' @param ncond1 
#' @param ncond2 
#'
#' @return
#' @export
#'
#' @examples
wilcoxDEG <- function(data, n1, n2){
  wilcox <- data.frame()
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test1 <- wilcox.test(as.numeric(x),as.numeric(y))
    test2 <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "less")
    test3 <- wilcox.test(as.numeric(x),as.numeric(y), alternative = "greater")
    pval <- c(test1$p.value,test2$p.value,test3$p.value)
    wilcox <- rbind(wilcox,pval)
  }
  colnames(wilcox)[1] <- "wilcox two-sided"
  colnames(wilcox)[2] <- "wilcox less"
  colnames(wilcox)[3] <- "wilcox greater"
  row.names(wilcox) = row.names(data)
  return(wilcox)
}