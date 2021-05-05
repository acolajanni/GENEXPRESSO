###########################
########    DEG   #########
###########################


library(limma)
source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))


#' Creates a Model matrix.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param cond1 Name of the first experimental condition.
#' @param cond2 Name of the second experimental condition.
#' @param ncond1 Number of sample in the first experimental condition.
#' @param ncond2 Number of sample in the the second experimental condition.
#'
#' @return Model matrix.
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

#' Differentially expressed genes analysis through the limma package.
#'
#' Uses the lm.fit(),eBayes() and topTable() from the limma package to automate the steps of the analysis.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param design vector of 0 and 1 values. 0 for the first experimental condition, 1 for the second one.
#' @param contrast.matrix design matrix like one produced by the make_designMatrix() function.
#'
#' @return
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL)
#' @export
#'
#' @examples
DEG_limma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  # To fit the model to the data, it needs the model matrix
  fit <- lmFit(dataset,design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  # Classifying genes through their pvalue being differentially expressed
  res.diff <- topTable(fit2, coef="diff",genelist=row.names(dataset), number=Inf)
  res.diff_limma <- data.frame(FoldChange = (res.diff$logFC) ,PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
  
  colnames(res.diff_limma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_limma)
}

#' Differentially expressed genes analysis through the limma package, with the GEOlimma method.
#'
#' Uses the lm.fit(),eBayes() and topTable() from the limma package to automate the steps of the analysis.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param design vector of 0 and 1 values. 0 for the first experimental condition, 1 for the second one.
#' @param contrast.matrix design matrix like one produced by the make_designMatrix() function.
#' 
#' @return
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL) 
#' @export
#' @examples
#' 
#' 
DEG_GEOlimma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  # To fit the model to the data, it needs the model matrix
  fit <- lmFit(dataset,design)
  fit2  <- contrasts.fit(fit, cm)
  # GEOlimma_probabilities.rda is a file that contains real probabilities of a gene being differentially expressed
  load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
  # Those data are used to help the model to fit
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  # Classifying genes through their pvalue being differentially expressed
  de <- topTable(fit22, number = nrow(dataset))
  res.diff_geolimma <- data.frame(FoldChange = (de$logFC), PValue=(de$adj.P.Val),genes=row.names(de))
  
  colnames(res.diff_geolimma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_geolimma)
}

#' Takes the results of DEG_GEOlimma and DEG_limma to extract the pvalues of alternative hypothesis.
#'
#' @param res.diff_limma 
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue),
#' the last column corresponds to the gene names (SYMBOL).
#'
#' @return
#' A dataframe with 3 columns is returned. The two columns corresponds to the pvalue
#' of a gene being respectively, more expressed and less expressed in the second condition
#' @export
#'
#' @examples
#' #get the results of limma analysis
#' res.diff = DEG_limma(dataset,design, contrast.matrix = cm)
#' #compute the results of alternative hypothesis
#' res.diff.alternative = DEG_alternative(res.diff)
DEG_alternative <- function(res.diff_limma){
  # two dataframes are created from the original given in argument
  res.up = copy(res.diff_limma)
  res.down = copy(res.diff_limma)
  # If the pvalue is significant, and depending on the sign of the logFoldChange, we can set the probability to one for one hypothesis
  # The idea is that, if the fist condition has significantly higher expression values than the second one, the gene can't be more expressed in the second condition, and vice versa.
  for (row in 1:nrow(res.diff_limma)){
    if (res.diff_limma[['logFC']][row] > 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.down[['PValue']][row] = 1 
      
    }else if (res.diff_limma[['logFC']][row] < 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.up[['PValue']][row] = 1
    }
  }
  # This way, two columns are created, one for the probability of a gene being more expressed in the second condition (Pvalue_Up), and less expressed (Pvalue_Down)
  res = data.frame(PValue_Up = res.up$PValue, Pvalue_Down = res.down$PValue, SYMBOL = res.diff_limma$SYMBOL)
  return(res)
}


#' Multiple Wilcoxons tests for each row of a dataframe
#'
#' Testing, for one row at the time, if the first series of values are different, greater or less than the values of the second condition.
#'
#' @param data dataframe of gene expression levels: Gene names in rows, samples in columns.
#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return
#' Dataframe with three columns, each corresponding to a hypothesis : 
#' weak hypothesis, and two for the strong hypothesis "alternative = "less" ", and "greater"
#' For the second condition has values that are less or greater than values of the first one.
#' @export
#'
#' @examples
wilcoxDEG <- function(data, n1, n2){
  wilcox <- data.frame()
  # for each row in the dataset, three Wilconxon test are made, depending on the hypothesis.
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test1 <- wilcox.test(as.numeric(x),as.numeric(y))
    test2 <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "less")
    test3 <- wilcox.test(as.numeric(x),as.numeric(y), alternative = "greater")
    pval <- c(test1$p.value,test2$p.value,test3$p.value)
    wilcox <- rbind(wilcox,pval)
  }
  # The dataframe named "wilcox" contains 3 columns  that contains the pvalue for the three possible hypothesis
  colnames(wilcox)[1] <- "wilcox two-sided"
  colnames(wilcox)[2] <- "wilcox less"
  colnames(wilcox)[3] <- "wilcox greater"
  row.names(wilcox) = row.names(data)
  return(wilcox)
}