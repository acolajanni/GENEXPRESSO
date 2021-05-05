###########################
###### CoExpression #######
###########################

library(reshape2)
library(WGCNA)


#' Compute a square matrix of correlation coefficients
#' 
#' 
#'
#' @param data Dataset of gene expression levels with genes in row and samples in columns
#' @param method Methods to compute the correlation coefficients.
#' "spearman" compute the Spearman's rho, "kendall", the Kendall's tau and
#' "pearson" the Pearson's product moment correlation coefficient.
#' These functions are called via the cor() function in the stats package.
#'
#' @return square matrix of correlation coefficients
#' @export
#'
#' @examples
#' # Creating a dataset
#' A = runif(50,5,100)
#' B = runif(50,5,100)
#' C = runif(50,5,100)
#' D = runif(50,5,100)
#' E = data.frame(A,B,C,D)
#' # Computing square matrix with Spearman's rho
#' Cor = Cor.square.matrix(E,"spearman")
Cor.square.matrix <- function(data, method){
  
  if (!method%in%c("pearson","kendall","spearman")){
    stop("Enter something that switches me!")
  }
  
  print(c(method, "square matrix"))
  data_expr=as.data.frame(t(data))
  # Computing the mutiple correlation coefficients with the methods in arguments
  result_correlation=cor(data_expr, method = method)
  return(result_correlation)
}

#' Compute a square matrix of Topological Overlap
#' 
#' Uses the TOMsimilarityFromExpr() function from the WGCNA package
#'
#' @param data Dataset of gene expression levels with genes in row and samples in columns
#' 
#' @return square matrix of TOM similarity measurements for each genes
#' @export
#'
#' @examples
#' # Creating a dataset
#' A = runif(50,5,100)
#' B = runif(50,5,100)
#' C = runif(50,5,100)
#' D = runif(50,5,100)
#' E = data.frame(A,B,C,D)
#' # Computing TOM similarity
#' TOM = TOM.square.matrix(E)
TOM.square.matrix <- function(data){
  # arrange the dataset in order to use the TOMsimilarityFromExpr() function without problem
  genes = row.names(data)
  data = t(as.data.frame(data))
  row.number = nrow(data)
  col.number = ncol(data)
  obs.number = col.number*row.number
  # Converting every values of the dataframe into numeric values
  data[1:obs.number] = sapply(data[1:obs.number], as.numeric) 
  TOM = TOMsimilarityFromExpr(data)
  # Needs to rename the columns and rows
  row.names(TOM) = genes
  colnames(TOM) = genes
  return(TOM)
}

#' Compute an adjacency table
#'
#' This function is similar to Cor.square.matrix() and TOM.square.matrix() as it converts the square matrix into a 
#' 3 columns dataframe with no repeated pair of variables. 
#'
#' @param data 
#' Dataset of gene expression levels with genes in row and samples in columns
#' 
#' @param method Methods to compute the correlation coefficients.
#' "spearman" compute the Spearman's rho, "kendall", the Kendall's tau and
#' "pearson" the Pearson's product moment correlation coefficient.
#' These functions are called via the cor() function in the stats package.
#' "TOM" uses the TOMsimilarityFromExpr() function from the WGCNA package. 
#' 
#' @return 
#' Datraframe of three columns. Two first contains pair of genes, and the third one the correlation coefficient or the TOM similarity.
#' @export
#'
#' @examples
#' #' # Creating a dataset
#' A = runif(50,5,100)
#' B = runif(50,5,100)
#' C = runif(50,5,100)
#' D = runif(50,5,100)
#' E = data.frame(A,B,C,D)
#' # computing correlation for pairs of genes with Spearman's rho
#' Adj = Make.adjacency.table(E,method = "spearman")
Make.adjacency.table <- function(data, method){
  
  tools.coexpr <- switch(method,
                         
                         kendall = {
                           result_correlation = Cor.square.matrix(data, method = "kendall")
                           method = paste0("cor.",method)
                           print("done")
                         },
                         
                         spearman = {
                           result_correlation = Cor.square.matrix(data, method = "spearman")
                           method = paste0("cor.",method)
                           print("done")
                         },
                         
                         TOM = {
                           result_correlation = TOM.square.matrix(data)
                           method = "TOM.similarity"
                         },
                         
                         stop("Enter something that switches me!")
  )
  
  list_correlation=melt(result_correlation)
  filtre = as.character(list_correlation[,1])<as.character(list_correlation[,2])
  interaction = list_correlation[filtre,]
  
  
  colnames(interaction) = c("Var1", "Var2", method)
  return(interaction)
}
