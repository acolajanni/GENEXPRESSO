###########################
###### CoExpression #######
###########################

library(reshape2)
library(WGCNA)
library(igraph)


#' Computing a square matrix of correlation coefficients
#'
#' @param data Dataset of gene expression levels with genes in row and samples in columns
#' @param method Methods to compute the correlation coefficients.
#' "spearman" compute the Spearman's rho, 
#' "kendall" uses the Kendall's tau and
#' "pearson" the Pearson's product moment correlation coefficient.
#' These functions are called via the cor() function in the stats package.
#'
#' @return square matrix of correlation coefficients
#' @export
#'
#' @examples
#' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' 
#' # Computing square matrix with Spearman's rho
#' Cor = Cor.square.matrix(df,"spearman")
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

#' Computing a square matrix of Topological Overlap
#' 
#' Uses the TOMsimilarityFromExpr() function from the WGCNA package
#'
#' @param data Dataset of gene expression levels with genes in row and samples in columns
#' 
#' @return square matrix of TOM similarity measurements for each genes
#' 
#' @import "WGCNA"
#' @export
#'
#' @examples
#' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' 
#' # Computing TOM similarity
#' TOM = TOM.square.matrix(df)
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

#' Computing an adjacency table
#'
#' This function is similar to Cor.square.matrix() and TOM.square.matrix() as it converts the square matrix into a 
#' 3 columns dataframe with no repeated pair of variables. 
#'
#' @param data 
#' Dataset of gene expression levels with genes in row and samples in columns
#' 
#' @param method Methods to compute the correlation coefficients.
#' "spearman" compute the Spearman's rho, 
#' "kendall" uses the Kendall's tau and
#' "pearson" the Pearson's product moment correlation coefficient.
#' These functions are called via the cor() function in the stats package.
#' "TOM" uses the TOMsimilarityFromExpr() function from the WGCNA package. 
#' 
#' @return 
#' Datraframe of three columns. Two first contains pair of genes, and the third one the correlation coefficient or the TOM similarity.
#' 
#' @import "reshape2"
#' @export
#'
#' @examples
#' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' 
#' # computing correlation for pairs of genes with Spearman's rho
#' Adj = Make.adjacency.table(df,method = "spearman")
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

#' Computing an adjacency table with Pvalues associated to a coefficient of correlation 
#'
#' This function is similar to Cor.square.matrix() and TOM.square.matrix() as it converts the square matrix into a 
#' 4 columns dataframe with no repeated pair of variables and Pvalue associated with a correlation coefficient.
#'
#' @param data 
#' Dataset of gene expression levels with genes in row and samples in columns
#' @param Fast Logical value. 
#' if Fast = TRUE, only the pearson correlation coefficients will be measured with a fast one-step method. It uses the corAndPvalue() function from the WGCNA package.
#' if Fast = FALSEALSE, the pvalue of the correlation of the gene pair will be computed. 
#' Carefull, if Fast = F, it may take a while depending on number of genes in the dataset.
#'  
#' 
#' @param method Methods to compute the correlation coefficients.
#' "spearman" compute the Spearman's rho, 
#' "kendall" uses the Kendall's tau and
#' "pearson" the Pearson's product moment correlation coefficient.
#' These functions are called via the cor() function in the stats package.
#'
#' @return 
#' Dataframe with 4 columns. the two first are the pairs of genes,
#' the other two are the correlation coefficients and the associated pvalue.
#'
#' @import "reshape2"
#' @export
#'
#' @examples
#' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' 
#' # computing fast correlation and pvalues 
#' Adj = Make.adjacencyPVal(df,Fast = FALSE ,method = "spearman")
Make.adjacencyPVal <-function(data, Fast = FALSE, method){
  if (Fast == F){
    #First, we retrive the adjacency table depending on the method
        tools.coexpr <- switch(method,
                           spearman = {
                             interaction = Make.adjacency.table(data, method = "spearman")
                           },
                           
                           kendall = {
                             interaction = Make.adjacency.table(data, method = "kendall")
                           },
                           
                           stop("Enter something that switches me!")
                           
    )
    col.name = paste0("Pval.",method)
    cor = NULL
    # Then the pvalue is computed for each pair of gene
    for (i in 1:nrow(interaction)){
      # Retrieving first gene name, and the associated values
      A = toString(interaction[["Var1"]][i])
      A = as.vector(data[A,])
      # Retrieving second gene name and the associated values
      B = toString(interaction[["Var2"]][i])
      B = as.vector(data[B,])
      # Testing the correlation between two genes
      test = cor.test(as.numeric(A),as.numeric(B), method = method)
      cor = c(cor, test$p.value)
      
    }
    # Adding a Pvalue column that corresponds to the correlation coefficients
    interaction = cbind(interaction, Pval = cor )
    colnames(interaction)[colnames(interaction) == "Pval"] = col.name 
  } 
  # For fast computation time, pearson pvalue is used via corAndPvalue() function
  else if (Fast == T){
    data = as.matrix(t(data))
    # Computing the pvalue and associated correlations values
    Cor = corAndPvalue(data) 
    # retrieving pvalues and corelation coefficients
    Pval = melt(Cor$p)
    Cor = melt(Cor$cor)
    interaction = cbind(Cor,PValue = Pval$value)
    colnames(interaction) = c("Var1","Var2","cor.pearson","PVal.pearson")
    # remove repetitions of pairs
    filtre = as.character(interaction[,1])<as.character(interaction[,2])
    filtre2 = as.character(interaction[,1])<as.character(interaction[,2])
    interaction = interaction[filtre,]
  }
  
  return(interaction)
}

#' Computing an adjacency table with or without Pvalues, with Spearman, Kendall correlation coefficient and TOM similarity.
#'
#' Calling Make.adjacencyPVal() with all the implemented methods. 
#'
#' @param data 
#' Dataset of gene expression levels with genes in row and samples in columns
#' @param PValue Logical value. 
#' if PValue = TRUE, all the pvalues associated with the corresponding coeffient of correlation will be computed. Could be very long on large dataframe
#' if PValue = FALSE, only the coefficient of correlation, and TOM similarity will be computed.
#' @return
#' Dataframe with several columns. the two first are the pairs of genes,
#' the other are the coefficient of correlation depending on the used method, and maybe the PValues depending on the argument.
#' @export
#'
#' @examples
#' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' 
#' # computing fast correlation and pvalues 
#' Adj = Make.full.adjacency(df,PValue = TRUE)
Make.full.adjacency <- function(data, PValue = T){
  
  if (PValue == T){
    # Computing all the Pvalues (and correlation coefficient) for all the methods
    interaction_spearman = Make.adjacencyPVal(data, Fast = F, method = "spearman")
    interaction_kendall = Make.adjacencyPVal(data, Fast = F, method = "kendall")
    interaction_Fast = Make.adjacencyPVal(data, Fast = T)
    # Merging them to have all the pvalues and coefficients in the same dataframe
    interaction_Cor = merge(interaction_spearman,interaction_kendall,by = c("Var1","Var2"))
    interaction_Cor = merge(interaction_Cor,interaction_Fast,by = c("Var1","Var2"))                        
  }
  else {
    # Since PValue = F, we only compute all the correlation coefficient for all the methods
    interaction_spearman =  Make.adjacency.table(data, method = "spearman")
    interaction_kendall =  Make.adjacency.table(data, method = "kendall")
    # Merging them to have all the corrleation coefficients in the same dataframe
    interaction_Cor = merge(interaction_spearman,interaction_kendall,by = c("Var1","Var2"))
  }
  # TOM similarity has no pvalue to test their signifiance
  interaction_TOM = Make.adjacency.table(data, method = "TOM")
  # Adding the TOM similarity to the dataframe
  interaction = merge(interaction_Cor, interaction_TOM, by = c("Var1","Var2"))
  
  return(interaction)
}


#' Computing a dataframe of usable format for Network plotting with iGraph 
#' 
#' Remove all the pairs that has a correlation coefficient below the cor.threshold, or the Pvalue.threshold given in argument.
#' Only the pairs that are sufficiently correlated or similar enough will be kept to produce an object of igraph class.
#'
#' @param data 
#' Dataset of gene expression levels with genes in row and samples in columns.
#' @param cor.threshold 
#' Threshold to apply for minimal correlation or similarity to keep the pair of gene.
#' @param Pvalue.threshold Logical value.
#' If TRUE, all the pairs not significantly correlated will be removed. Could take a while since all the pvalues needs to be computed.
#' If FALSE, only a thresholding through the correlation coefficient will be apply.
#' @param method Methods to compute the correlation coefficients and the p-values (if Pvalue.threshold = T)
#' "spearman" compute the Spearman's rho, 
#' "kendall" uses the Kendall's tau and
#' These functions are called via the cor() function in the stats package.
#' "TOM" uses the TOMsimilarityFromExpr() function from the WGCNA package.
#'
#' @return igraph class object
#' @import "igraph"
#' @export
#'
#' @examples
#' #' # Creating a dataset
#' df = matrix(runif(500, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=1))
#' colnames(df) = group
#' row.names(df) = genes
#' # Compute the graph in an usable format for igraph
#' Graph = Make.df.graph(df, cor.threshold = 0.5, Pvalue.threshold = FALSE, method = "spearman")
#' # Plotting it with igraph
#' plot(Graph)
Make.df.graph<-function(data, cor.threshold, Pvalue.threshold = FALSE, method ){
  
  if (Pvalue.threshold == FALSE){
    # In this case, no need to compute pvalues
    tools.graph <- switch(method,
                          spearman = {
                            # Only correlation coefficients are produced
                            relations = Make.adjacency.table(data, method = "spearman")
                            cor = 'cor.spearman' 
                          },
                          
                          kendall = {
                            relations = Make.adjacency.table(data, method = "kendall")
                            cor = "cor.kendall"
                          },
                          
                          TOM =  {
                            relations = Make.adjacency.table(data, method = "TOM")
                            cor = "TOM.similarity"
                          },
                          
                          stop("Enter something that switches me!")
    )
  }
  else { 
    # In this case pvalues and correlation coefficients are computed with Make.adjacencyPval()
    tools.graph <- switch(method,
                          spearman = {
                            relations = Make.adjacencyPVal(data, method = "spearman")
                            cor = 'cor.spearman'
                            pvalue = "Pval.spearman"
                          },
                          
                          kendall = {
                            relations = Make.adjacencyPVal(data, method = "kendall")
                            cor = "cor.kendall"
                            pvalue = "Pval.kendall"
                          },
                          
                          TOM =  {
                            relations = Make.adjacency.table(data, method = "TOM")
                            cor = "TOM.similarity"
                            # No pvalues could be associated with TOM similarity
                            Pvalue.threshold = FALSE
                          },
                          
                          stop("Enter something that switches me!")
    )
  }  
  # TOM similarity varies between 0 and 1, whereas coefficients varies between -1 and 1.
  # We are not interested in the sign of the correlation, so absolute values are used only for corrlation coefficients
  if (method == "TOM"){
    df.graph = subset(relations, relations[[cor]] >= cor.threshold)
  }
  else {
    relations[[cor]] = abs(relations[[cor]])
  }
  # When the pvalues are used the threshold is 0.05
  # subset is used to only kept the part of a dataframe that respect the condition
  if (Pvalue.threshold == T){
    df.graph = subset(relations, (relations[[cor]] >= cor.threshold) & (relations[[pvalue]] <= 0.05) )
  }
  # When we don't want to look at pvalues, only correlation coefficients or TOM similarity are used.
  else {
    df.graph = subset(relations, relations[[cor]] >= cor.threshold)
  }
  
  # To compute igraph class object, only the pair of variables are needed
  # This way, we remove all the coefficients and Pvalues associated columns
  df.graph = df.graph[-c(3:ncol(df.graph))]
  
  # If the threshold is too high, an error is raised
  # A too high threshold produce a dataframe with no rows (no pair of gene is correlated enough)
  if (nrow(df.graph) != 0){
    df.graph = graph.data.frame(df.graph, directed = FALSE)
  }
  else{
    message = paste("No variable have been found having a", cor, "this high")
    stop(message)
  }
  
  return(df.graph)
}