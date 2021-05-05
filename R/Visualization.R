############################
###### Visualization #######
############################

library(factoextra)
library(FactoMineR)
library(UpSetR)
library(igraph)
library(data.table)

#' Compute a PCA of the different tool used in normalization and DEG analysis step to discriminate each methods by their found p-values on each genes 
#'
#' @param data.to.comp datarame containing pvalues of genes being differentially expressed, 
#' with tools used in row, and genes in columns
#'
#' @return Principal component analysis plot with cos2 values for each methods
#' @export
#'
#' @examples
PCA_tools <- function(data.to.comp){
  # Compute a PCA
  res.pca <- PCA(data.to.comp,graph=F)
  # Display the PCA plot with different parameters
  fviz_pca_ind (res.pca, col.ind = "cos2",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,
                ggrepel.max.overlaps = 3)
  return(res.pca)
}


#' Compute a usable matrix to make an Upset plot
#' 
#' Filling the dataframe with binary values. 0 for non differentially expressed genes and 1 for the other. 
#'
#' @param data.to.comp Dataframe of DEG pvalues with genes in columns, and methods in rows.
#' @param threshold Threshold value to fill the dataframe with binary values. By default threshold = 0.05
#'
#' @return dataframe filled with binary values
#' 
#' @export
#'
#' @examples
Upset.Binary.Dataframe <- function(data.to.comp, threshold){
  if(missing(threshold)){
    threshold = 0.05
  }
  # à refaire avec la fonction ifelse()
  Upset <- copy(data.to.comp)
  for (i in names(Upset)){
    for (u in 1:nrow(Upset)){
      if (log == F){
        if(data.to.comp[[i]][u] > threshold){
          Upset[[i]][u] = 0}
        else{
          Upset[[i]][u] = 1
        }
      }
    }
  }
  Upset = as.data.frame(t(Upset))
  return(Upset)
}