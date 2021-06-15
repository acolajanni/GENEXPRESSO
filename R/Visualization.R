############################
###### Visualization #######
############################

library(factoextra)
library(FactoMineR)
library(UpSetR)
library(igraph)


#' Compute a PCA of the different tool used in normalization and DEG analysis step to discriminate each methods by their found p-values on each genes 
#'
#' @param data.to.comp dataframe containing pvalues of genes being differentially expressed, 
#' with tools used in row, and genes in columns
#'
#' @return Principal component analysis plot with cos2 values for each methods
#' 
#' @import "factoextra" "FactoMineR"
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute Pvalues for all the methods 
#' # res.DEG = tools.DEG.RNAseq.merge(Data)
#' # Plotting PCA on methods
#' # PCA_tools(res.DEG)
PCA_tools <- function(data.to.comp){
  # Compute a PCA
  res.pca <- PCA(data.to.comp,graph=F)
  # Display the PCA plot with different parameters
  fviz_pca_ind (res.pca, col.ind = "cos2",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,
                ggrepel.max.overlaps = 3)
}


#' Compute a usable matrix to make an Upset plot
#' 
#' Fill the dataframe with binary values. 0 for non differentially expressed genes and 1 for the other. 
#'
#' @param data.to.comp Dataframe of DEG pvalues with genes in columns, and methods in rows.
#' @param threshold Integer. 
#' For the whole dataframe, if value <= threshold, value = 1, else value = 0.
#'
#' @return dataframe filled with binary values
#' 
#' @import "UpSetR"
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute Pvalues for all the methods 
#' # res.DEG = tools.DEG.RNAseq.merge(Data)
#' # Make a binary matrix to construct an Upset plot
#' # Upset = Upset.Binary.Dataframe(res.DEG)
#' #upset(Upset, sets = names(Upset), 
#' #          sets.bar.color = "#56B4E9", 
#' #          order.by = "freq", 
#' #          empty.intersections = NULL )
Upset.Binary.Dataframe <- function(data.to.comp, threshold){
  if(missing(threshold)){
    threshold = 0.05
  }
  # Replacing values of data.to.comp
  # 1 if below the threshold (significantly differentially expressed), else 0
  Upset = ifelse(data.to.comp <= threshold, 1, 0)
  Upset = as.data.frame(t(Upset))
  
  return(Upset)
}

#' Produce an igraph network with different colors to compare two networks with igraph
#'
#' For two graphs, red edges will be the common ones. It is used to compare two methods with one another, or two experimental conditions of a same set of genes for example. 
#'
#' @param g1 igraph class object.
#' @param g2 igraph class object.
#' @param g1.name Character string to be displayed on the legend of the graph for the first graph.
#' @param g2.name Character string to be displayed on the legend of the graph for the second graph.
#' @param display Logical value.
#' @param color.g1 Character string for a color to be displayed onto the g1 edges.
#' By default the color is "blue".
#' @param color.g2 Character string for a color to be displayed onto the g2 edges.
#' By default the color is "darkgreen".
#' 
#' TRUE if you want to to display the graph, otherwise it will be stocked as an object.
#' @return igraph class object with colored edges.
#' 
#' @import "igraph"
#' @export
#'
#' @examples
#' # Creating two datasets
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' # Computing relation network
#' #G.spearman = Make.df.graph(Data, cor.threshold = 0.25, Pvalue.threshold = F, 
#' #                                                      method = "spearman")
#' #G.kendall = Make.df.graph(Data, cor.threshold = 0.25, Pvalue.threshold = F, 
#' #                                                  method = "kendall")
#' # Comparing both relation network
#' #G.comp = relations.comparison(g1 = G.spearman, g2 = G.kendall,  
#' #                              g1.name = "Spearman", 
#' #                              g2.name = "kendall", 
#' #                              display = T)
relations.comparison <- function(g1,g2,g1.name,g2.name, display, color.g1, color.g2){
  # by default parameters
  if (missing(g1.name)){
    g1.name = "g1"
  }
  if (missing(g2.name)){
    g1.name = "g2"
  }
  if(missing(display)){
    display = TRUE
  }
  if(missing(color.g1)){
    color.g1 = "blue"
  }
  if(missing(color.g2)){
    color.g2 = "darkgreen"
  }
  # finding g1 exclusive edges
  diffg1 <- graph.difference(g1, g2)
  # finding g2 exclusive edges
  diffg2 <- graph.difference(g2, g1)
  # finding common edges
  interg1 <- graph.intersection(g1,g2, keep.all.vertices = T)
  # converting igraph class object in dataframe to add a "color" associated column 
  G1_edges = as_data_frame(diffg1,what = "edges")
  # by default, G1 exclusive edes, will be blue
  if (nrow(G1_edges) != 0){
    G1_edges$color = color.g1
  }
  # G2, exclusive will be green
  G2_edges = as_data_frame(diffg2,what = "edges")
  if (nrow(G2_edges) != 0){
    G2_edges$color = color.g2
  }
  
  # Common edges will be red
  Common_edges = as_data_frame(interg1, what="edges")
  if (nrow(Common_edges) != 0){
    
    Common_edges = Common_edges[-c(3,4)]
    Common_edges$color = "red"
  }
  # merging the three graph 
  # This object is equivalent to the union of both graph, but with specified colored edges 
  graph = rbind(G1_edges,G2_edges, Common_edges)
  # Converting back the dataframe into an igraph class object
  g = graph.data.frame(graph, directed = F)
  # Display or not the graph in the plot window.
  if(display){
    plot(g, #layout = lay,
         edge.width = 2,
         vertex.size = 2,
         vertex.label = NA,
         vertex.color = "white",
         label.color = "black",
         label.font = 2)
    
    inter = paste(g1.name, "and", g2.name, "intersection")
    
    legend(#"bottomleft",
      #"bottom",
      x=-1.5, y=-1.1,
      c(g1.name,g2.name,inter), 
      pch=18, 
      col=c(color.g1,color.g2,"red"), 
      pt.cex=0, #size of dots in the legend
      cex=1.2, #font size
      lty=c(1,1,1),
      lwd = 3,
      bty="n", #no frame around the legend
      ncol=1)
  }
  return(g)
}



#' Produce an Heatmap with pheatmap
#'
#' The dataset is subsetted to only keep DE genes in rows. Then, an heatmap from those genes is produced on the two groups within the whole sample
#' 
#' @param DEG 
#' Vector of SYMBOLs of DE genes
#' @param dataset 
#' Dataframe with sample in columns and genes in rows
#' @param col.group1 
#' colour for the annotation of the first group
#' @param col.group2 
#' colour for the annotation of the second group
#' @param n.breaks 
#' Number of breaks for the colour palette of the heatmap
#' @param is.log2 
#' Logical value. TRUE if the dataset is on log2 scale
#' @param sample.condition
#' Vector of the sample condition in the same order as sample names
#' @param samples 
#' Vector of the samples in the same order as sample conditions
#' @param main
#' Title of the heatmap
#' 
#' @return
#' @import "pheatmap"
#' @export
#'
#' @examples
#'# Import the dataset
#'Data = matrix(runif(5000, 10, 100), ncol=20)
#'group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#'genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#'colnames(Data) = group
#'row.names(Data) = genes 
#'
#'# Get a set of particular genes
#'DEG = genes <- paste0(rep(LETTERS[1:25], each=5), rep(c(1:10),each = 1))
#'
#'# Get the samples 
#'sample.condition = c(rep("control", times = 10), rep("case", times = 10) )
#'samples = colnames(Data)
#'
#'DEG.pheatmap(DEG, Data, is.log = FALSE, sample.condition = sample.condition, samples = samples )
#'
DEG.pheatmap = function(DEG, dataset, col.group1 = "cyan", col.group2 = "deeppink4", n.breaks, is.log2, sample.condition,samples,main ){
  # By default values
  if(missing(col.group1)){
    col.group1 = "deeppink4"
  }
  if(missing(col.group2)){
    col.group2 = "cyan"
  } 
  
  if(missing(n.breaks)){
    n.breaks = 8
  }
  
  # Subsetting the dataset to keep only the DE genes in rows
  dataset = dataset[ row.names(dataset)%in%DEG , ]
  
  
  # log2 scale allows a nicer visualization on heatmaps
  if(!is.log2){
    dataset = log2(dataset)
  }
  
  # Centering by the median
  MedianCentering <- function(x){
    (x - median(x))
  }
  dataset = t(apply(dataset, 1, MedianCentering))
  
  
  # Creating sample annotation data frame
  sample.DF = data.frame("sample" = sample.condition, row.names = samples)
  
  # Configuring colors
  annotation_colors = list(sample = c(sample1 = col.group2, sample2 = col.group1))
  names(annotation_colors$sample) = levels(sample.DF$sample)
  
  heatmap = pheatmap(dataset,
           annotation_col = sample.DF,
           annotation_colors = annotation_colors,
           cutree_cols = 2,
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           show_rownames = FALSE
           , inferno(n.breaks)
           , drop_levels = TRUE
           , main = main
  )
  
}
