# Nom               : coexpr.R
# Type              : Programme
# Objet             : Analyse coexpression Pearson + WGCNA
# Input             : Dataset (microarray, RNA-seq, Nanostring)
# Output            : iplot
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 16.04.2021
#______________________________________________________________________________
library("reshape2")
#install.packages('igraph')
library('igraph')
library(data.table)

data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)

data_expr <- read.table("~/GIT/CPRD/DATA/TEST_COEXPR/data_expression.csv", 
                            row.names = 1, quote="\"",header = T)

# la fonction cor() calcul les corrélation entre les colonnes et pas les lignes : 
# il faut transposer la matrice et pour cela il existe la fonction t()
data_exprT=t(data_expr)
# calculer la corrélation entre toutes les colonnes d'une matrice
result_correlation=cor(data_exprT)

# on veut passer d'une matrice à une liste de paires: on va utiliser la fonction "melt" 
# qui est disponible dans la library  reshape2 du package reshape2: il faut donc l'installer
list_correlation=melt(result_correlation)
# Matrice symétrique : supprimer les paires redondantes (ex : geneA - geneB et geneB - geneA)
# il faut utiliser la fonction as.charater() pour comparer des chaines
filtre = as.character(list_correlation[,1])<as.character(list_correlation[,2])
interaction = list_correlation[filtre,]

Matrix = copy(result_correlation)
for (col in 1:ncol(Matrix)){
  for (row in 1:nrow(Matrix)){
    if (Matrix[col,row] >= 0.8 || Matrix[col,row] <= -0.8){
      Matrix[col,row] = 1
    }else{
      Matrix[col,row] = 0
    }
  }
}
A = graph_from_adjacency_matrix(Matrix, mode='undirected', diag=F)
Isolated = which(degree(A)==0)
G2 = delete.vertices(A, Isolated)
plot(G2,
     edge.color = "black",
     edge.width = 3,
     frame = T)



