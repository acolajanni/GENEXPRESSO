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


data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)

data_expr <- read.table("~/GIT/CPRD/DATA/TEST_COEXPR/data_expression.csv", 
                            row.names = 1, quote="\"",header = T)

# la fonction cor() calcul les corrélation entre les colonnes et pas les lignes : 
# il faut transposer la matrice et pour cela il existe la fonction t()
data_exprT=as.data.frame(t(data_expr))
# calculer la corrélation entre toutes les colonnes d'une matrice
result_correlation=cor(data_exprT)

# on veut passer d'une matrice à une liste de paires: on va utiliser la fonction "melt" 
# qui est disponible dans la library  reshape2 du package reshape2: il faut donc l'installer
list_correlation=melt(result_correlation)
# Matrice symétrique : supprimer les paires redondantes (ex : geneA - geneB et geneB - geneA)
# il faut utiliser la fonction as.charater() pour comparer des chaines
filtre = as.character(list_correlation[,1])<as.character(list_correlation[,2])
interaction = list_correlation[filtre,]

# Calcul des pvalue de corrélation pour le graphe de coexpression
data_expr = as.matrix(data_expr)
Cor_spe <- NULL
Cor_ken <- NULL
for (i in 1:nrow(interaction)){
  A = toString(interaction[["Var1"]][i])
  A = as.vector(data_expr[A,])
  
  B = toString(interaction[["Var2"]][i])
  B = as.vector(data_expr[B,])
  
  test_spe = cor.test(A,B, method = "spearman")
  Cor_spe = c(Cor_spe, test_spe$p.value)
  
  test_ken = cor.test(A,B, method = "kendall")
  Cor_ken = c(Cor_ken, test_ken$p.value)
}

interaction = cbind(interaction,PVal.Kendall = Cor_ken, Pval.Spearman = Cor_spe )


Matrix = copy(result_correlation)
for (col in 1:ncol(Matrix)){
  for (row in 1:nrow(Matrix)){
    if (Matrix[col,row] >= 0.75 || Matrix[col,row] <= -0.75){
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
     edge.width = 5,
     vertex.size = 10,
     vertex.color = "white",
     label.color = "black",
     label.font = 2)



