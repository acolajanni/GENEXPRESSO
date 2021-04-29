# Nom               : coexpr.R
# Type              : Programme
# Objet             : Analyse coexpression Pearson + WGCNA
# Input             : Dataset (microarray, RNA-seq, Nanostring)
# Output            : iplot
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 16.04.2021
#______________________________________________________________________________

#library("reshape2")
#install.packages('igraph')
#library('igraph')
#library(data.table)
source(file.path("./SCRIPT","functions.R"))

data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)


# la fonction cor() calcul les corrélation entre les colonnes et pas les lignes : 
# il faut transposer la matrice et pour cela il existe la fonction t()
data_exprT=as.data.frame(t(data_expr))
# calculer la corrélation entre toutes les colonnes d'une matrice
result_correlation=cor(data_exprT)
A = corAndPvalue(result_correlation) ## Package WGCNA

# on veut passer d'une matrice à une liste de paires: on va utiliser la fonction "melt" 
# qui est disponible dans la library  reshape2 du package reshape2: il faut donc l'installer
list_correlation=melt(result_correlation)

# Matrice symétrique : supprimer les paires redondantes (ex : geneA - geneB et geneB - geneA)
# il faut utiliser la fonction as.charater() pour comparer des chaines
filtre = as.character(list_correlation[,1])<as.character(list_correlation[,2])
filtre2 = as.character(A[,1])<as.character(A[,2])
interaction1 = list_correlation[filtre,]
interaction2 = A[filtre2,]

interaction = merge(interaction,interaction2, by = c("Var1","Var2"))
colnames(interaction) = c("Gene1","Gene2","Statistic","PValue")

###
data = as.matrix(t(data_expr))
Cor = corAndPvalue(data) 
Pval = melt(Cor$p)
Cor = melt(Cor$cor)
interaction = cbind(Cor,PValue = Pval$value)
colnames(interaction) = c("Gene1","Gene2","cor.level","PValue")

filtre = as.character(interaction[,1])<as.character(interaction[,2])
filtre2 = as.character(A[,1])<as.character(A[,2])
interaction = interaction[filtre,]
###


# Calcul des pvalue de corrélation pour le graphe de coexpression
data_expr = as.matrix(data_expr)
Cor_spe <- NULL
Cor_ken <- NULL
for (i in 1:nrow(interaction)){
  print(i)
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
    #print(Matrix[col,row])
    if (Matrix[col,row] >= 0.75 || Matrix[col,row] <= -0.75){
      Matrix[col,row] = 1
    }else{
      Matrix[col,row] = 0
    }
  }
}

# Affichage avec igraph
set.seed(12)
G = graph_from_adjacency_matrix(Matrix, mode='undirected', diag=F)
LO = layout_with_fr(G)

Isolated = which(degree(G)==0)
G2 = delete.vertices(G, Isolated)
LO2 = LO[-Isolated,]
plot(G2, layout=LO2,
     edge.color = "black",
     edge.width = 5,
     vertex.size = 10,
     vertex.color = "white",
     label.color = "black",
     label.font = 2)

############################################################################
# Fonctions
############################################################################

Cor.square.matrix <- function(data, method){
  
  if (missing(method)){
    method = "kendall"
  }
  
  if (!method%in%c("pearson","kendall","spearman")){
    stop("Enter something that switches me!")
  }
  
  
  data_expr=as.data.frame(t(data))
  result_correlation=cor(data_expr, method = method)
  return(result_correlation)
}

TOM.square.matrix <- function(data){
  data = t(as.data.frame(data))
  row.number = nrow(data)
  col.number = ncol(data)
  obs.number = col.number*row.number
  data[1:obs.number] = sapply(data[1:obs.number], as.numeric) 
  TOM = TOMsimilarityFromExpr(data)
  row.names(TOM) = row.names(RNA)
  colnames(TOM) = row.names(RNA)
  return(TOM)
}


Make.adjacency.table <- function(data, method){
  if (missing(method)){
    method = "kendall"
  }
  
  if (method == "TOM"){
    result_correlation = TOM.square.matrix(data)
    method = "TOM.similarity"
  }
  
  else{
    result_correlation = Cor.square.matrix(data, method)
    method = paste0("cor.",method)
  }

  list_correlation=melt(result_correlation)
  filtre = as.character(list_correlation[,1])<as.character(list_correlation[,2])
  interaction = list_correlation[filtre,]

  
  colnames(interaction) = c("Var1", "Var2", method)
  return(interaction)
}

Make.adjacencyPVal <-function(data, Fast = F){
  
  if (Fast == F){
    interaction2 = Make.adjacency.table(data, method = "spearman")
    interaction1 = Make.adjacency.table(data, method = "kendall")
    interaction = merge(interaction1,interaction2, by = c("Var1","Var2"))
    
    Cor_spe <- NULL
    Cor_ken <- NULL
    for (i in 1:nrow(interaction)){
      #print(i)
      A = toString(interaction[["Var1"]][i])
      A = as.vector(data[A,])

    
      B = toString(interaction[["Var2"]][i])
      B = as.vector(data[B,])
    
      test_spe = cor.test(as.numeric(A),as.numeric(B), method = "spearman")
      Cor_spe = c(Cor_spe, test_spe$p.value)
    
      test_ken = cor.test(as.numeric(A),as.numeric(B), method = "kendall")
      Cor_ken = c(Cor_ken, test_ken$p.value)
    }
    interaction = cbind(interaction,PVal.Kendall = Cor_ken, Pval.Spearman = Cor_spe )
  } 
  else if (Fast == T){
    data = as.matrix(t(data))
    Cor = corAndPvalue(data) 
    Pval = melt(Cor$p)
    Cor = melt(Cor$cor)
    interaction = cbind(Cor,PValue = Pval$value)
    colnames(interaction) = c("Var1","Var2","cor.pearson","PVal.Pearson")
    
    filtre = as.character(interaction[,1])<as.character(interaction[,2])
    filtre2 = as.character(interaction[,1])<as.character(interaction[,2])
    interaction = interaction[filtre,]
  }
  
  return(interaction)
}

Make.full.adjacency <- function(data, PValue = T){
  
  if (PValue == T){
    interaction_Cor = Make.adjacencyPVal(data, Fast = F)
    interaction_Fast = Make.adjacencyPVal(data, Fast = T)
    interaction_Cor = merge(interaction_Cor,interaction_Fast, by = c("Var1","Var2"))
  }
  else {
    interaction_spearman =  Make.adjacency.table(data, method = "spearman")
    interaction_kendall =  Make.adjacency.table(data, method = "kendall")
    interaction_pearson = Make.adjacency.table(data, method = "pearson")
    interaction_Cor = merge(interaction_spearman,interaction_kendall,by = c("Var1","Var2"))
    interaction_Cor = merge(interaction_Cor,interaction_pearson,by = c("Var1","Var2"))
    
  }
  
  interaction_TOM = Make.adjacency.table(data, method = "TOM")
  interaction = merge(interaction_Cor, interaction_TOM, by = c("Var1","Var2"))
  
  return(interaction)
}

Make.relation.matrix <- function(data){
  data = as.matrix(data)
  Matrix = Cor.square.matrix(data)
  for (col in 1:ncol(Matrix)){
    for (row in 1:nrow(Matrix)){
      if (is.na(Matrix[col,row])){
        Matrix[col,row] = 0
      }
      else if (Matrix[col,row] >= 0.75 || Matrix[col,row] <= -0.75){
        Matrix[col,row] = 1
      }else{
        Matrix[col,row] = 0
      }
    }
  }
  return(Matrix)
}

Make.df.graph<-function(data, threshold, method ){
  
  relations = Make.full.adjacency(data,PValue = F)

  if (method == "spearman"){
    cor = 'cor.spearman'
  }
  else if (method == "kendall"){
    cor = "cor.kendall"
  }
  else if (method == "pearson"){
    cor = "cor.pearson"
  }
  else if (method == "TOM"){
    cor = "TOM.similarity"
  }
  else {
    stop("Enter something that switches me!")
  }
  
  if (!method%in%c("TOM")){
    relations[[cor]] = abs(relations[[cor]])
  }
  
  var1 = NULL
  var2 = NULL
  
  for (row in nrow(relations):1){
    if (relations[[cor]][row] >= threshold) {
      var1 = c(var1, toString(relations[["Var1"]][row]))
      var2 = c(var2, toString(relations[["Var2"]][row]))
    }
  }
  df.graph = data.frame(Var1 = var1, Var2 = var2)
  df.graph = graph.data.frame(df.graph, directed = FALSE)
  return(df.graph)
}


Plot.relation.graph <-function(data){
  Matrix = Make.relation.matrix(data)
  set.seed(12)
  G = graph_from_adjacency_matrix(Matrix, mode='undirected', diag=F)
  LO = layout_with_fr(G)
  
  Isolated = which(degree(G)==0)
  G2 = delete.vertices(G, Isolated)
  LO2 = LO[-Isolated,]
  plot = plot(G2, 
       edge.color = "black",
       edge.width = 3,
       vertex.size = 10,
       vertex.color = "white",
       label.color = "black",
       label.font = 2)
  return(G2)
}

########### Appel de fonctions
RNA = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ10000x30.csv", header = TRUE,row.names = 1)
RNA = RNA[1:60,]


spearman = Make.df.graph(RNA,0.83,"spearman")
kendall = Make.df.graph(RNA,0.65,"kendall")
TOM = Make.df.graph(RNA,0.25,'TOM')
plot(spearman)
plot(kendall)
plot(TOM)

######################################################################
#Total_RNA = Make.adjacencyPVal(RNA)
Total_RNA = Make.adjacency.graph(RNA)

N = ncol(RNA)
RNA_control = RNA[1:(N/2)]
RNA_test = RNA[((N/2)+1):N]

par(mfrow=c(1,3))
# Plot total
RNA_tot = Make.adjacency.graph(RNA)
Plot.relation.graph(RNA)
#Plot gr1
RNA_C = Make.adjacency.graph(RNA_control)
A = Plot.relation.graph(RNA_control)
#Plot gr2
RNA_T = Make.adjacency.graph(RNA_test)
B = Plot.relation.graph(RNA_test)

A$densite = graph.density(A)
B$densite = graph.density(B)

V(A)$degre = degree(A)
E(A)$betw = edge.betweenness(A)



########### Appel de fonctions
data_expr = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
data_expr = data_expr[1:100,]
## Co-expression totale
Total_Micro = Make.adjacencyPVal(data_expr, Fast = F)
Plot.relation.graph(data_expr)

#Co-expression grp1 - grp 2
# Formation des groupes
N = ncol(data_expr)
data.Grp1 = data_expr[1:(N/2)]
data.Grp2 = data_expr[((N/2)+1):N]

#Graph d'adjacence
Grp1_Micro = Make.adjacencyPVal(data.Grp1)
Grp2_Micro = Make.adjacencyPVal(data.Grp2)

par(mfrow=c(1,3))
# Plot total
Plot.relation.graph(data_expr)
#Plot gr1
Plot.relation.graph(data.Grp1)
#Plot gr2
Plot.relation.graph(data.Grp2)

### Nanostring
raw.data = readRDS(file = "./DATA/NANOSTRING/Nanostring_Data.rds" )
group = table(raw.data$samples.IDs$tp53.status)
n1 = as.integer(group[1])
n2 = as.integer(group[2])
Nano = raw.data[["rcc.df"]]

#Co-expression grp1 - grp 2
# Formation des groupes
Nano_Mutated = Nano[1:n1]
Nano_WildType = Nano[(n1+1):(n1+n2)]

# Plot total
Total_Nano = Make.adjacency.graph(Nano)
Plot.relation.graph(Nano)
#Plot gr1
Nano_Mutated = Make.adjacency.graph(Nano_Mutated)
Plot.relation.graph(Nano_Mutated)
#Plot gr2
WT_Nano = Make.adjacency.graph(Nano_WildType)
Plot.relation.graph(Nano_WildType)



