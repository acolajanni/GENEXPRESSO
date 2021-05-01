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


############################################################################
# Fonctions
############################################################################

Cor.square.matrix <- function(data, method){
  
  
  if (!method%in%c("pearson","kendall","spearman")){
    stop("Enter something that switches me!")
  }
  
  print(c(method, "square matrix"))
  data_expr=as.data.frame(t(data))
  result_correlation=cor(data_expr, method = method)
  return(result_correlation)
}

TOM.square.matrix <- function(data){
  genes = row.names(data)
  data = t(as.data.frame(data))
  row.number = nrow(data)
  col.number = ncol(data)
  obs.number = col.number*row.number
  data[1:obs.number] = sapply(data[1:obs.number], as.numeric) 
  TOM = TOMsimilarityFromExpr(data)

  row.names(TOM) = genes
  colnames(TOM) = genes
  return(TOM)
}


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

Make.adjacencyPVal <-function(data, Fast = F, method){
  
  
  if (Fast == F){
    
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
    
    for (i in 1:nrow(interaction)){
      A = toString(interaction[["Var1"]][i])
      A = as.vector(data[A,])

      B = toString(interaction[["Var2"]][i])
      B = as.vector(data[B,])
    
      test = cor.test(as.numeric(A),as.numeric(B), method = method)
      cor = c(cor, test$p.value)

    }
    interaction = cbind(interaction, Pval = cor )
    colnames(interaction)[colnames(interaction) == "Pval"] = col.name 
  } 
  
  else if (Fast == T){
    data = as.matrix(t(data))
    Cor = corAndPvalue(data) 
    Pval = melt(Cor$p)
    Cor = melt(Cor$cor)
    interaction = cbind(Cor,PValue = Pval$value)
    colnames(interaction) = c("Var1","Var2","cor.pearson","PVal.pearson")
    
    filtre = as.character(interaction[,1])<as.character(interaction[,2])
    filtre2 = as.character(interaction[,1])<as.character(interaction[,2])
    interaction = interaction[filtre,]
  }
  
  return(interaction)
}

Make.full.adjacency <- function(data, PValue = T){
  
  if (PValue == T){
    interaction_spearman = Make.adjacencyPVal(data, Fast = F, method = "spearman")
    interaction_kendall = Make.adjacencyPVal(data, Fast = F, method = "kendall")
    interaction_Fast = Make.adjacencyPVal(data, Fast = T)
    
    interaction_Cor = merge(interaction_spearman,interaction_kendall,by = c("Var1","Var2"))
    interaction_Cor = merge(interaction_Cor,interaction_Fast,by = c("Var1","Var2"))                        
  }
  else {
    interaction_spearman =  Make.adjacency.table(data, method = "spearman")
    interaction_kendall =  Make.adjacency.table(data, method = "kendall")
    
    interaction_Cor = merge(interaction_spearman,interaction_kendall,by = c("Var1","Var2"))
  }
  
  interaction_TOM = Make.adjacency.table(data, method = "TOM")

  
  interaction = merge(interaction_Cor, interaction_TOM, by = c("Var1","Var2"))

  return(interaction)
}


Make.df.graph<-function(data, cor.threshold, Pvalue.threshold = F, method ){
  
  if (Pvalue.threshold == F){
    tools.graph <- switch(method,
                          spearman = {
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
                            Pvalue.threshold = F
                          },
                          
                          stop("Enter something that switches me!")
    )
  }  
   
  if (method == "TOM"){
    df.graph = subset(relations, relations[[cor]] >= cor.threshold)
  }
  else {
    relations[[cor]] = abs(relations[[cor]])
  }
  
  if (Pvalue.threshold == T){
    df.graph = subset(relations, (relations[[cor]] >= cor.threshold) & (relations[[pvalue]] <= 0.05) )
  }
  else {
    df.graph = subset(relations, relations[[cor]] >= cor.threshold)
  }
  
  df.graph = df.graph[-c(3:ncol(df.graph))]
  
  if (nrow(df.graph) != 0){
    df.graph = graph.data.frame(df.graph, directed = FALSE)
  }
  else{
    message = paste("No variable have been found having a", cor, "this high")
    stop(message)
  }
  
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

relations.comparison <- function(g1,g2,g1.name,g2.name, diplay){
  
  if (missing(g1.name)){
    g1.name = "g1"
  }
  if (missing(g2.name)){
    g1.name = "g2"
  }
  if(missing(diplay)){
    diplay = TRUE
  }
  
  diffg1 <- graph.difference(g1, g2)
  diffg2 <- graph.difference(g2, g1)
  interg1 <- graph.intersection(g1,g2, keep.all.vertices = T)
  
  G1_edges = as_data_frame(diffg1,what = "edges")
  
  if (nrow(G1_edges) != 0){
    G1_edges$color = "blue"
  }
  
  G2_edges = as_data_frame(diffg2,what = "edges")
  if (nrow(G2_edges) != 0){
    G2_edges$color = "darkgreen"
  }
  
  
  Common_edges = as_data_frame(interg1, what="edges")
  if (nrow(Common_edges) != 0){
    
    Common_edges = Common_edges[-c(3,4)]
    Common_edges$color = "red"
  }
  
  graph = rbind(G1_edges,G2_edges, Common_edges)

  g = graph.data.frame(graph, directed = F)
  
  if(diplay){
    plot(g, #layout = lay,
         edge.width = 2,
         vertex.size = 2,
         vertex.label = NA,
         vertex.color = "white",
         label.color = "black",
         label.font = 2)
    
    inter = paste(g1.name, g2.name, "intersection")
    
    legend("bottomleft",
           #x=-1.5, y=-1.1,
           c(g1.name,g2.name,inter), 
           pch=18, 
           col=c("blue","darkgreen","red"), 
           pt.cex=0, #taille des bulles légendes 
           cex=.8, #taille de la police légende
           lty=c(1,1,1),
           lwd = 3,
           bty="n", #absence de cadre autour de la légende 
           ncol=1)
  }
  
  return(g)
}

########################################################################################

########### Appel de fonctions
RNA = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ10000x30.csv", header = TRUE,row.names = 1)
RNA = RNA[1:1000,]

## Co-expression totale
Total_RNAseq = Make.full.adjacency(RNA, PValue = F)

spearman = Make.df.graph(RNA, cor.threshold = 0.9,Pvalue.threshold = F ,method = "spearman")
#plot(spearman)

TOM = Make.df.graph(RNA, cor.threshold = 0.33,method = "TOM")
#plot(TOM)

kendall = Make.df.graph(RNA,  cor.threshold = 0.75,Pvalue.threshold = F,method = "kendall")
#plot(kendall)

par(mfrow = c(1,3))
CompGraph_TOM_kendall = relations.comparison(TOM, kendall, "TOM", "Kendall")
CompGraph_TOM_spearman = relations.comparison(TOM, spearman, "TOM", "Spearman")
CompGraph_spearman_kendall = relations.comparison(spearman, kendall, "Spearman", "Kendall")


#Co-expression grp1 - grp 2
# Formation des groupes
N = ncol(RNA)
data.Grp1 = RNA[1:(N/2)]
data.Grp2 = RNA[((N/2)+1):N]

#Graph d'adjacence

#Méthode : Spearman
dev.off()
spearmanG1 = Make.df.graph(data.Grp1, cor.threshold = 0.85,Pvalue.threshold = F ,method = "spearman")
spearmanG2 = Make.df.graph(data.Grp2, cor.threshold = 0.85,Pvalue.threshold = F ,method = "spearman")

CompGraph_total = relations.comparison(spearmanG1, spearmanG2, "Control", "Case")

#Méthode : Kendall
kendallG1 = Make.df.graph(data.Grp1, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")
kendallG2 = Make.df.graph(data.Grp2, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")

CompGraph_total = relations.comparison(kendallG1, kendallG2, "Control", "Case")

#Méthode : TOM
TomG1 = Make.df.graph(data.Grp1, cor.threshold = 0.05,Pvalue.threshold = F ,method = "TOM")
TomG2 = Make.df.graph(data.Grp2, cor.threshold = 0.05,Pvalue.threshold = F ,method = "TOM")

CompGraph_total = relations.comparison(TomG1, TomG2, "Control", "Case")

# Tableaux d'adjacences
Grp1_RNAseq = Make.full.adjacency(data.Grp1, PValue = F)
Grp2_RNAseq = Make.full.adjacency(data.Grp2, PValue = F)


########################################################################################

########### Appel de fonctions
micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
#data_expr = data_expr[1:100,]

## Co-expression totale

Total_Micro = Make.full.adjacency(micro, PValue = F)

spearman = Make.df.graph(micro, cor.threshold = 0.8,Pvalue.threshold = F ,method = "spearman")
plot(spearman)

TOM = Make.df.graph(micro, cor.threshold = 0.2,method = "TOM")
plot(TOM)

kendall = Make.df.graph(micro,  cor.threshold = 0.65,Pvalue.threshold = F,method = "kendall")
plot(kendall)

par(mfrow = c(1,3))
CompGraph_TOM_kendall = relations.comparison(TOM, kendall, "TOM", "Kendall")
CompGraph_TOM_spearman = relations.comparison(TOM, spearman, "TOM", "Spearman")
CompGraph_spearman_kendall = relations.comparison(spearman, kendall, "Spearman", "Kendall")

#Co-expression grp1 - grp 2
# Formation des groupes
N = ncol(micro)
data.Grp1 = micro[1:(N/2)]
data.Grp2 = micro[((N/2)+1):N]

#Graph d'adjacence

#Méthode : Spearman
spearmanG1 = Make.df.graph(data.Grp1, cor.threshold = 0.85,Pvalue.threshold = F ,method = "spearman")
spearmanG2 = Make.df.graph(data.Grp2, cor.threshold = 0.85,Pvalue.threshold = F ,method = "spearman")

CompGraph_total = relations.comparison(spearmanG1, spearmanG2, "Control", "Case")

#Méthode : Kendall
kendallG1 = Make.df.graph(data.Grp1, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")
kendallG2 = Make.df.graph(data.Grp2, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")

CompGraph_total = relations.comparison(kendallG1, kendallG2, "Control", "Case")

#Méthode : TOM
TomG1 = Make.df.graph(data.Grp1, cor.threshold = 0.15,Pvalue.threshold = F ,method = "TOM")
TomG2 = Make.df.graph(data.Grp2, cor.threshold = 0.15,Pvalue.threshold = F ,method = "TOM")

CompGraph_total = relations.comparison(TomG1, TomG2, "Control", "Case")

# Tableaux d'adjacences
Grp1_Micro = Make.full.adjacency(data.Grp1, PValue = F)
Grp2_Micro = Make.full.adjacency(data.Grp2, PValue = F)


########################################################################################

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


## Co-expression totale
Total_nano = Make.full.adjacency(Nano, PValue = F)

spearman = Make.df.graph(Nano, cor.threshold = 0.85,Pvalue.threshold = F ,method = "spearman")
plot(spearman)

TOM = Make.df.graph(Nano, cor.threshold = 0.55,method = "TOM")
plot(TOM)

kendall = Make.df.graph(Nano,  cor.threshold = 0.75,Pvalue.threshold = F,method = "kendall")
plot(kendall)

par(mfrow = c(1,3))
CompGraph_TOM_kendall = relations.comparison(TOM, kendall, "TOM", "Kendall")
CompGraph_TOM_spearman = relations.comparison(TOM, spearman, "TOM", "Spearman")
CompGraph_spearman_kendall = relations.comparison(spearman, kendall, "Spearman", "Kendall")

#Co-expression grp1 - grp 2


#Méthode : Spearman
dev.off()
spearmanG1 = Make.df.graph(Nano_Mutated, cor.threshold = 0.9,Pvalue.threshold = F ,method = "spearman")
spearmanG2 = Make.df.graph(Nano_WildType, cor.threshold = 0.9,Pvalue.threshold = F ,method = "spearman")

CompGraph_total = relations.comparison(spearmanG1, spearmanG2, "Control", "Case")

#Méthode : Kendall
kendallG1 = Make.df.graph(Nano_Mutated, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")
kendallG2 = Make.df.graph(Nano_WildType, cor.threshold = 0.75,Pvalue.threshold = F ,method = "kendall")

CompGraph_total = relations.comparison(kendallG1, kendallG2, "Control", "Case")

#Méthode : TOM
TomG1 = Make.df.graph(Nano_Mutated, cor.threshold = 0.55,Pvalue.threshold = F ,method = "TOM")
TomG2 = Make.df.graph(Nano_WildType, cor.threshold = 0.55,Pvalue.threshold = F ,method = "TOM")

CompGraph_total = relations.comparison(TomG1, TomG2, "Control", "Case")

# Tableaux d'adjacences
Grp1_Micro = Make.full.adjacency(Nano_Mutated, PValue = F)
Grp2_Micro = Make.full.adjacency(Nano_WildType, PValue = F)









######################### recherches
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


