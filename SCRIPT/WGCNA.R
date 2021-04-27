library(WGCNA)
library(cluster)

## Microarray
data = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)
design = as.numeric(rep(c(1,2),each = 6)) #y
data = as.data.frame(t(data))
## Nanostring
raw.data = readRDS(file = "./DATA/NANOSTRING/Nanostring_Data.rds" )
#data.dir <- "./DATA/NANOSTRING"
#RCC.dir <- file.path(data.dir,"GSE146204_RAW")
#raw <- RCC.dir
samples.IDs <- raw.data$samples.IDs
rcc.samples <- raw.data$rcc.df
annots.df <- raw.data$annots.df
samples.IDs <- raw.data$samples.IDs
# design matrix for nanostring
group = table(samples.IDs$tp53.status)
n1 = as.integer(group[1])
n2 = as.integer(group[2])
design = as.numeric(c(rep(1,n1),rep(2,n2)))
data <- as.data.frame(rcc.samples[grepl("Endog",annots.df$CodeClass),])

data = as.data.frame(t(data))
#data[1:784] = lapply(data[1:784], as.numeric) 
data[1:750] = lapply(data[1:750], as.numeric) 

##


Names = colnames(data)

plotClusterTreeSamples(datExpr=data, y=design)

cor.matrix= as.numeric(cor(design, data, use="p")) #GS1
cor.matrix2 = cor(design,data,use='p')

P=corPvalueFisher(cor.matrix, nSamples =length(design) )

P2=P
P2[is.na(P)]=1
Q.std=qvalue(P2)$qvalues     


# Pas possible de déterminer le bruit on utilise des données traitées
GeneScreening=data.frame(Names,PearsonCorrelation=cor.matrix, P, Q.std)
head(GeneScreening)

# here we define the adjacency matrix using soft thresholding with beta=6
adjacency=abs(cor(data,use="p"))^6
K=softConnectivity(datE=data,power=6) 

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(K)
scaleFreePlot(K, main="Check scale free topology\n")



data=data[, rank(-K,ties.method="first" )<=50]
# Turn adjacency into a measure of dissimilarity
diss.Adj=1-adjacency
diss.tom=TOMdist(diss.Adj)


###

cmd=cmdscale(as.dist(diss.tom),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd,  main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

## Heatmap : TOM = topological overlap matrix similarity and dissimilarity
power=6
diss=1-TOMsimilarityFromExpr( data)
hier=hclust(as.dist(diss), method="average" )
plot(hier)
diag(diss) = NA;
sizeGrWindow(7,7)
TOMplot(diss^4, hier, 
        main = "TOM heatmap plot, module genes" )
help(TOMplot)
# Heatmap 2
power=6
diss=1-adjacency( data, power = 6 )
hier=hclust(as.dist(diss), method="average" )
diag(diss) = NA;
sizeGrWindow(9,7)
TOMplot(diss^4, hier,
        main = "Adjacency heatmap plot, module genes" )

#heatmap 3
#topList=rank(NS1$p.Weighted,ties.method="first")<=150
top = rank(GeneScreening$PearsonCorrelation,ties.method = "first")<=750
topGene= names(data)[top]

#methode 1 : Nous intéresse c'est la corrélation, pas forcément le signe  
# On peut comparer la méthode utilisant des coeff de corr. de Pearson 
# Avec TOM Topological overlap matrix
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="unsigned", useTOM=FALSE,
                   power=6, main="unsigned correlations")



# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="TOM in an unsigned network")

# _________________________________________________________________________

Coexpression<-function(datatype, plot, gene.number){
  if(datatype == "Nanostring"){
    raw.data = readRDS(file = "./DATA/NANOSTRING/Nanostring_Data.rds" )
    samples.IDs <- raw.data$samples.IDs
    rcc.samples <- raw.data$rcc.df
    annots.df <- raw.data$annots.df
    samples.IDs <- raw.data$samples.IDs
    # design matrix for nanostring
    group = table(samples.IDs$tp53.status)
    n1 = as.integer(group[1])
    n2 = as.integer(group[2])
    design = as.numeric(c(rep(1,n1),rep(2,n2)))
    data <- as.data.frame(rcc.samples[grepl("Endog",annots.df$CodeClass),])
    
    data = as.data.frame(t(data))
    #data[1:784] = lapply(data[1:784], as.numeric) 
    size = dim(data)[2]

  } 
  else if (datatype == "Microarrays"){
    data = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
    design = as.numeric(rep(c(1,2),each = 12)) #y
    data = as.data.frame(t(data))
    size = dim(data)[2]
  }
  else if (datatype == "RNAseq"){
    data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ1000x30.csv", header = TRUE,row.names = 1)
    design = as.numeric(rep(c(1,2),each = 15)) #y
    data = as.data.frame(t(data))
    size = dim(data)[2]
    
  }
  else{
    stop("Enter something that switches me!")   
  }
  data[1:size] = lapply(data[1:size], as.numeric) 
  gene.names = colnames(data)
  cor.matrix= as.numeric(cor(design, data, use="p"))
  cor.matrix2 = cor(design,data,use='p')
  P=corPvalueFisher(cor.matrix, nSamples =length(design) )
  P2=P
  P2[is.na(P)]=1
  Q.std=qvalue(P2)$qvalues     
  GeneScreening=data.frame(gene.names,PearsonCorrelation=cor.matrix, P, Q.std)
  
  # 
  top = rank(GeneScreening$PearsonCorrelation,ties.method = "first")<=gene.number
  topGene= names(data)[top]
  
  if (plot == "TOM module genes"){
    K=softConnectivity(datE=data,power=6) 
    data=data[, rank(-K,ties.method="first" )<=gene.number]
    
    diss=1-TOMsimilarityFromExpr( data)
    hier=hclust(as.dist(diss), method="average" )
    plot(hier)
    diag(diss) = NA;
    sizeGrWindow(7,7)
    TOMplot(diss^4, 
            hier, 
            main = "TOM heatmap plot, module genes" )
  }
  else if (plot == "Adj module genes"){
    K=softConnectivity(datE=data,power=6) 
    data=data[, rank(-K,ties.method="first" )<=gene.number]
    
    diss=1-adjacency( data, power = 6 )
    hier=hclust(as.dist(diss), method="average" )
    diag(diss) = NA;
    sizeGrWindow(7,7)
    TOMplot(diss^4, hier,
            main = "Adjacency heatmap plot, module genes" )
  }
  else if (plot == "Adjacency network"){
    plotNetworkHeatmap(data, plotGenes = topGene,
                       networkType="unsigned", useTOM=FALSE,
                       power=6, main="unsigned correlations")
  }
  else if (plot == "TOM network"){
    plotNetworkHeatmap(data, plotGenes = topGene,
                       networkType="unsigned", useTOM=TRUE,
                       power=6, main="TOM in an unsigned network")
  }
  else{
    stop("Enter something that switches me!")   
  }
  
}
sizeGrWindow(10,5)
par(mfrow=c(1,2))
Coexpression("Nanostring","TOM module genes", 500)
Coexpression("Nanostring","Adj module genes", 500)
Coexpression("Nanostring","TOM network", 500)
Coexpression("Nanostring","Adjacency network", 500)
