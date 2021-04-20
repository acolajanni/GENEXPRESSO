library(WGCNA)
library(cluster)

data = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)
data = as.data.frame(t(data))

Names = colnames(data)

design = as.numeric(rep(c(1,2),each = 6)) #y
plotClusterTreeSamples(datExpr=data, y=design)

cor.matrix= as.numeric(cor(design, data, use="p")) #GS1
P=corPvalueFisher(cor.matrix, nSamples =length(design) )

P2=P
P2[is.na(P)]=1
Q.std=qvalue(P2)$qvalues     

# Pas possible de déterminer le bruit on utilise des données traitées
GeneScreening=data.frame(Names,PearsonCorrelation=cor.matrix, P, Q.std)

# here we define the adjacency matrix using soft thresholding with beta=6
adjacency=abs(cor(data,use="p"))^6
K=softConnectivity(datE=data,power=6) 

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(K)
scaleFreePlot(K, main="Check scale free topology\n")



data=data[, rank(-K,ties.method="first" )<=100]

# Turn adjacency into a measure of dissimilarity
diss.Adj=1-adjacency
diss.tom=TOMdist(diss.Adj)


# Suite infaisable : pas de données d'activitation en fonctions des couleurs (microarray)

###

cmd=cmdscale(as.dist(diss.tom),2)
sizeGrWindow(7, 6)
par(mfrow=c(1,1))
plot(cmd,  main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

## Heatmap : TOM = topological overlap matrix similarity and dissimilarity
power=6
diss=1-TOMsimilarityFromExpr( data, power = 6 )
hier=hclust(as.dist(diss), method="average" )
diag(diss1) = NA;
sizeGrWindow(7,7)
TOMplot(diss^4, hier, 
        main = "TOM heatmap plot, module genes" )

# Heatmap 2
power=6
diss=1-adjacency( data, power = 6 )
hier=hclust(as.dist(diss), method="average" )
diag(diss) = NA;
sizeGrWindow(7,7)
TOMplot(diss^4, hier, 
        main = "Adjacency heatmap plot, module genes" )

#heatmap 3
topList=rank(NS1$p.Weighted,ties.method="first")<=30
top = rank(GeneScreening$PearsonCorrelation,ties.method = "first")<=30
topGene= names(data)[top]

#methode 1
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")
#Méthode 2
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="signed", useTOM=FALSE,
                   power=1, main="signed correlations")

# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="signed", useTOM=TRUE,
                   power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(data, plotGenes = topGene,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="D. TOM in an unsigned network")
