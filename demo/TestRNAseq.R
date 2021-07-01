library(edgeR)
library(DESeq)
library(DESeq2)
library(limma)
library(RankProd)

#########################
#########################
#########################
#########################
# EdgeR avec la normalisation de DESeq2 : 
data = Simul.data("RNAseq",n.cond1 = 15, n.cond2 = 15,nb.genes = 1000)

condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,type,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

# Calcul des size factors de DESeq2
SizeF = estimateSizeFactors(dds)
SizeF = sizeFactors(SizeF)
names(SizeF) = NULL

# Analyse de EdgeR
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2),
              norm.factors = SizeF) # functions for model with n1 = n2 


Disp = estimateCommonDisp(DGE)
Disp = estimateTagwiseDisp(Disp)

# Application de size factors de DESeq2
Disp$samples$norm.factors = SizeF
# Vérification
Disp$samples #ça fonctionne

pvalue = exactTest(Disp)
Norm.DESeq.ExactTest = topTags(pvalue)
Norm.DESeq.ExactTest = Norm.DESeq.ExactTest$table

#########################
#########################
#########################
#########################
# DESeq2 avec la normalisation de TMM de edgeR 

data = Simul.data("RNAseq",n.cond1 = 15, n.cond2 = 15,nb.genes = 1000)

# Dataset into DEGlist (edgeR class)
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 
# Normalization 
res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")
NormF = res.norm$samples$norm.factors

condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,type,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

dds = estimateSizeFactors(dds)
##
avant = sizeFactors(dds)
sizeFactors(dds) = NormF
##
apres = sizeFactors(dds)
# On vérifie les changements : 
avant
apres # ça fonctionne


DEG <- DESeq(dds, test = "Wald")
DEG <- results(DEG)

DEG= data.frame(DEG$log2FoldChange, DEG$padj, row.names = row.names(DEG))
Norm.TMM.Wald = DEG[order(DEG$DEG.padj),]
Norm.TMM.Wald = Norm.TMM.Wald[1:10,]



######## Méthodes standard : 
#### EdgeR TMM

# Dataset into DEGlist (edgeR class)
DGE = DGEList(count = data, 
              group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 
# Normalization 
res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")

Disp = estimateCommonDisp(res.norm)
Disp = estimateTagwiseDisp(Disp)
# Retrieving pvalues
pvalue = exactTest(Disp)
TMM = topTags(pvalue)$table



#### DESeq2 Wald

#data into DESeqDataSet (DESeq2 class) :
condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )

type = factor(rep("single-read",ncol(data)))
colData = data.frame(condition,row.names=colnames(data))
dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

# Retrieve DEG analysis (test Wald)
DEG <- DESeq(dds, test = "Wald")
DEG <- results(DEG)

DEG= data.frame(DEG$log2FoldChange, DEG$padj, row.names = row.names(DEG))
Wald = DEG[order(DEG$DEG.padj),]
Wald = Wald[1:10,]



#### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### 
#load("./data/RNAseqTEST.RData")

## Load data
load("./data/CountTable.RData")
TCGA.dir = "/autofs/unitytravail/travail/acolajanni/BLCA_TCGA/"
TCGA.table = paste0(TCGA.dir,"gdc_sample_sheet.2021-06-04.tsv")

tab = read.delim(TCGA.table,check.names=FALSE,as.is=TRUE, header = T, fill = TRUE)




keep = rowSums(table.count) >= (1*ncol(table.count))
table.count = table.count[keep,]




#count.matrix = table.count.subset

design = factor(tab$`Sample Type`)
class = ifelse(design == "Primary Tumor", yes = 1, no = 2)
#class = design

#Norm = tools.norm.RNAseq(table.count, tool = "voom", design = class)
#DEG = tools.DEG.RNAseq(table.count, nf = Norm,
#                     tool.norm = "vst", 
#                     tool.DEG = "limma", 
#                     design = class)



A = tools.norm.RNAseq(table.count, tool = "TMM" , design = class)
tmp = tools.DEG.RNAseq(table.count, nf = A,
                       tool.norm = "TMM", 
                       tool.DEG = "nbinom", 
                       design = class)



tools.norm = c("TMM", "TMMwsp", "RLE", "Upperquartile", "voom", "vst", "vst2")
tools.DEG = c("ExactTest","GLM",
              #"nbinom",
              "nbinom.Wald","nbinom.LRT")


data.to.comp = data.frame(NULL)

for (norm in tools.norm){
  
  Norm = tools.norm.RNAseq(table.count, tool = norm , design = class)
  
  print(norm)
  
  if (norm %in% c("vst2","voom")){
    tmp = tools.DEG.RNAseq(table.count, nf = Norm,
                           tool.norm = norm, 
                           tool.DEG = "limma", 
                           design = class)
    
    data.to.comp = merge(data.to.comp, tmp, by = "SYMBOL")
  }
  else{
    
    for (deg in tools.DEG){
      
      print(deg)
      tmp = tools.DEG.RNAseq(table.count, nf = Norm,
                             tool.norm = norm, 
                             tool.DEG = deg, 
                             design = class)
      
      if (dim(data.to.comp)[1] == 0){
        data.to.comp = tmp
      }
      else{
        data.to.comp = merge(data.to.comp, tmp, by = "SYMBOL")
      }
      
    }

  }
}  
  
  
row.names(data.to.comp) = data.to.comp$SYMBOL
data.to.comp$SYMBOL = NULL  
  
save(data.to.comp, file = "./data/RNAseqDataToComp.RData")
load(file = "./data/RNAseqDataToComp.RData")

data.to.comp = as.data.frame(t(data.to.comp))

PCA_tools(data.to.comp)
    
# Liste de toute les méthodes
methods = row.names(data.to.comp)
# Dataframe rempli de valeur binaire (0/1)
upset = Upset.Binary.Dataframe(data.to.comp)


upset(upset,
      sets = methods,
      sets.bar.color = "#56B4E9", 
      order.by = "freq", 
      text.scale = 1.2,
      mb.ratio = c(0.6,0.4),
      set_size.show = TRUE,
      set_size.scale_max = 25000
      )


d <- dist(data.to.comp, method = "euclidean") # distance matrix
fit <- hclust(d, method="complete") 

subset_cluster = cutree(fit, k=3)

colors = paletteer_d("yarrr::xmen")[1:3]
bars = as.data.frame(subset_cluster)
bars$subset_cluster[bars$subset_cluster == 1] = colors[1]
bars$subset_cluster[bars$subset_cluster == 2] = colors[2]
bars$subset_cluster[bars$subset_cluster == 3] = colors[3]
bars$subset_cluster[bars$subset_cluster == 4] = colors[4]
bars$subset_cluster[bars$subset_cluster == 5] = colors[5]
bars = as.matrix(bars)


# allow content to go into outer margin 
par(mar = c(15, 4, 4, 3) + 0.1,
    xpd = TRUE) 
plot(fit)
fit = as.dendrogram(fit)
rect.dendrogram(fit, k = 3, border = c(colors[3],colors[1],colors[2]))

colored_bars(colors = bars, dend = fit, rowLabels = "Cluster")

legend("topright", 
       legend = c('edgeR', "DESeq2" ,"limma"), 
       pch = 15, 
       pt.cex = 3, 
       cex = 1,
       bty = "o",
       inset = c(0, -0.49), 
       title = "Cluster", 
       col = colors,
       xpd = TRUE
       #horiz = TRUE
)









A = tools.DEG.RNAseq.merge(Data,tools.norm = c("TMM","TMMwsp"),  tools.DEG = "all", design = design)



#' Merge the DEG Pvalues for each tool used by tools.DEG.Microarrays in one dataframe
#'
#' @param data Dataframe with genes in row, and methods used in columns. It contains the differentially expressed p-values for each gene.
#' @param tools Different methods are used to compute pvalues of differentially expressed genes for microarrays
#' "Wilcox" uses the wilcoxDEG() function implemented in this very same package 
#' "limma" and "GEOlimma uses respectively the functions DEG_limma() and DEG_GEOlimma() that come from the limma package
#' "RankProduct","RankProduct.log" perform a Rank Product analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 
#' "RankSum","RankSum.log" perform a Rank Sum analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 

#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return Dataframe of pvalues of genes being differentially expressed with genes in columns and methods in rows
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
#' tools = c("limma", "Wilcox","RankProduct","RankProduct.log","RankSum","RankSum.log")
#' res.DEG = tools.DEG.Microarrays.merge(Data,tools,10,10)
tools.DEG.Microarrays.merge <-function(data,tools,n1,n2){
  # By default, all the implemented methods are used
  if (missing(tools)){
    tools = c("limma", "GEOlimma", "Wilcox","RankProduct","RankProduct.log","RankSum","RankSum.log")
  }
  if (missing(n1)){
    n1=ncol(data)/2
  }
  if (missing(n2)){
    n2=ncol(data)/2
  }
  # Creating the dataframe with only genes in it to merge further results with it
  data_to_comp = data.frame(Gene.ID = row.names(data))
  
  for (tool in tools){
    print(tool)
    tmp = tools.DEG.Microarrays(data,tool,n1, n2)
    data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
  }
  # adding the genes as row names
  row.names(data_to_comp) <- data_to_comp$Gene.ID
  data_to_comp <- data_to_comp[,-1]
  # Genes in columns, methods in rows
  data_to_comp = as.data.frame(t(data_to_comp))
  return(data_to_comp)
}
