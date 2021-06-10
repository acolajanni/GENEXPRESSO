library(pheatmap)
library(hgu133plus2.db)


# https://slowkow.com/notes/pheatmap-tutorial/ 
################################################################################
################################################################################
# Exemple : 
example_file <- system.file ("extra/TagSeqExample.tab", package="DESeq")
data <- read.delim(example_file, header=T, row.names="gene")
data_subset <- as.matrix(data[rowSums(data)>50000,])

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

data_subset_norm <- t(apply(data_subset, 1, cal_z_score))
my_hclust_gene <- hclust(dist(data_subset), method = "complete")
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
set.seed(1984)
my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
my_gene_col$random <- my_random
my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(4,2)))
row.names(my_sample_col) <- colnames(data_subset)




my_heatmap <- pheatmap(data_subset,
                       annotation_row = my_gene_col,
                       annotation_col = my_sample_col,
                       cutree_cols = 2)


data_subset

################################################################################
################################################################################
library(RColorBrewer)
library(viridis)

load("./bgmasinter.RData")
load("./tab.RData")

test = bg.mas.inter
test2 = log2(test)

standardise <- function(x){
  (x - median(x))
}




my_sample_col = data.frame(sample = subset(tab$PreOpClinStage, paste0(tab$GEO,'.CEL.gz') %in% colnames(test)) )
row.names(my_sample_col) <- colnames(test)

test2 = t(apply(test2, 1, standardise))

my_heatmap <- pheatmap(test2,
                       annotation_col = my_sample_col,
                       cutree_cols = 2,
                       clustering_distance_cols	= "correlation",
                       clustering_method = "average",
                       show_rownames = FALSE
                       #, breaks = -20:20
                       , magma(length(mat_breaks) - 1)
                       , drop_levels = TRUE
                       )


mat_breaks <- seq(min(test2), max(test2), length.out = 10)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(test2, n = 10)

mat_breaks


################################################################################
################################################################################

load(file = "./data/inter&union.RData")
load(file = "./data/NormCompTotal.RData")


########### Comparaison background :


upset.bg = upset[, str_detect(colnames(upset) ,"background")]
upsetInter = Get.DEG.2(upset.bg, alternative=TRUE, method = "intersect")
upsetInter = c(upsetInter$Upregulated, upsetInter$Downregulated)

# récup / construction des jeux de données :
bg.mas = as.data.frame(dataset.bg$mas)
bg.rma = as.data.frame(dataset.bg$rma)

mapped.bg.mas = mapping.affymetrix.probe(bg.mas )
mapped.bg.rma = mapping.affymetrix.probe(bg.rma ) 

bg.rma.inter = mapped.bg.rma[ row.names(mapped.bg.rma) %in% UpsetGenes.intersect , ]
bg.rma.union = mapped.bg.rma[ row.names(mapped.bg.rma) %in% UpsetGenes.union , ]
  
bg.mas.inter = mapped.bg.mas[ row.names(mapped.bg.mas) %in% UpsetGenes.intersect , ]
bg.mas.union = mapped.bg.mas[ row.names(mapped.bg.mas) %in% UpsetGenes.union , ]



bg.mas.inter = mapped.bg.mas[ row.names(mapped.bg.mas) %in% upsetInter , ]
bg.mas.inter
###########################################

# Construction des heatmap
library(dendextend)
celpath = file.path("/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684")
txt.dir = paste0(celpath,"/GSE31684_table_of_clinical_details.txt")
tab = read.delim(txt.dir,check.names=FALSE,as.is=TRUE, header = T, fill = TRUE)
samples = subset(tab, tab$PreOpClinStage == 'T1' | tab$PreOpClinStage == 'T2') 
## grp T1 :
T1 = subset(samples, samples$PreOpClinStage == 'T1')
T1 = T1$GEO
T1 = paste0(T1,".CEL.gz")
## grp T2
T2 = subset(samples, samples$PreOpClinStage == 'T2')
T2 = T2$GEO
T2 = paste0(T2,".CEL.gz")




indiv <-  hclust(dist(t(bg.mas.inter)), method = "complete")
plot(indiv)
rect.hclust(indiv, k=2, border="red") 
cluster = cutree(indiv, k=2)
cluster
by(cluster, cluster, length)

indiv.clust <- data.frame(cluster = ifelse(test = cluster == 1, yes = "T1", no = "T2"))
#head(indiv.clust)

Original = c(T1,T2)
groups = c(rep("T1",19),rep("T2",55))
Orignal = data.frame(row.names = Original, CancerState = groups)


#total = merge(indiv.clust, Orignal, by = "row.names")
#row.names(total) = total$Row.names
#total$Row.names = NULL

my_data %>% arrange(Sepal.Length)


pheatmap(bg.mas.inter, annotation_col = Original, annotation_colors = c("red","blue") )


length(Orignal$CancerState)
ncol(bg.mas.inter)


####
Separe = by(indiv.clust, indiv.clust$cluster, row.names)
T1%in%Separe$T1
T2%in%Separe$T2
####

bg.mas.norm <- t(apply(bg.mas.inter, 1, cal_z_score))
pheatmap(bg.mas.norm,cutree_cols = 2)




