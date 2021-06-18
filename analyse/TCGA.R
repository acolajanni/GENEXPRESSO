############################################
######## Load  (GenomicDataCommons) ########
############################################
library(GenomicDataCommons)

TCGA.dir = "/autofs/unitytravail/travail/acolajanni/BLCA_TCGA/"
TCGA.table = paste0(TCGA.dir,"gdc_sample_sheet.2021-06-04.tsv")
counts.dir = paste0(TCGA.dir, "gdc_download_20210604_095100.950304/")

tab = read.delim(TCGA.table,check.names=FALSE,as.is=TRUE, header = T, fill = TRUE)

# Récupération des ID et noms de fichiers
directories = tab$`File ID`
file.names = tab$`File Name`
# Initialisation du dataframe de counts
table.count = data.frame()
for (i in 1:length(directories)){
  file.dir = paste0(counts.dir,directories[i],"/")
  count = paste0(file.dir,file.names[i])
  # lecture du DF de count pour un échantillon
  tmp = readHTSeqFile(count, samplename = directories[i])
  
  if(dim(table.count)[1] == 0){
    table.count = tmp
  }
  # Fusion des colonnes de comptes pour chaque échantillon
  else {
    table.count = merge(table.count, tmp, by = "feature")
  }
} 
# Certaines lignes ne correspondent pas à des comptes, on les retire :
to.skip = c("__alignment_not_unique","__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")
table.count = subset(table.count, !table.count$feature %in% to.skip )

save(table.count, file = "./data/CountTable.RData")
load("./data/CountTable.RData")

############################################
############## Load  (DESeq2) ##############    
############################################
library(DESeq2)


TCGA.dir = "/autofs/unitytravail/travail/acolajanni/BLCA_TCGA/"
TCGA.table = paste0(TCGA.dir,"gdc_sample_sheet.2021-06-04.tsv")
counts.dir = paste0(TCGA.dir, "gdc_download_20210604_095100.950304/")

#list.files(counts.dir)
tab = read.delim(TCGA.table,check.names=FALSE,as.is=TRUE, header = T, fill = TRUE)
file.names = paste0(tab$`File ID`,"/",tab$`File Name`)
sampleTable <- data.frame(sampleName = tab$`File ID`,
                          fileName = file.names,
                          condition = tab$`Sample Type`)


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = counts.dir,
                                       design= ~ condition)

table.count = as.data.frame(counts(ddsHTSeq))


ens = row.names(table.count)


############################################
###########        ENSGR         ###########
############################################
library(stringr)
ENSGR = ens[str_detect(ens, "ENSGR")]
write.csv(ENSGR, file="./ENSGR", row.names = FALSE)

############################################
###########        Mapping       ###########
############################################

library(EnsDb.Hsapiens.v79)
ENSEMBL_EnsDb <- keys(EnsDb.Hsapiens.v79, keytype = "GENEID")
test = as.data.frame(ENSEMBL_EnsDb)

gene_dataframe_EnsDb <- ensembldb::select(EnsDb.Hsapiens.v79, keys=ENSEMBL_EnsDb, 
                                          columns=c("SYMBOL", "GENEBIOTYPE"), keytype="GENEID")

ProteinCoding = subset(gene_dataframe_EnsDb, gene_dataframe_EnsDb$GENEBIOTYPE == "protein_coding")
CodingProt = ProteinCoding[,c("GENEID","SYMBOL")]
# Nb d'identifiants uniques :
sapply(CodingProt, function(x) length(unique(x)))

############################################
# Test de filtre (d'après l'API de DESeq2) #
############################################
table.count2 = table.count
row.names(table.count2) = table.count$feature
#table.count2$feature = NULL

keep = rowSums(table.count2) >= (10*ncol(table.count2))
table.count2 = table.count2[keep,]
dim(table.count2) #32.568 ID (toujours trop)


############################################
#  Création du sous-jeu de données (test)  #
############################################

normal.ID = tab$`File ID`[tab$`Sample Type` == "Solid Tissue Normal"]
case.ID = tab$`File ID`[tab$`Sample Type` == "Primary Tumor"][1:19]

tab.subset = tab[tab$`File ID` %in% c(normal.ID,case.ID),]

table.count.subset = table.count[colnames(table.count) %in% tab.subset$`File ID`]
table.count.subset = table.count.subset[1:2000,]

save(tab.subset, table.count.subset, file="./data/RNAseqTEST.RData")
load(file = "./data/RNAseqTEST.RData")


library(GENEXPRESSO)
library(edgeR)
library(limma)
count.matrix = table.count.subset




nf = calcNormFactors(count.matrix, method = "TMM")

voom.data = voom(count.matrix, 
                 design = model.matrix(~factor(tab.subset$`Sample Type`)),
                 lib.size = colSums(count.matrix) * nf)


voom.data$genes = rownames(count.matrix)
voom.fitlimma = lmFit(voom.data, design = model.matrix(~factor(tab.subset$`Sample Type`)))
voom.fitbayes = eBayes(voom.fitlimma)
voom.pvalues = voom.fitbayes$p.value[, 2]
voom.adjpvalues = p.adjust(voom.pvalues, method = "BH")
voom.adjpvalues

