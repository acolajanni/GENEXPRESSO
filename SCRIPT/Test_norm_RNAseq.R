# Nom               : Test_norm_RNAseq
# Type              : Programme
# Objet             : Tester la normalisation (les fonctions iront dans functions.R)
# Input             : dataset RNA-seq
# Output            : ?
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 15.03.2021
#______________________________________________________________________________

library(edgeR)
library(DEFormats)
library(DESeq)
library("DESeq2")



#help("read.csv")

# importation de notre jeu de données
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
#Besoin de stocker le jeu de données sous forme de matrice
data.matrix(data)

#help("readDGE")
#help('DGEList')
#help('simulateRnaSeqData')

# On créer une DGElist pour pouvoir faire d'autres analyses
data = DGEList(count = data, 
        group = rep(1:2,each=ncol(data)/2)) #2 groupes (control / test) ==> 2 rep(1:2) / chaque grp a autant d'indiv chacun (ncol(data)/2)

#help("calcNormFactors.DGEList")

# On peut maintenant calculer des facteurs de normalisation : (on les calcule juste, on les applique pas donc voir comment faire ?)



# Méthode par défaut : TMM  + comparaison avec / sans cpm
TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")
TMM_norm
A = cpm(TMM_norm)
A
A = DGEList(count = A, 
               group = rep(1:2,each=ncol(data)/2))

Disp = estimateCommonDisp(A)
Disp = estimateTagwiseDisp(Disp)
A = exactTest(Disp)
topTags(A)


B <- calcNormFactors.DGEList(data, method = "TMM")
B = estimateCommonDisp(B)
B = estimateTagwiseDisp(B)
B = exactTest(B)
topTags(B)

# Methode utiles si le jeu de données est rempli de 0
TMMwsp_norm <- calcNormFactors.DGEList(data, method = "TMMwsp") 
TMMwsp_norm

# Relative log expression ==> moyenne geometric
RLE_norm <- calcNormFactors.DGEList(data, method = "RLE")
RLE_norm

# méthode normalisation par rapport au 75% quantile
upperquartile_norm <- calcNormFactors.DGEList(data, method = "upperquartile")
upperquartile_norm

#Maintenant, besoin d'estimer la dispersion puis tagwise dispersion
# Dispersion : permet d'établir un critère de dispersion des données de reads mappé sur un gène 
# Sans ce paramètre on peut pas établir si la différence de read mappé est due à une différence de longueur de gène ou a une expression différentielle
# Ensuite : Tagwise : établit la dispersion des données pour chacune des valeurs du dataset
Disp = estimateCommonDisp(TMM_norm)
Disp = estimateTagwiseDisp(Disp)

#help("estimateCommonDisp") 
#help('estimateTagwiseDisp')
#help(exactTest)
#help(topTags)

# Test des DEG avec une méthode proche du Test exact de fisher
DEG = exactTest(Disp)
DEG

####
#Récupérer les données normalisées en log2 : 
# A = cpm(DEGlist, log = TRUE)

# De cette manière on peut renvoyer un jeu de données normalisé analysable avec d'autres outils que edgeR

# on affiche par classement, les gènes les plus différentiellement exprimés
topTags(DEG, sort.by = 'PValue')


tools_norm_RNAseq.inspect <- function(raw.data,tool){
  data = as.matrix(raw.data)
  DEG = data.frame("genes" = row.names(data))
  
  tools_norm_RNAseq.fnc <- switch(tool,
                                  
        edgeR_TMM = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")
          tool_name = "TMM"
          #tmp = data.frame("TMM_norm" = TMM_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          #res.norm = cpm(TMM_norm)
        },
        
        edgeR_RLE = {
           DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
           res.norm <- calcNormFactors.DGEList(DGE, method = "RLE")
           tool_name = "RLE"
           #tmp = data.frame("RLE_norm" = RLE_norm$samples$norm.factors)
           #Norm_factors = cbind(Norm_factors, tmp)
           #res.norm = cpm(RLE_norm)
        },
        
        edgeR_TMMwsp = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          res.norm <- calcNormFactors.DGEList(DGE, method = "TMMwsp")
          tool_name = "TMMwsp"
          #tmp = data.frame("TMMwsp" = RLE_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          #res.norm = cpm(TMM_wsp)
        },
        
        edgeR_upperquartile = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          res.norm <- calcNormFactors.DGEList(DGE, method = "upperquartile")
          tool_name = "Upperquartile"
          #tmp = data.frame("upperquartile" = RLE_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          #res.norm = cpm(Upper)
        },
        
        deseq2.Wald = {
          conditions <-factor(c("condition1","condition1","condition1","condition1","condition1","condition1","condition2", "condition2","condition2","condition2","condition2","condition2"))
          type = factor(rep("single-read",12))
          colData = data.frame(conditions,type,row.names=colnames(data))
          dds<-DESeqDataSetFromMatrix(data,colData,design=~conditions)
          #Normalisation (median of ratio) + DEG
          #dds = estimateSizeFactors(dds)
          #dds = estimateDispersions(dds)
          DEG <- DESeq(dds, test = "Wald")
          DEG <- results(DEG)
          #DEG <- data.frame(deseq2_Wald=-log10(DEG$padj),SYMBOL=row.names(DEG))
          # à voir : Quel format pour les pvalues : log2, -log10 ...?
          DEG <- data.frame(deseq2_Wald=(DEG$padj),genes=row.names(DEG))
            
        },
        deseq2.LRT = {
          conditions <-factor(c("condition1","condition1","condition1","condition1","condition1","condition1","condition2", "condition2","condition2","condition2","condition2","condition2"))
          type = factor(rep("single-read",12))
          colData = data.frame(conditions,type,row.names=colnames(data))
          dds<-DESeqDataSetFromMatrix(data,colData,design=~conditions)
            #Normalisation (median of ratio) + DEG
          #dds = estimateSizeFactors(dds)
          #dds = estimateDispersions(dds)
          DEG <- DESeq(dds, test = "LRT",reduced = ~1)
          DEG <- results(DEG)
          #DEG <- data.frame(deseq2_Wald=-log10(DEG$padj),SYMBOL=row.names(DEG))
          # à voir : Quel format pour les pvalues : log2, -log10 ...?
          DEG <- data.frame(deseq2_LRT=(DEG$padj),genes=row.names(DEG))
          
        },
        deseq = {
          condition <-factor(c("condition 1","condition 1","condition 1","condition 1","condition 1","condition 1","condition 2", "condition 2","condition 2","condition 2","condition 2","condition 2"))
          colnames(data)
          type = factor(rep("single-read",12))
          design = data.frame(condition,type,row.names=colnames(data))
          singleSamples = design$type == 'single-read'
          countTable = data[,singleSamples]
          conds = design$condition[singleSamples]
          
          cds <- newCountDataSet(data, conds)
          cds = estimateSizeFactors(cds)
          cds = estimateDispersions(cds)
          DEG = nbinomTest(cds, condA = "condition 1", condB = "condition 2")
          DEG <- data.frame(deseq = (DEG$padj),genes= (DEG$id))
            
        },
        stop("Enter something that switches me!") 
        
  )
  if (!tool%in%c("deseq2.Wald","deseq2.LRT", "deseq")){
   # méthode DEG : Exact Test (edgeR)
    colname1 = paste(tool_name,"ExactTest")
    Disp = estimateCommonDisp(res.norm)
    Disp = estimateTagwiseDisp(Disp)
    pvalue = exactTest(Disp)$table$PValue
    DEG = data.frame(DEG,pvalue)
    
    #Méthode DEG : GLM (edgeR)
    colname2 = paste(tool_name,"GLM")
    design = model.matrix(~0+group, data = res.norm$samples)
    colnames(design) <- c("Control","Test")
    y = estimateDisp(res.norm,design)
    fit <- glmQLFit(y, design)
    BvsA <- makeContrasts(Control-Test, levels=design)
    qlf <- glmQLFTest(fit, contrast=BvsA)
    pvalue = qlf[["table"]][["PValue"]]
    DEG = data.frame(DEG,pvalue)
    
    names(DEG)[-1][-2] = colname1
    names(DEG)[-1][-1] = colname2
    DEG
  }
  return(DEG)
  
}

# importation de notre jeu de données
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
data.matrix(data)
data_to_comp = tools_norm_RNAseq.inspect(data,tool = 'edgeR_TMM')
head(data_to_comp)
################################################################ 
# De cette manière, on calcule les pvalue pour chacune de ces méthodes de normalisation

tools = c("edgeR_RLE","edgeR_upperquartile","edgeR_TMMwsp","deseq2.Wald","deseq2.LRT", "deseq")
for (tool in tools){
  print(tool)
  tmp = tools_norm_RNAseq.inspect(data,tool)
  data_to_comp = merge(data_to_comp,tmp,by = "genes",all=T)  
}

row.names(data_to_comp) <- data_to_comp$genes
data_to_comp <- data_to_comp[,-1]
data_to_comp = as.data.frame(t(data_to_comp))
head(data_to_comp)
head(data_to_comp)

################################################################

# Comparaison des jeux de données : fonction qui renvoit cpm() vs ne renvoit pas cpm()
data_to_comp = tools_norm_RNAseq.inspect(data,tool = 'edgeR_TMM')
data_to_comp

data_to_comp = DGEList(count = data_to_comp, 
               group = rep(1:2,each=ncol(data)/2))

Disp = estimateCommonDisp(data_to_comp)
Disp = estimateTagwiseDisp(Disp)
DEG = exactTest(Disp)
topTags(DEG, sort.by = 'PValue')

# Problème : avec cpm on estrait les données en log2 (comme annoncé mais les DEG ne sont plus significatifs)
# Pour comparer avec la même normalisation mais sans cpm(): 
# (pour les calculs j'ai enlevé le log2)

TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")
Disp_comp = estimateCommonDisp(TMM_norm)
Disp_comp = estimateTagwiseDisp(Disp_comp)
DEG_comp = exactTest(Disp_comp)
topTags(DEG_comp, sort.by = 'PValue')

#Pour comparer sans normalisation:
No_norm = estimateCommonDisp(data)
No_norm = estimateTagwiseDisp(No_norm)
DEG_no_norm = exactTest(No_norm)
topTags(DEG_no_norm, sort.by = 'PValue')

# Conclusion : n'utiliser cpm() que pour comparer les autres méthodes de DEG
# Les jeux de données cpm sont pareils entre les deux premiers
# Hypothèse : réutiliser la fonction DGElist() sur les données de cpm() applique des changements


# page 22/122 sur la doc de EdgeR
# GLM 
# Création de la matrice de design
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
data = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))


TMM_test <- calcNormFactors.DGEList(data, method = "TMM")


design = model.matrix(~0+group, data = data$samples)
design = model.matrix(~0+group, data = TMM_test$samples)
colnames(design) <- c("A","B")

design

#compare treatments
y_tmm = estimateDisp(TMM_test,design)
y = estimateDisp(data,design)

fit <- glmQLFit(y, design)
fit_tmm <- glmQLFit(y_tmm, design)

# 2 conditions (Control = A / test = B) : Test de DEG = makecontrast
# Utilise un modèle de log linéaire négative binomial pour observer des différences d'expr de gènes
BvsA <- makeContrasts(B-A, levels=design)
AvsB <- makeContrasts(A-B, levels=design)
qlf <- glmQLFTest(fit, contrast=BvsA)
qlf2 <- glmQLFTest(fit, contrast=AvsB)
# A vs B = B vs A
qlf
topTags(qlf)
topTags(qlf2)

### Au final : on peut appliquer le GLM sur n'importe quel facteur de normalisation si c'est une DEGlist
qlf_tmm <- glmQLFTest(fit_tmm, contrast=BvsA)
qlf <- glmQLFTest(fit, contrast=BvsA)
topTags(qlf_tmm)
topTags(qlf)
qlf_tmm$table

################# ################# #################  DESeq(1)
BiocManager::install("DESeq")
library(DESeq)
### Pareil que DESeq2
SimulRNASEQ = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
condition <-factor(c("condition 1","condition 1","condition 1","condition 1","condition 1","condition 1","condition 2", "condition 2","condition 2","condition 2","condition 2","condition 2"))
type = factor(rep("single-read",12))
design = data.frame(condition,type,row.names=colnames(SimulRNASEQ))

singleSamples = design$type == 'single-read'
countTable = SimulRNASEQ[,singleSamples]
conds = design$condition[singleSamples]
conds
### Absolument équivalent à DESeq2
cds <- newCountDataSet(SimulRNASEQ, conds)
cds = estimateSizeFactors(cds)

# Normalized = True : Divise chaque valeur par le sizeFactor associé
counts(cds, normalized = TRUE)
# établir la dispersion
cds = estimateDispersions(cds)

cds
res = nbinomTest(cds, condA = "condition 1", condB = "condition 2")
res
help("nbinomTest")
# most downregulated :
head( res[ order( res$foldChange, -res$baseMean ), ] )

# most upregulated
head( res[ order( -res$foldChange, -res$baseMean ), ] )
# On remarque que très peu de p-value sont significative
# MAIS on a travaillé sur un jeu de donnée peu transformé
# essayons avec le GLM
# GLM : Ne fonctionne qu'en 2 facteurs (= ici qu'un facteur puisqu'on a réglé tout sur "single reads")
# Le facteur uniuqe étant la condition : test vs control

# Conclusion : DESeq sert pas à grand chose
