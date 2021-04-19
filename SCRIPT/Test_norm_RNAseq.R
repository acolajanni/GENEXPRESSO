# Nom               : Test_norm_RNAseq
# Type              : Programme
# Objet             : Normaliser + Analyser data RNAseq
# Input             : dataset RNA-seq
# Output            : data_to_comp 
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 15.03.2021
#______________________________________________________________________________

source(file.path("./SCRIPT","functions.R"))


# Fonction combinant la normalation + l'analyse DEG
tools_norm_RNAseq.inspect <- function(data,tool){
  # Besoin de mettre le jeu de données dans les bons formats
  data = as.matrix(data)
  # Créeation d'un tableau à une colonne : "Gènes" qui liste tous les gènes
  res.diff = data.frame("genes" = row.names(data))
  
  tools_norm_RNAseq.fnc <- switch(tool,
                                  
        edgeR_TMM = {
          # Jeu de données dans la classe de EdgeR (DGElist)
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          # Normalisation :
          res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")
          tool_name = "TMM"
        },
        
        edgeR_RLE = {
           DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
           res.norm <- calcNormFactors.DGEList(DGE, method = "RLE")
           tool_name = "RLE"
        },
        
        edgeR_TMMwsp = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          res.norm <- calcNormFactors.DGEList(DGE, method = "TMMwsp")
          tool_name = "TMMwsp"
        },
        
        edgeR_upperquartile = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          res.norm <- calcNormFactors.DGEList(DGE, method = "upperquartile")
          tool_name = "Upperquartile"
        },
        
        deseq2.Wald = {
          # Mise en forme des données pour créer un DESeqDataSet :
          condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
          type = factor(rep("single-read",ncol(data)))
          colData = data.frame(condition,type,row.names=colnames(data))
          dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)
          
          # Récupération des résultats avec l'analyse DEG (test Wald)
          DEG <- DESeq(dds, test = "Wald")
          DEG <- results(DEG)
          
          #DEG <- data.frame(deseq2_Wald=-log10(DEG$padj),SYMBOL=row.names(DEG))
          # à voir : Quel format pour les pvalues : log2, -log10 ...?
          # Rajout d'une colonne au tableau 
          res.diff <- data.frame(deseq2_Wald=(DEG$padj),genes=row.names(DEG))
            
        },
        
        deseq2.LRT = {
          condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
          type = factor(rep("single-read",ncol(data)))
          colData = data.frame(condition,type,row.names=colnames(data))
          dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)

          DEG <- DESeq(dds, test = "LRT",reduced = ~1)
          DEG <- results(DEG)
          
          res.diff <- data.frame(deseq2_LRT=(DEG$padj),genes=row.names(DEG))
        },
        
        deseq = {
          condition <- rep(c("control","case"),each = (ncol(data)/2) )
          type = factor(rep("single-read",ncol(data)))
          design = data.frame(condition,type,row.names=colnames(data))
          singleSamples = design$type == 'single-read'
          countTable = data[,singleSamples]
          conds = design$condition[singleSamples]
          
          cds <- newCountDataSet(data, conds)
          cds = estimateSizeFactors(cds)
          cds = estimateDispersions(cds)
          DEG = nbinomTest(cds, condA = "control", condB = "case")
          res.diff <- data.frame(deseq = (DEG$padj),genes= (DEG$id))
            
        },
        stop("Enter something that switches me!") 
        
  )
  # Tests DEG pour les données normalisées non traités (on exclut donc DESeq)
  if (!tool%in%c("deseq2.Wald","deseq2.LRT", "deseq")){
   # méthode DEG : Exact Test (edgeR)
    colname1 = paste(tool_name,"ExactTest")
    # Calculs des dispersions
    Disp = estimateCommonDisp(res.norm)
    Disp = estimateTagwiseDisp(Disp)
    # Excat test + récupération des PValue + gènes correspondants
    pvalue = exactTest(Disp)$table[3]
    
    # tableau 2 colonnes : Gènes - pvalue
    res.diff1 = data.frame(genes = row.names(pvalue),pvalue = pvalue$PValue)

    #Méthode DEG : GLM (edgeR)
    colname2 = paste(tool_name,"GLM")
    # Matrice de design : 
    design = model.matrix(~0+group, data = res.norm$samples)
    colnames(design) <- c("Control","Test")
    # caluls des dispersions
    y = estimateDisp(res.norm,design)
    fit <- glmQLFit(y, design)
    BvsA <- makeContrasts(Control-Test, levels=design)
    # Calcul des DEG :
    qlf <- glmQLFTest(fit, contrast=BvsA)
    # On récupère les PValue + gènes correspondants
    pvalue = qlf[["table"]][4]
    res.diff2 = data.frame(genes = row.names(pvalue),pvalue = pvalue$PValue)
    res.diff = merge(res.diff1,res.diff2,by = "genes",all=T)
    
    # On renomme les deux dernières colonnes créées
    names(res.diff)[-1][-2] = colname1
    names(res.diff)[-1][-1] = colname2
  }
  return(res.diff)
  
}
################################################################ 
# Appel de la fonction : 
# importation de notre jeu de données
#data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ1000x30.csv", header = TRUE,row.names = 1)


# Essai de la fonction sur un paramètre
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
#on récupère le tableau semblable à MakeComparisonTable.R
data_to_comp = as.data.frame(t(data_to_comp))

############################################################
# PCA : 

PCA_tools(data_to_comp)

############################################################
# UpsetPlot : 

UpsetPlot(data.to.comp = data_to_comp,threshold = 0.05)#, empty.intersections = "yes")

#____________________________________________________________
# Pipeline d'analyse en fonction des packages :

################# ################# #################  EdgeR :

# importation de notre jeu de données
data = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)

# On créer une DGElist pour pouvoir faire d'autres analyses
data = DGEList(count = data, group = rep(1:2,each=ncol(data)/2)) 
#2 groupes (control / test) ==> 2 rep(1:2) / chaque grp a autant d'indiv chacun (ncol(data)/2)

################################ Méthodes normalisation EdgeR :

# Méthode TMM
TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")
# Methode utiles si le jeu de données est rempli de 0
TMMwsp_norm <- calcNormFactors.DGEList(data, method = "TMMwsp") 
# Relative log expression ==> moyenne geometric
RLE_norm <- calcNormFactors.DGEList(data, method = "RLE")
# méthode normalisation par rapport au 75% quantile
upperquartile_norm <- calcNormFactors.DGEList(data, method = "upperquartile")

################################ Méthodes analyse DEG EdgeR :
## Exact Test :
#--------------

# Analyse DEG necessite l'estimation des paramètres de dispersion : 
Disp = estimateCommonDisp(TMM_norm)
Disp = estimateTagwiseDisp(Disp)

# Test des DEG avec une méthode proche du Test exact de fisher
DEG = exactTest(Disp)
# on affiche par classement, les gènes les plus différentiellement exprimés
topTags(DEG, sort.by = 'PValue')

##    GLM     :
#--------------

# a besoin d'un autre paramètre : matrice de "design"
design = model.matrix(~0+group, data = TMM_norm$samples)
colnames(design) <- c("Control","Test")

# Estimation de la dispersion à partir du design
y = estimateDisp(TMM_norm,design)
fit <- glmQLFit(y, design)
BvsA <- makeContrasts(Control-Test, levels=design)
# Calcul des DEG
qlf <- glmQLFTest(fit, contrast=BvsA)
# on affiche par classement, les gènes les plus différentiellement exprimés
topTags(qlf)

#________________________________________________________
# Test de la fonction cpm pour exportation des données normalisées
# Méthode par défaut : TMM  + comparaison avec / sans cpm

# Avec cpm() : 
TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")

A = cpm(TMM_norm)
A = DGEList(count = A, group = rep(1:2,each=ncol(data)/2))
Disp = estimateCommonDisp(A)
Disp = estimateTagwiseDisp(Disp)
A = exactTest(Disp)
topTags(A)

# Avec cpm()
B <- calcNormFactors.DGEList(data, method = "TMM")
B = estimateCommonDisp(B)
B = estimateTagwiseDisp(B)
B = exactTest(B)
topTags(B)
#____________________________________________________


################# ################# #################  DESeq(1)

### On récupère le jeu de données :
SimulRNASEQ = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)
# Besoin de créer un data frame de design :
condition <-factor(c("condition 1","condition 1","condition 1","condition 1","condition 1","condition 1","condition 2", "condition 2","condition 2","condition 2","condition 2","condition 2"))
type = factor(rep("single-read",12))
design = data.frame(condition,type,row.names=colnames(SimulRNASEQ))
singleSamples = design$type == 'single-read'
countTable = SimulRNASEQ[,singleSamples]
conds = design$condition[singleSamples]

### à partir du dataframe, on peut stocker le jeu de donnée dans la classe propre à DESeq
cds <- newCountDataSet(SimulRNASEQ, conds)
## Estimation des size factors (nécessaire pour l'analyse DEG)
cds = estimateSizeFactors(cds)
# établir la dispersion (nécessaire pour l'analyse DEG)
cds = estimateDispersions(cds)

# Récupérer le jeu de données normalisé :
# Normalized = True : Divise chaque valeur par le sizeFactor associé
counts(cds, normalized = TRUE)

# Analyse DEG : 
res = nbinomTest(cds, condA = "condition 1", condB = "condition 2")
res
# most downregulated :
head( res[ order( res$foldChange, -res$baseMean ), ] )
# most upregulated
head( res[ order( -res$foldChange, -res$baseMean ), ] )

################# ################# #################  DESeq2

### On récupère le jeu de données :
set.seed(20210403) #permet d'obtenir toujours les memes donnees
counts <- simulateRnaSeqData (n=1000,m=12)
counts2=counts[sort(row.names(counts)),] 
SimulRNASEQ <- as.data.frame(counts)
for (i in 1:ncol(SimulRNASEQ)) {
  if (i <= ncol(SimulRNASEQ)/2){
    names(SimulRNASEQ)[i] = paste("control",i,sep='_')    
  }else{
    names(SimulRNASEQ)[i] = paste("test",i-ncol(SimulRNASEQ)/2,sep="_")
  }
}

SimulRNASEQ = read.csv("~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", header = TRUE,row.names = 1)

# Besoin d'énnoncer certains paramètres pour importer dans la classe propore à DESeq2 :
conditions <-factor(c("condition1","condition1","condition1","condition1","condition1","condition1","condition2", "condition2","condition2","condition2","condition2","condition2"))
type = factor(rep("single-read",12))
colData = data.frame(conditions,type,row.names=colnames(SimulRNASEQ))
# Importation du jeu de données dans la classe DESeqDataSet (nécessaire à l'analyse)
dds<-DESeqDataSetFromMatrix(SimulRNASEQ,colData,design=~conditions)

################################ Normalisation + analyse DEG :
# Normalisation = Negative binomiale

## Wald test for GLM coefficients :
#--------------------------------

#Normalisation + calcul des statistiques
DEG_Wald <- DESeq(dds, test = "Wald")
# Affichage des résultats DEG
results(DEG_Wald)

## Likelihood ratio test (LRT) Chi2 pour GLM :
#---------------------------------------------
#Normalisation + calcul des statistiques
DEG_LRT <- DESeq(dds, test = "LRT",reduced = ~1)
# Affichage des résultats DEG
results(DEG_LRT)
