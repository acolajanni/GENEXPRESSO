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



# Méthode par défaut : TMM
TMM_norm <- calcNormFactors.DGEList(data, method = "TMM")
TMM_norm

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
#AveLogCPM ==> + ou - la meme chose
# De cette manière on peut renvoyer un jeu de données normalisé analysable avec d'autres outils que edgeR

# on affiche par classement, les gènes les plus différentiellement exprimés
topTags(DEG, sort.by = 'PValue')


tools_norm_RNAseq.inspect <- function(raw.data,tool){
  data = as.matrix(raw.data)
  Norm_factors = data.frame("Samples" = colnames(data))
  
  tools_norm_RNAseq.fnc <- switch(tool,
                                  
        edgeR_TMM = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          TMM_norm <- calcNormFactors.DGEList(DGE, method = "TMM")
          #tmp = data.frame("TMM_norm" = TMM_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          res.norm = cpm(TMM_norm)
        },
        
        edgeR_RLE = {
           DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
           RLE_norm <- calcNormFactors.DGEList(DGE, method = "RLE")
           #tmp = data.frame("RLE_norm" = RLE_norm$samples$norm.factors)
           #Norm_factors = cbind(Norm_factors, tmp)
           res.norm = aveLogCPM(RLE_norm)
        },
        
        edgeR_TMMwsp = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          TMM_wsp <- calcNormFactors.DGEList(DGE, method = "TMMwsp")
          #tmp = data.frame("TMMwsp" = RLE_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          res.norm = aveLogCPM(TMM_wsp)
        },
        
        edgeR_upperquartile = {
          DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
          Upper <- calcNormFactors.DGEList(DGE, method = "upperquartile")
          #tmp = data.frame("upperquartile" = RLE_norm$samples$norm.factors)
          #Norm_factors = cbind(Norm_factors, tmp)
          res.norm = aveLogCPM(Upper)
        }
        
        )

  return(res.norm)
}
# La fonction ci dessus retoure le jeu de données normalisé( ==> faire une fonction qui reprend la suite : calcul DEG)
data_to_comp = tools_norm_RNAseq.inspect(data,tool = 'edgeR_TMM')

################################################################ Pour l'instant ça c'est pas pas bon (pas besoin de merge les jeux de données)
# Réutiliser les bouts de code pour appliquer les calculs DEG

tools = c("edgeR_RLE","edgeR_upperquartile","edgeR_TMMwsp")
for (tool in tools){
  print(tool)
  tmp = tools_norm_RNAseq.inspect(data,tool)
  data_to_comp = merge(data_to_comp,tmp,by = "Samples",all=T)  
}
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
data = DGEList(count = data, 
               group = rep(1:2,each=ncol(data)/2))


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
topTags(qlf)
topTags(qlf2)

### Au final : on peut appliquer le GLM sur n'importe quel facteur de normalisation si c'est une DEGlist
qlf_tmm <- glmQLFTest(fit_tmm, contrast=BvsA)
qlf <- glmQLFTest(fit, contrast=BvsA)
topTags(qlf_tmm)
topTags(qlf)

