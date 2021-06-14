library(GENEXPRESSO)
## Jeu de données
# save(abatch, file = "abatch.RData")
load("./data/abatch.RData")
## PCFIXE
celpath = "./data"
## Cremi
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


?tools.norm.Microarray

######################
normalize.AffyBatch.methods()
bgcorrect.methods()
pmcorrect.methods()
express.summary.stat.methods()
###################### 

#---------------------------#
# Comparaisons des méthodes #  
#---------------------------#

#### BG correct
# Fixation des méthodes autres à :
norm = "constant"
pm.cor = "mas"
stat = "mas"
bg.cor = c("mas","rma")

DEG.method = "RankSum.log"

dataset = NULL
for (bg in bg.cor){
  tmp = tools.norm.Microarray(GEOiD = abatch, tools = "custom", 
                              tools.normalize = norm,
                              tools.bgcorrect = bg,
                              tools.pmcorrect = pm.cor,
                              tools.express.summary.stat = stat)
  
  
  
  dataset[[bg]] = tmp
}

dataset.bg = dataset


#### PM correct
pm.cor = c("mas","pmonly","subtractmm")
bg.cor = "mas"

dataset.pm = NULL
pm.cor = "subtractmm"
for (pm in pm.cor){
  tmp = tools.norm.Microarray(GEOiD = abatch, tools = "custom", 
                              tools.normalize = norm,
                              tools.bgcorrect = bg.cor,
                              tools.pmcorrect = pm,
                              tools.express.summary.stat = stat)
  
  
  
  dataset.pm[[pm]] = tmp
}

#### SumStat method
sumstat = c("avgdiff","mas","medianpolish","playerout")
norm = "constant"
pm.cor = "mas"
bg.cor = "mas"
stat = "liwong"
dataset.sumstat = NULL
for (stat in sumstat){
  tmp = tools.norm.Microarray(GEOiD = abatch, tools = "custom", 
                              tools.normalize = norm,
                              tools.bgcorrect = bg.cor,
                              tools.pmcorrect = pm.cor,
                              tools.express.summary.stat = stat)
  
  
  dataset.sumstat[[stat]] = tmp
}

#### Normalization method
pm.cor = "mas"
bg.cor = "mas"
stat = "mas"
norm = c("constant","contrasts","invariantset","loess", "qspline", "quantiles", "quantiles.robust")
norm = c("contrasts",
         "invariantset"
         ,"loess"
         , "qspline"
         , "quantiles"
         , "quantiles.robust")

dataset.norm = NULL
for (n in norm){
  tmp = tools.norm.Microarray(GEOiD = abatch, tools = "custom", 
                              tools.normalize = n,
                              tools.bgcorrect = bg.cor,
                              tools.pmcorrect = pm.cor,
                              tools.express.summary.stat = stat)
  
  
  dataset.norm[[n]] = tmp
}

save( dataset.bg, dataset.norm, dataset.pm, dataset.sumstat, file = "./data/NormCompTotal.RData")

############################################################
bg.mas = as.data.frame(dataset.bg$mas)
bg.rma = as.data.frame(dataset.bg$rma)

map.bg.mas = mapping.affymetrix.probe(bg.mas, hgu133plus2SYMBOL)
map.bg.mas = relocate(map.bg.mas, T1,T2)

map.bg.rma = mapping.affymetrix.probe(bg.rma, hgu133plus2SYMBOL)
map.bg.rma = relocate(map.bg.rma, T1,T2)

DEG.method = "RankSum"
DEG.bg.rma = tools.DEG.Microarrays(map.bg.rma, DEG.method, length(T1), length(T2) )

############################################################


  # un jeu de donnée + le nom de la méthode utilisée + T1 et T2 (vecteur des noms de colonnes appartenant au grp 1 / 2)
Expresso.comp.methods = function(dataset, T1, T2, DEG.method){
  
  data = as.data.frame(dataset)
  data = mapping.affymetrix.probe(data) #, hgu133plus2SYMBOL)
  data = relocate(data, T1,T2)
  
  DEG = tools.DEG.Microarrays(data, DEG.method, length(T1), length(T2) )
  return(DEG)
  
}
#test.fnc = Expresso.comp.methods(test1, T1, T2, "RankSum")

# Suppression of NA filled dataset
dataset.pm$subtractmm = NULL
# One object with all datasets
dataset.list = list(background = dataset.bg, pm.cor = dataset.pm, sumstat = dataset.sumstat, norm = dataset.norm)
# Counting all the existing dataset
nb.dataset = sum(sapply(dataset.list,length))
# listing parameters of expresso function
params = names(dataset.list)
params = "norm"


# initializing the dataframe to compare methods
data.to.comp = data.frame(NULL)

library(dplyr)
library(hgu133plus2.db)

params="background"
tmp.params = "rma"

for (param in params){
  print(" ....... ")
  print(param)
  
  tmp.list = dataset.list[[param]]
  #tmp.params = names(tmp.list)
 
  for (tmp.param in tmp.params){
    print(tmp.param)
    
    tmp.dataset = tmp.list[[tmp.param]]
    if (tmp.param == "medianpolish"){
      tmp = Expresso.comp.methods(tmp.dataset, T1, T2, "RankSum.log")
    }else{
      tmp = Expresso.comp.methods(tmp.dataset, T1, T2, "RankSum")
    }
    
    
    method = paste(param, tmp.param)
    colnames(tmp) = c( paste(method,"Up"), paste(method,"Down"), "SYMBOL" )
    
    if( dim(data.to.comp)[1] == 0){
      data.to.comp = tmp
    }
    else{
      data.to.comp = merge(data.to.comp, tmp, by = "SYMBOL", all=T )
    }
  }
}


row.names(data.to.comp) = data.to.comp$SYMBOL
data.to.comp$SYMBOL = NULL
data.to.comp = as.data.frame(t(data.to.comp))
pca = PCA_tools(data.to.comp)
pca

load("./dataToComp.RData")
# Séparation du jeu de données en up/down regulated
data.to.comp = as.data.frame(t(data.to.comp))
methods = row.names(data.to.comp)

meth.bg.cor = methods[grepl("background",methods)]
meth.bg.cor = meth.bg.cor[grepl("Up",meth.bg.cor)]

meth.bg.up = data.to.comp[meth.bg.cor,]

Upreg = data.to.comp[grepl("Up|less",methods)]
Downreg = data.to.comp[grepl("Down|greater",methods)]

# PCA pour UP et Down
Upreg = as.data.frame(t(Upreg))
Downreg = as.data.frame(t(Downreg))

pca_up = PCA_tools(t(Upreg))
pca_down = PCA_tools(Downreg)
pca_down


#save(data.to.comp, file = "./dataToComp.RData")
#load("./dataToComp.RData")

library(UpSetR)
upsetDown = Upset.Binary.Dataframe(Downreg)
methods = row.names(upsetDown)
upset(upsetDown, 
      sets = methods, 
      sets.bar.color = "#56B4E9", 
      order.by = "freq"
      )

###########################
########  Heatmap  ########
###########################
library(stringr)
load("./dataToComp.RData")
load("./data/NormCompTotal.RData")
data.to.comp = as.data.frame(t(data.to.comp))

# Liste de toute les méthodes
methods = row.names(data.to.comp)
# Dataframe rempli de valeur binaire (0/1)
upset = Upset.Binary.Dataframe(data.to.comp)

#upset.bg =

upsetUnion = Get.DEG.2(upset, alternative=FALSE, method = "union")
upsetInter = Get.DEG.2(upset, alternative=TRUE, method = "intersect")


UpsetGenes.intersect = c(upsetInter$Upregulated, upsetInter$Downregulated)
UpsetGenes.intersect = data.frame(DEG = UpsetGenes.intersect)

UpsetGenes.union = c(upsetUnion$Upregulated, upsetUnion$Downregulated)

save(UpsetGenes.intersect, UpsetGenes.union, file = "./data/inter&union.RData")


write.csv(UpsetGenes.intersect, file = "./data/GSE31684_Intersect_genes.csv",row.names = FALSE)






Get.DEG <- function(binary.matrix, method, alternative){
  # By default argument : 
  # We get the results for both upregulated and down regulate genes
  if (missing(alternative)){
    alternative = TRUE
  }
  
  if (alternative){
    # Series of logical value wether or not the "method" argument is contained in the names of columns 
    pattern = str_detect(colnames(binary.matrix),method)
    if (!TRUE%in% pattern){
      stop("Unknown method")
    }
    # actual colnames corresponding to the pattern of "method" is extracted
    # example : 
    # method = "DESeq2"
    # method.name = "DESeq2_Up", "DESeq2_Down"
    method.name = colnames(binary.matrix)[pattern]
  }
  else{
    method.name = method[method %in% colnames(binary.matrix)]
    
    if (length(method.name) != length(method)){
      notIN = method[!method %in% colnames(binary.matrix)]
      warning("One or more method doesn't match : ", notIN)
    }
  }
  
  #######################################################
  
  binary.matrix = as.data.frame(binary.matrix)
  # Removing all the other methods
  result = binary.matrix[ , colnames(binary.matrix) %in% method.name ]
  result$genes = row.names(result)
  
  
  #return(result)
  
  #######################################################
  
  # if alternative hypothesis does not interest us
  if (!alternative){
    # Only one column (given in argument) is retrieved
    DEG = subset(result, result[[method.name]] == 1)
    DEG = list(DEG = DEG[["genes"]])
  }
  else{
    # If we want alternative hypothesis, two dataframe are needed (up and downregulated)
    DEG_up = data.frame(Upregulated = result[[ method.name[1] ]], genes = result$genes)
    # only genes with "1" are differentially expressed, other does not interest us
    DEG_up = filter(DEG_up, Upregulated == 1 )
    
    DEG_down = data.frame(Downregulated = result[[ method.name[2] ]], genes = result$genes)
    DEG_down = filter(DEG_down, Downregulated == 1 )
    
    # A list of 2 vectors is returned : for both up and downregulated hypothesis
    DEG = list(Upregulated = DEG_up$genes , Downregulated = DEG_down$genes)
  }
  return(DEG)
}






