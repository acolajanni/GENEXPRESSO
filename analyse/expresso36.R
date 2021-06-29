parameters = c(#"RMA.invariantset.medianpolish" ,
               #"RMA.constant.medianpolish" ,
               #"RMA.robust",
               #"RMA.invariantset.mas" ,
               #"RMA.constant.mas" ,
               #"RMA.quantile.mas" ,
               #"RMA.invariantset.liwong" ,
               #"RMA.constant.liwong" ,
               #"RMA.quantile.liwong" ,
               #"MAS.quantile" ,
               #"MAS5" ,
               #"MAS.invariantset",
               #"MAS2.quantile.medianpolish" ,
               #"MAS.medianpolish" ,
               #"MAS2.invariantset.medianpolish" ,
               #"MAS2.quantile.liwong" ,
               #"MAS.liwong" ,
               #"MAS2.invariantset.liwong" ,
               #"MASpm.quantile.mas" ,
               #"MASpm.constant.mas" ,
               #"MASpm.invariantset.mas" ,
               #"MASpm.quantile.medianpolish" ,
               #"MASpm.constant.medianpolish" ,
               #"MASpm.invariantset.medianpolish" ,
               #"MASpm.quantile.liwong" ,
               #"MASpm.constant.liwong" ,
               "MASpm.invariantset.liwong",
               "rma",
               "gcrma",
               "mas5",
               "liwong")


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


library(GENEXPRESSO)
library(affy)
library(gcrma)
library(hgu133plus2.db)
library(dplyr)


list.expresso = NULL
i = i-1
nb.param = 31
for (param in parameters){
  i = i + 1
  message("Normalisation : ",param)
  message("Normalisation ", i,  " sur ", nb.param)
  tmp = tools.norm.Microarray(GEOiD = abatch, 
                              tools = param)
  
  list.expresso[[param]] = tmp
  message("Normalisation : ",param," done !")
  }

save(list.expresso, file = "./data/listexpresso.RData")
#load(file = "./data/listexpresso.RData")
### log 2 transform values :

methods = names(list.expresso)
need.log2 = methods[! methods%in% c("RMA.invariantset.medianpolish", "RMA.constant.medianpolish",
                                    "RMA.robust", "MAS2.quantile.medianpolish",
                                    "MAS.medianpolish","MAS2.invariantset.medianpolish",
                                    "MASpm.quantile.medianpolish","MASpm.constant.medianpolish",
                                    "MASpm.invariantset.medianpolish", "rma", "gcrma")]

for (method in need.log2){
  list.expresso[[method]] = log2(list.expresso[[method]])
}

save(list.expresso, file = "./data/listexpressoLOGGED.RData")

gcrma = as.data.frame(list.expresso$gcrma)
gcrma = mapping.affymetrix.probe(gcrma)

rma = as.data.frame(list.expresso$rma)
rma = mapping.affymetrix.probe(rma)

save(gcrma, file = "./data/gcrmaMAPPED.RData")
save(rma, file = "./data/rmaMAPPED.RData")

### Mapping + DEG
Expresso.comp.methods = function(dataset, T1, T2, DEG.method){
  
  data = as.data.frame(dataset)
  data = mapping.affymetrix.probe(data) #, hgu133plus2SYMBOL)
  

  data = relocate(data, T1,T2)
  DEG = tools.DEG.Microarrays(data, DEG.method, length(T1), length(T2) )
  return(DEG)
}


# listing parameters of expresso function
params = names(list.expresso)
params = parameters

# initializing the dataframe to compare methods
data.to.comp = data.frame(NULL)

for (param in params){
  print(param)

  tmp.dataset = list.expresso[[param]]
  
  tmp = Expresso.comp.methods(tmp.dataset, T1, T2, "RankSum.log")

  colnames(tmp) = c( paste(param,"Up"), paste(param,"Down"), "SYMBOL" )
  
  if( dim(data.to.comp)[1] == 0){
    data.to.comp = tmp
  }
  else{
    data.to.comp = merge(data.to.comp, tmp, by = "SYMBOL", all=T )
  }
  
}
row.names(data.to.comp) = data.to.comp$SYMBOL
data.to.comp$SYMBOL = NULL

save(data.to.comp, file = "./data/expresso31.RData")

data.to.comp = as.data.frame(t(data.to.comp))
pca = PCA_tools(as.data.frame(t(data.to.comp)))
pca
