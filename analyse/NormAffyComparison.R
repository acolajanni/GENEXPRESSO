library(GENEXPRESSO)
## Jeu de données
# save(abatch, file = "abatch.RData")
load("./data/abatch.RData")
txt.dir = paste0(celpath,"/GSE31684_table_of_clinical_details.txt")
tab = read.delim(txt.dir,check.names=FALSE,as.is=TRUE, header = T)
samples = subset(tab, tab$PreOpClinStage == 'T1' | tab$PreOpClinStage == 'T2') 



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

dataset.norm = NULL
for (n in norm){
  tmp = tools.norm.Microarray(GEOiD = abatch, tools = "custom", 
                              tools.normalize = n,
                              tools.bgcorrect = bg.cor,
                              tools.pmcorrect = pm.cor,
                              tools.express.summary.stat = stat)
  
  
  dataset.norm[[n]] = tmp
}



library(hgu133plus2.db)
library(dplyr)
map = as.data.frame(dataset$mas)
map2 = as.data.frame(dataset$rma)
map = mapping.affymetrix.probe(map, hgu133plus2SYMBOL )
map2 = mapping.affymetrix.probe(map2, hgu133plus2SYMBOL )

