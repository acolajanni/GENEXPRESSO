parameters = c("RMA.invariantset.medianpolish" ,
               "RMA.constant.medianpolish" ,
               "RMA.robust",
               "RMA.invariantset.mas" ,
               "RMA.constant.mas" ,
               "RMA.quantile.mas" ,
               "RMA.invariantset.liwong" ,
               "RMA.constant.liwong" ,
               "RMA.quantile.liwong" ,
               "MAS.quantiles" ,
               "MAS5" ,
               "MAS.invariantset",
               "MAS2.quantile.medianpolish" ,
               "MAS.medianpolish" ,
               "MAS2.invariantset.medianpolish" ,
               "MAS2.quantile.liwong" ,
               "MAS.liwong" ,
               "MAS2.invariantset.liwong" ,
               "MASpm.quantile.mas" ,
               "MASpm.constant.mas" ,
               "MASpm.invariantset.mas" ,
               "MASpm.quantile.medianpolish" ,
               "MASpm.constant.medianpolish" ,
               "MASpm.invariantset.medianpolish" ,
               "MASpm.quantile.liwong" ,
               "MASpm.constant.liwong" ,
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


list.expresso = NULL
i = 0
nb.param = length(parameters)
for (param in parameters){
  print(param)
  tmp = tools.norm.Microarray(GEOiD = abatch, 
                              tools = param)
  
  list.expresso[[param]] = tmp
  i = i + 1
  
  message("Normalisation ", i, " sur ", nb.param)
  }



