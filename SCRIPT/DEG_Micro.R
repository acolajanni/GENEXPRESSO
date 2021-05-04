# Nom               : DEG_Micro.R
# Type              : Programme
# Objet             : Analyse DEG microarray
# Input             : dataset Microarray
# Output            : data_to_comp 
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 16.03.2021
#______________________________________________________________________________
source(file.path("./SCRIPT","functions.R"))
source(file.path("./SCRIPT","functionwilcoxtest.R"))



micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)
micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/SimulmicroarraysBIG.csv", header = TRUE,row.names = 1)
micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", header = TRUE,row.names = 1)
#head(micro)

data_to_comp = tools.microarrays.inspect(micro,"limma", n1 = 12, n2 = 12)
tools = c("GEOlimma", "Wilcox","RankProduct.param1","RankProduct.param2","RankProduct.param3","RankProduct.param4")
for (tool in tools){
  print(tool)
  tmp = tools.microarrays.inspect(micro,tool,n1 = 12, n2 = 12)
  data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
}


row.names(data_to_comp) <- data_to_comp$Gene.ID
data_to_comp <- data_to_comp[,-1]
data_to_comp = as.data.frame(t(data_to_comp))
#head(data_to_comp)
## test log10
#data_to_comp = -log10(data_to_comp)


#Compute PCA
PCA_tools(data_to_comp)

## Upsetplot : si erreur (pas d'intesection car très peu de gènes DE)
Upset = UpsetPlot(data.to.comp = data_to_comp,threshold = 0.05)
methods = colnames(Upset)

Upreg = methods[grepl("Up|less",methods)]
Downreg = methods[grepl("Down|greater",methods)]

upset(Upset, sets = Upreg, sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = NULL )


upset(Upset, sets = Downreg, sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = NULL )



#UpsetPlot(data.to.comp = data_to_comp_Up,threshold = 0.05)
