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
#head(data_to_comp)

data_to_comp_Up = copy(data_to_comp)
data_to_comp_Down = copy(data_to_comp)

for (i in ncol(data_to_comp):1){
  if (i%%2 == 0){
    # upregulated : condition A < condition B (donc upregulated = wilcox.test(A,B,"Less")) 
    # Downregulated : condition A > condition B (si on prend A comme référence)
    data_to_comp_Up = data_to_comp_Up[,-i]
  }else{
    data_to_comp_Down = data_to_comp_Down[,-i]
  }
}

#head(data_to_comp_Down)
#head(data_to_comp_Up)


#on récupère le tableau semblable à MakeComparisonTable.R
data_to_comp_Up = as.data.frame(t(data_to_comp_Up))
data_to_comp_Down = as.data.frame(t(data_to_comp_Down))

## PCA : 
#data_to_comp = t(data_to_comp)
#PCA_tools(data_to_comp)
PCA_tools(data_to_comp_Up)
PCA_tools(data_to_comp_Down)

## Upsetplot : si erreur (pas d'intesection car très peu de gènes DE)
par(mfrow=c(1,2))
A = UpsetPlot(data.to.comp = data_to_comp_Down,threshold = 0.05)
names = colnames(A)

upset(A, sets = names, sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on" )

B = UpsetPlot(data.to.comp = data_to_comp_Up,threshold = 0.05)
names = colnames(B)

upset(B, sets = names, sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = NULL )



#UpsetPlot(data.to.comp = data_to_comp_Up,threshold = 0.05)
