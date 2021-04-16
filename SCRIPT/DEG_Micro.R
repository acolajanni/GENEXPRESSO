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
micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", header = TRUE,row.names = 1)


tools.microarrays.inspect <- function(data,tool){
  
  design = make_designMatrix(dataset = Micro)

  DEG_Microarrays_tools.fnc <- switch(tool,
                                      GEOlimma = {
                                        res.diff= DEG_GEOlimma(data,design)
                                      },
                                      
                                      limma = {
                                        res.diff = DEG_limma(data,design)
                                      },
                                      
                                      stop("Enter something that switches me!")
  )
  return(res.diff)
}

data_to_comp = tools.microarrays.inspect(micro,"limma")
tools = c("GEOlimma")
for (tool in tools){
  print(tool)
  tmp = tools.microarrays.inspect(micro,tool)
  data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
}
row.names(data_to_comp) <- data_to_comp$Gene.ID
data_to_comp <- data_to_comp[,-1]
#on récupère le tableau semblable à MakeComparisonTable.R
data_to_comp = as.data.frame(t(data_to_comp))


## PCA : erreur (que deux méthodes)
PCA_tools(data_to_comp)

## Upsetplot : erreur (que deux méthodes)
UpsetPlot(data.to.comp = data_to_comp,threshold = 0.05)

