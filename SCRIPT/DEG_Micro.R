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


tools.microarrays.inspect <- function(data,tool,n1,n2){
  
  design = make_designMatrix(dataset = data)

  DEG_Microarrays_tools.fnc <- switch(tool,
                                      GEOlimma = {
                                        res.diff_up= DEG_GEOlimma(data,design, comp = "up")
                                        res.diff_down= DEG_GEOlimma(data,design, comp = "down")
                                        res.diff = merge(res.diff_up,res.diff_down,by = "Gene.ID",all=T)
                                      },
                                      
                                      limma = {
                                        res.diff_up = DEG_limma(data,design, comp = "up")
                                        res.diff_down = DEG_limma(data,design, comp = "down")
                                        res.diff = merge(res.diff_up,res.diff_down,by = "Gene.ID",all=T)
                                      },
                                      
                                      Wilcox = {
                                        res.diff = wilcoxDEG4(micro, n1, n2)
                                        res.diff$Gene.ID = row.names(res.diff)
                                        res.diff = res.diff[,-1]
                                      },
                                      
                                      stop("Enter something that switches me!")
  )
  return(res.diff)
}

data_to_comp = tools.microarrays.inspect(micro,"limma", n1 = 6, n2 = 6)
tools = c("GEOlimma", "Wilcox")
for (tool in tools){
  print(tool)
  tmp = tools.microarrays.inspect(micro,tool,n1 = 6, n2 = 6)
  data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
}
row.names(data_to_comp) <- data_to_comp$Gene.ID
data_to_comp <- data_to_comp[,-1]
head(data_to_comp)

data_to_comp_Up = copy(data_to_comp)
data_to_comp_Down = copy(data_to_comp)

for (i in ncol(data_to_comp):1){
  if (i%%2 == 0){
    # uprelated : condition B > condition A (donc upregulated = wilcox.test(A,B,"Less"))
    # Downregulated : condition B < condition A (si on prend A comme référence)
    data_to_comp_Up = data_to_comp_Up[,-i]
  }else{
    data_to_comp_Down = data_to_comp_Down[,-i]
  }
}

head(data_to_comp_Down)
head(data_to_comp_Up)

#on récupère le tableau semblable à MakeComparisonTable.R
data_to_comp_Up = as.data.frame(t(data_to_comp_Up))
data_to_comp_Down = as.data.frame(t(data_to_comp_Down))

## PCA : erreur (que deux méthodes)
PCA_tools(data_to_comp_Up)
PCA_tools(data_to_comp_Down)

## Upsetplot : erreur (pas d'intesection car très peu de gènes DE)
UpsetPlot(data.to.comp = data_to_comp_Down,threshold = 0.05)

