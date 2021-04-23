data_to_comp = tools.microarrays.inspect(DataIn,"limma", n1, n2)
for (tool in tools){
  print(tool)
  tmp = tools.microarrays.inspect(DataIn,tool,n1, n2)
  data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
}

row.names(data_to_comp) <- data_to_comp$Gene.ID
data_to_comp <- data_to_comp[,-1]

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

#Recuperation du tableau
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