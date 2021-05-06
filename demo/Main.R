###########################
###   Demonstration   #####
###########################

# Tout d'abord, il faut choisir un type de données parmis les 3 suivants :
type = "RNAseq"
type = "Microarrays"
type = "Nanostring"

# Un jeu de donnée va être simulé pour Microarrays et RNAseq, 
# Il sera récupéré dans le dossier /DATA pour Nanostring
data = Simul.data(type,n.cond1 = 15, n.cond2 = 15,nb.genes = 1000)

# Afin d'analyser les jeux de données pour trouver les gènes étant différentiellement exprimé,
# il faut d'abord traiter les données. On va les "normaliser".
# Elle pourra s'appliquer pour Nanostring et RNAseq.
# Les données Microarrays sont simulées telles qu'il est inutile de les normaliser

# on va utiliser un certain nombre de technique de traitement de données propre à chacune des technique
# Ces méthodes seront stockées dans "tools"
if (type == "RNAseq"){
  # La particularité des pacakges d'analyse de RNAseq est qu'ils utilisent des normalisations
  # intégrés dans l'analyse statistique. C'est pour cela qu'un seul objet "tools" est créer
  tools = c("edgeR_TMM", "deseq2.Wald","deseq2.LRT","deseq")
  data.to.comp = tools.DEG.RNAseq.merge(data,tools)
}else if(type == "Microarrays"){
  # Un seul objet tools est créer car seul l'analyse statistique nous intéresse
  tools = c("limma", "Wilcox","RankProduct","RankSum")
  data.to.comp = tools.DEG.Microarrays.merge(data,tools,n1=15,n2=15)
  
}else if(type == "Nanostring"){
  tools_DEG = c("limma","RankProduct","RankSum" )
  tools_norm = c("nappa.NS","nanoR.total")
  RCC.dir <- file.path("./DATA/NANOSTRING","GSE146204_RAW")
  data.to.comp = tools.DEG.Nanostring.merge()
  res.DEG = tools.DEG.Nanostring.merge(data,tools_DEG,tools_norm,DESeq=T,dir = RCC.dir)
}

# Selon le type choisi, un tableau contenant les pvalues qu'un gène soit différentiellement exprimé pour la RNAseq
# Et Sur/sous exprimé dans la deuxième condition expérimental pour les données Nanostring et Microarrays.

# Ces tableaux pourront permettre de faire un graphique d'analyse en composante principale pour comparer les méthodes employées
PCA = PCA_tools(data.to.comp)

# Une seconde information pourra être apportée par les upsets plots :
# Le nombre de gènes repéré comme différentiellement exprimé en communs entre chaque méthode


