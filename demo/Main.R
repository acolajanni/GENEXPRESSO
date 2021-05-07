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
if (type == "Nanostring"){
  raw.data = data
  data = data$rcc.df
}

# Afin d'analyser les jeux de données pour trouver les gènes étant différentiellement exprimé,
# il faut d'abord traiter les données. On va les "normaliser".
# Elle pourra s'appliquer pour Nanostring et RNAseq.
# Les données Microarrays sont simulées telles qu'il est inutile de les normaliser

###############
# Analyse DEG #
###############

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
  tools_norm = c("nappa.NS","nanostringnorm.default")
  data.to.comp = tools.DEG.Nanostring.merge()
  res.DEG = tools.DEG.Nanostring.merge(data,tools_DEG,tools_norm,DESeq=T)
}

### Visualisation

# Selon le type choisi, un tableau contenant les pvalues qu'un gène soit différentiellement exprimé pour la RNAseq
# Et Sur/sous exprimé dans la deuxième condition expérimental pour les données Nanostring et Microarrays.

# Ces tableaux pourront permettre de faire un graphique d'analyse en composante principale pour comparer les méthodes employées
PCA = PCA_tools(data.to.comp)

# Une seconde information pourra être apportée par les upsets plots :
# Le nombre de gènes repéré comme différentiellement exprimé en communs entre chaque méthode
# On construit d'abord un matrice constitutée de 0 et de 1 :
# 1 pour les gènes différentiellement exprimés, 0 pour les autres
Binary.Matrix = Upset.Binary.Dataframe(data.to.comp)
# On peut donc maintenant construire l'upset plot
upset(Upset, sets = names(Upset), sets.bar.color = "#56B4E9", order.by = "freq",empty.intersections = NULL )



########################
# Analyse coExpression #
########################

if (type == "Nanostring"){
  # Les gènes servant de contrôles positifs et négatifs ne nous intéresse pas pour cette étape
  # Ils ne sont utile que lors de la Normalisation
  # Ainsi on ne va garder que les gènes à l'origine de signaux "Endogènes"
  data <- as.data.frame(data[grepl("Endog",raw.data$annots.df$CodeClass),])
}

# On va d'abord récupérer les données de corrélations et/ou de similarité par paire de gène :
# De manière rapide :
Relations = Make.full.adjacency(data,PValue = FALSE)
# De manière (beaucoup) plus lente :
Relations = Make.full.adjacency(data,PValue = TRUE)

# Ces informations sont récupérés pour l'intégralité du jeu de données.

# On peut comparer les différentes méthodes entre elles :

# Les graphes sont contruits avec les trois méthodes implémentées
spearman = Make.df.graph(data, cor.threshold = 0.8,Pvalue.threshold = F ,method = "spearman")
TOM = Make.df.graph(data, cor.threshold = 0.5,method = "TOM")
kendall = Make.df.graph(data,  cor.threshold = 0.65,Pvalue.threshold = F,method = "kendall")
# Affichage en comparaisondeux à deux des trois méthodes 
par(mfrow = c(1,3))
CompGraph_TOM_kendall = relations.comparison(TOM, kendall, "TOM", "Kendall",color.g1 = "blue", color.g2 = "orange")
CompGraph_TOM_spearman = relations.comparison(TOM, spearman, "TOM", "Spearman",color.g1 = "blue", color.g2 = "darkgreen")
CompGraph_spearman_kendall = relations.comparison(spearman, kendall, "Spearman", "Kendall", color.g1 = "darkgreen", color.g2 = "orange")

# On peut aussi comparer deux à deux les données au sein des deux conditions expérimentales

# On divise le jeu de données en deux groupes correspondant aux deux conditions expérimentales
if (type == "Nanostring"){
  # On récupère les données des groupes
  group = table(raw.data$samples.IDs$tp53.status)
  n1 = as.integer(group[1])
  n2 = as.integer(group[2])
  
}else {
  # 15 individus dans chaque groupe pour les autres données
  n1 = 15
  n2 = 15
}
# Division du jeu de données en deux groupes
data_gr1 = data[1:n1]
data_gr2 = data[(n1+1):(n1+n2)]

# On peut comparer les deux groupes selon le critère de similarité "TOM" par exemple
TomG1 = Make.df.graph(data_gr1, cor.threshold = 0.15,Pvalue.threshold = F ,method = "TOM")
TomG2 = Make.df.graph(data_gr2, cor.threshold = 0.15,Pvalue.threshold = F ,method = "TOM")
# Afficahge de la comparaison entre les deux groupes
dev.off()
CompGraph_total = relations.comparison(TomG1, TomG2, "groupe 1", "groupe 2")

