# Nom               : geolimma.R
# Type              : Programme
# Objet             : test l'analyse DEG de geolimma
# Input             : Dataset normalisé de Nanostring + dataset microarray
# Output            : analyse DEG 
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 14.04.2021
#______________________________________________________________________________

source(file.path("./SCRIPT","functions.R"))


# Problème 1 : GEOlimma ??? Gene Expression Omnibus - Articles récents mais pas de package dispo
Micro = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv", header = TRUE,row.names = 1)

Nano = read.csv("~/GIT/CPRD/DATA/NANOSTRING/NanoNormNappa.csv", header = TRUE,row.names = 1)

#________________________________________________ 
#Limma : 
# Nano : Déja implémenté

#________________________________________________ 
#Limma : 
# Microarray : (extrait de functions.R)

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

#Récupération des pvalues
res.diff_limma = DEG_limma(Micro,design)
head(res.diff_limma)
tail(res.diff_limma)
#________________________________________________
# GEOlimma sur Micro array

# construction de la matrice de design
design = make_designMatrix(dataset = Micro)

# récupération des pvalues
res.diff_GEOlimma= DEG_GEOlimma(Micro,design)
head(res.diff_GEOlimma)
tail(res.diff_GEOlimma)



  
