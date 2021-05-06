# Nom               : Main.R
# Type              : Programme
# Objet             : Proposition de programme principal appelant tous les autres
# Input             : 
# Output            : 
# Auteur            : Juliette Casemajor
# R version         : 3.6
# Date de creation  : 23.04.21
#__________________________________________________________________Environnement
#Chemin du projet
.ROOT <- file.path("~/GIT/CPRD/")
#__________________________________________________________________________DEBUT
#Chargement des chemins par defaut
.PROGRAMS <- file.path(.ROOT,"SCRIPT")
setwd(.PROGRAMS)
source("setup.R")
source("functions.R")
#___________________________________________________Normalisation et analyse DEG
#Chargement des donnees brutes
.DATA <- file.path(.ROOT,"DATA/MICROARRAYS")
setwd(.DATA)
DataIn <- read.csv("Simulmicroarrays1000.csv", header = TRUE, row.names = 1)

#Type de donnees : "nanostring", "microarrays" ou "rna-seq"
type <- "microarrays"

#Indiquer nombre d'individus du groupe 1 et du groupe 2
n1 <- 6
n2 <- 6

#Pour nanostring et rna seq, outils de normalisation : mettre liste
#Outils normalisation : 
# microarrays : normalisation intégrée dans le package de production de données
# rna seq : "EdgeR", "DESeq" et "DESeq2"
# nanostring : "NanostringNorm", "NAPPA", "NanostringR", "NanoR" et "DESeq2"
norm <- "a completer"
#Outils DEG : 
# microarray : "Wilcox, "RankProduct.param1","RankProduct.param2","RankProduct.param3","RankProduct.param4", "GEOlimma" et/ou "limma"
# rna seq : "EdgeR", "DESeq" et "DESeq2"
# nanostring : "NanostringR", "nappa.NS", "nappa.default", "nappa.param1", "nappa.param2", "nappa.param3", "nanostringnorm.default", "nanostringnorm.param1", "nanostringnorm.param2",
#              "desq2", "nanostringDiff", "nanoR.top100" et/ou "nanoR.total"
tools = c("GEOlimma", "Wilcox","RankProduct.param1","RankProduct.param2","RankProduct.param3","RankProduct.param4")

setwd(.PROGRAMS)

if (type == "microarrays") {
  #Pas de normalisation, execution du programme d'analyse de microarrays
  source("testmainmicroarrays.R") #Je l'ai rajoute dans les scripts, a supprimer plus tard bien sur
  } else if (type == "nanostring") {
  #script d'analyse nanostring avec les variables deja implementees ici
  #le fichier source d'analyse ne devra contenir que norm, DataIn, tools, etc
  } else if (type == "rna-seq") {
  #Idem
}

#NB : les seules modifications que l'utilisateur devra faire seront dans ce fichier
#Elements a modifier par l'utilisateur :
#emplacement des fichiers (.ROOT, .PROGRAMS, .DATA)
#Fichiers de raw data (DataIn)
#Type de donnees (type)
#Nombre d'individus (n1 et n2)
#Outils de normalisation et de DEG utilisés (norm et tools)

#__________________________________________________________________Visualisation
#Peut être que la visualisation pourrait être integree dans l'analyse DEG ou 
#Bien la sortir pourrait être utile à l'utilisateur s'il veut faire lui-même
#sa visualisation a partir des p-value ?
#__________________________________________________________________Co-expression

#_________________________________________________________________________Output
#Enregistrement des fichiers .csv et img dans un dossier output
