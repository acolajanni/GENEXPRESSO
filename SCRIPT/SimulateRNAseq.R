# Nom               : simulateRNAseq
# Type              : Programme
# Objet             : Simule des données de type RNAseq
# Input             : None
# Output            : SimulRNASEQ.csv
# Auteur            : Juliette Casemajor / Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 04.03.2021
#______________________________________________________________________________

#Importation des packages

BiocManager::install(c("DEFormats","edgeR"))
library("DEFormats")
#Le data set de depart est un jeu de donnees simulees
#NB : count = nombre de reads / genes

set.seed(20210403) #permet d'obtenir toujours les memes donnees
counts <- simulateRnaSeqData (n=100,m=12)
#n correspond au nombre de genes et m au nombre de patients

#Tri des resultats en fonction du nom des genes (pour faciliter la comparaison)

counts2=counts[sort(row.names(counts)),] 
#attention, comme les noms sont gene1,gene2,gene3... Pas forcement besoin de les
#classer car sinon on a gene1,gene10,gene100,gene11...

head(counts) #affiche le debut du tableau de donnees

#Conversion en data.frame
SimulRNASEQ <- as.data.frame(counts)

#Renommer les colonnes 
# Ne fonctionne comme il faut que si on a un nombre paire de colonnes
for (i in 1:ncol(SimulRNASEQ)) {
  if (i <= ncol(SimulRNASEQ)/2){
    names(SimulRNASEQ)[i] = paste("control",i,sep='_')    
  }else{
    names(SimulRNASEQ)[i] = paste("test",i-ncol(SimulRNASEQ)/2,sep="_")
  }
}

#help("simulateRnaSeqData")

### Plus gros jeu de données
counts <- simulateRnaSeqData (n=1000,m=30, seed = 222)
group = paste0(rep(c("control", "case"), each = 15),rep(c(1:15),each = 1))
colnames(counts) = group
counts = as.data.frame(counts)

#Exportation du jeu de donnees dans ~/GIT/DATA/RNASEQ
write.csv(SimulRNASEQ,"~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ.csv", row.names = TRUE)
write.csv(counts,"~/GIT/CPRD/DATA/RNASEQ/SimulRNASEQ1000x30.csv", row.names = TRUE)
