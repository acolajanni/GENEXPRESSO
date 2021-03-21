# Nom               : SimulateNanostring.R
# Type              : Programme
# Objet             : Simule des données de type Nanostring
# Input             : None
# Output            : NANO.csv
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 14.03.2021
#______________________________________________________________________________

#Importation des packages
library('NanoStringNorm')

#Le data set de depart est un jeu de donnees simulees

Nano = NanoString.mRNA
head(Nano)

# Poser la question des gènes ==> Quel différence dans les classes "Positive" / "Negative" / "Endogeneous" 
#Quantifiation des niveaux de gènes : gènes témoins = Positive/négative (= gènes dont on sait qu'il sont faiblement ou fortement exprimé / Endogeneous : gènes qui nous intéressent) ??
# à voir quelle colonne garder là dedans, peut-être que l'on peut se débarasser des num accessions par ex
# on aurait donc en ligne les gènes et en colonne les individus