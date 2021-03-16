# Nom               : test_norm_microarray
# Type              : Programme
# Objet             : Tester la normalisation (les fonctions iront dans functions.R)
# Input             : dataset micro array
# Output            : ?
# Auteur            : Antonin COLAJANNI
# R version         : 3.6
# Date de creation  : 16.03.2021
#______________________________________________________________________________

library(limma)
library(madsim)

mydata1 <- madsim(mdata = NULL, n = 100, ratio = 0, fparams, dparams, sdn, rseed)
micro <- as.data.frame(mydata1$xdata)
micro

data = read.csv("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv",header = T, sep = ',')
data
targets = readTargets("~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv", row.names= 1, sep = ',')
targets



help(readTargets)
help(read.maimages)

RG = read.maimages(targets, source = 'agilent', wt.fun = F)


help("madsim")

