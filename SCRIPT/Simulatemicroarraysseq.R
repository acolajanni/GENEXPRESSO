# Nom               : SimulateRNAseq
# Type              : Programme
# Objet             : Simule des données de type Microarrays
# Input             : None
# Output            : Simulmicroarrays.csv
# Auteur            : Juliette Casemajor
# R version         : 3.6
# Date de creation  : 04.03.2021
#______________________________________________________________________________

#Importation des packages
install.packages("madsim")
library(madsim)

#Charge un echantillon de vraies donnees microarrays
data(madsim_test)

#Definition des paramètres
mdata <- madsim_test$V1;
fparams <- data.frame(m1 = 6, m2 = 6, shape2 = 4, lb = 4, ub = 14,pde=0.02,sym=0.5)
#m1 et m2 indiquent respectivement le nombre d'echantillons controle et d'echantillons tests
#on garde les autres parametres par defaut ? oui pour l'instant, ça a l'air trè spécifique d'après la doc
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
sdn <- 0.4
rseed <- 50

#Generation des donnees
#n indique le nombre de genes, par defaut il est a 10000, on le met ici a 100
mydata1 <- madsim(mdata = NULL, n = 100, ratio = 0, fparams, dparams, sdn, rseed)
mydata1
#On recupere les donnees generees
Simulmicro <- as.data.frame(mydata1$xdata)




#Generation de noms de gènes aléatoires
listenom <- paste0(rep(LETTERS[1:12], each=9), rep(1:9, 10))
listenom = listenom[1:100]
row.names(Simulmicro) <- listenom

<<<<<<< HEAD
#Voir si les deux prochaines lignes sont utiles +> oui je trouve que c'est plus propre sur les noms de colonnes
nomcolonne <- c("sample1","sample2","sample3","sample4","sample5","sample6","sample7","sample8","sample9","sample10","sample11","sample12")
=======
#Voir si les deux prochaines lignes sont utiles > Oui
nomcolonne <- c("control1","control2","control3","control4","control5","control6","test1","test2","test3","test4","test5","test6")
>>>>>>> 9cbb6688648df4f4ffd034d6c744a6316ea3e6e8
colnames(Simulmicro)<-nomcolonne

############## Gros jeu de données (10.000) : 
fparams <- data.frame(m1 = 12, m2 = 12, shape2 = 4, lb = 4, ub = 14,pde=0.06,sym=0.5)
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
sdn <- 0.4
rseed <- 50
Simulmicro.big <- madsim(mdata = NULL, n = 10000, ratio = 0, fparams, dparams, sdn, rseed)

by(Simulmicro.big,Simulmicro.big$xid,nrow) 
# 291 downreg
# 296 upreg
Simulmicro.big = Simulmicro.big$xdata

listenom <- paste0(rep(LETTERS[1:26], each=400), rep(1:400, 26))
listenom = listenom[1:10000]
row.names(Simulmicro.big) <- listenom

listecol1 = paste0(rep("Control",each = 12),rep(1:12, each = 1) )
listecol2 = paste0(rep("Test",each = 12),rep(1:12, each = 1) )
listecol = c(listecol1,listecol2)
colnames(Simulmicro.big)<-listecol

############## Plus petit jeu de données (1.000) : 
fparams <- data.frame(m1 = 12, m2 = 12, shape2 = 4, lb = 4, ub = 14,pde=0.06,sym=0.5)
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
sdn <- 0.4
rseed <- 50
Simulmicro.smaller <- madsim(mdata = NULL, n = 1000, ratio = 0, fparams, dparams, sdn, rseed)

by(Simulmicro.big,Simulmicro.smaller$xid,nrow) 
# 21 downreg
# 35 upreg
Simulmicro.smaller = Simulmicro.smaller$xdata

listenom <- paste0(rep(LETTERS[1:26], each=40), rep(1:40, 26))
listenom = listenom[1:1000]
row.names(Simulmicro.smaller) <- listenom

listecol1 = paste0(rep("Control",each = 12),rep(1:12, each = 1) )
listecol2 = paste0(rep("Test",each = 12),rep(1:12, each = 1) )
listecol = c(listecol1,listecol2)
colnames(Simulmicro.smaller)<-listecol

#Exportation des donnees
write.csv(Simulmicro,"~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarraysname.csv", row.names = TRUE)
write.csv(Simulmicro.big,"~/GIT/CPRD/DATA/MICROARRAYS/SimulmicroarraysBIG.csv", row.names = TRUE)
write.csv(Simulmicro.smaller,"~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays1000.csv", row.names = TRUE)
