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
fparams <- data.frame(m1 = 12, m2 = 12, shape2 = 4, lb = 4, ub = 14,pde=0.02,sym=0.5)
#m1 et m2 indiquent respectivement le nombre d'echantillons controle et d'echantillons tests
#on garde les autres parametres par defaut ?
dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
sdn <- 0.4
rseed <- 50

#Generation des donnees
#n indique le nombre de genes, par defaut il est a 10000, on le met ici a 100
mydata1 <- madsim(mdata = NULL, n = 100, ratio = 0, fparams, dparams, sdn, rseed)
mydata1

#On recupere les donnees generees
micro <- as.data.frame(mydata1)
Simulmicro<-micro[,13:24]

#Voir si les deux prochaines lignes sont utiles
nomcolonne <- c("control1","control2","control3","control4","control5","control6","test1","test2","test3","test4","test5","test6")
colnames(Simulmicro)<-nomcolonne

#Exportation des donnees
write.csv(Simulmicro,"~/GIT/CPRD/DATA/MICROARRAYS/Simulmicroarrays.csv", row.names = TRUE)
