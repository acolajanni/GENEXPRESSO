#install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz ',repos=NULL)
#BiocManager::install("oligoClasses", dependancies = TRUE)
#BiocManager::install("affycoretools")

# à regarder : 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

# PC perso
celpath = "/home/antonin/Bureau/Cours M1 S7/Stage/GSE31684"
# Crémi 
#export RSTUDIO_WHICH_R=/net/cremi/acolajanni/Bureau/espaces/travail/R-3.6.3/bin/R
setwd("~/GIT/GENEXPRESSO")
celpath = "/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684"


getGEOSuppFiles("GSE31684", fetch_files = TRUE, baseDir = "./data")
untar(tarfile = "/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684/GSE31684_RAW.tar", exdir = "/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684")
gunzip(filename = "/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684/GSE31684_table_of_clinical_details.txt.gz",
       destname = "/net/cremi/acolajanni/Bureau/espaces/travail/GSE31684/GSE31684_table_of_clinical_details.txt")
            
       

f <- list.files(path = celpath, pattern = "CEL.gz", full.names = TRUE)
txt.dir = paste0(celpath,"/GSE31684_table_of_clinical_details.txt")
tab = read.delim(txt.dir,check.names=FALSE,as.is=TRUE, header = T)

# Compter les classes de cancers
by(tab,tab$PreOpClinStage,nrow)

# Ne garder que les T1 et T2
samples = subset(tab, tab$PreOpClinStage == 'T1' | tab$PreOpClinStage == 'T2') 
files = samples$GEO
#files
files = paste0(celpath,"/", files,".CEL.gz")

abatch <- ReadAffy(filenames = files)

# fait planter :
eset.mas5 = mas5(abatch)
# Fonctionne : 
eset = rma(abatch)

load("./Reunion.RData")

expr.val = exprs(eset)
head(row.names(expr.val))




############################################
##              Avec biomaRt              ##
############################################

library(biomaRt)
library(hgu133plus2.db)
# Liste des ensemble de sondes dispo dans hgu133plus2 + autres packages
ls("package:hgu133plus2.db")

mart = useMart(biomart = "ensembl"
                  #, dataset = 
                  )
# Toutes les annotations disponibles
listDatasets(mart)

# On récupère toutes les annotations pour Homo spaiens
ensembl = useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
# notre jeu de données utilise hgu133plus2
annotation = eset@annotation
# chercher les attributs qui nous interessent :
attributes = listAttributes(ensembl)
attributes = attributes$name
# On cherche le symbole des gènes : 

### att = str_detect(attributes$name, pattern = "symbol")
att = grepl(pattern = "symbol", attributes, fixed = TRUE)
test = attributes[grepl("symbol",attributes)]
attributes = listAttributes(ensembl)
# Vérification que c'est bien ce qu'on veut :
subset(attributes, attributes$name == test[1] | attributes$name == test[2])
# symbol human gene nomenclature comitee / uniprot symbol

affy_ensembl = c("affy_hg_u133_plus_2", "ensembl_gene_id", "uniprot_gn_symbol") # regarder uniprot symbol


#################
#################
#################

# si erreur getnodeset() : https://www.ensembl.org/info/website/archives/index.html le site est down

BM = getBM(attributes = affy_ensembl, 
           mart = ensembl, uniqueRows = TRUE, useCache = FALSE
           )
# On récupère la liste des sondes réellement présentes dans la bdd
BM2 = subset(BM, BM$affy_hg_u133_plus_2 != "")

#################
#################
#################

# Au final on garde le max() ==> à voir ce qu'on garde (mean, median, ...)




############################################
##              avec hgu.db               ##
############################################

library(hgu133plus2.db)
# Liste des ensemble de sondes dispo dans hgu133plus2 + autres packages
ls("package:hgu133plus2.db")
# Tous les attributs du package
columns(hgu133plus2.db)
# description des attributs (SYMBOL nous interesse)
help("SYMBOL")
# L'attribut qui nous interesse est "SYMBOL"
# On vérifie: OK
head(keys(hgu133plus2.db, keytype="SYMBOL"))

# On veut les ID des sondes : 
# ça correcpond à nos noms de lignes
k = keys(hgu133plus2.db,keytype="PROBEID")

Probes = select(hgu133plus2.db, keys=k, columns="SYMBOL", keytype="PROBEID")

#### Test de mapping : 

# On récupère la liste des sondes de notre échantillon
probeName = data.frame(PROBEID = row.names(expr.val))

library("hgu133plus2.db")
x <- hgu133plus2SYMBOL
Llength(x)
# Nombre / sondes de sondes mappés : 
count.mappedkeys(x)
nhit(x)
# nb mappé / non mappé
# non mappé : 10344
# mappé 1fois : 41922
# mappé plus d'une fois :
table(nhit(toggleProbes(hgu133plus2SYMBOL, "all")))
table(nhit(x))

test = toggleProbes(hgu133plus2SYMBOL,"multiple")
ls(test)

A = nhit(test)


?toggleProbes

probes.ALL=row.names(expr.val)


# ALIAS2PROBE mapping


# Sauvegarder l'avancée actuelle
# save.image(file = "./GSE31684.RData")


###############"
hgu133plus2.annot = read.csv("./data/HG-U133_Plus_2.na36.annot.csv", skip = 25, header=T, na.strings ="---")


ID.ALL = subset(hgu133plus2.annot, select = c("Probe.Set.ID","Gene.Symbol") )



#Remove control probes
# Trouver les sondes controles
#Extract the control probes
expr.val = as.data.frame(expr.val)
ControlProbes <- grep("AFFX",row.names(expr.val))


#remove control : 
Without.control.probe=expr.val[-ControlProbes,]
# Nb de sondes contrôles : 
nrow(expr.val[ControlProbes,])

# On récupère nos noms de sondes
probes.ALL=row.names(expr.val)
# Les symboles qui vont avec : 
symbol.ALL = unlist(mget(probes.ALL, hgu133plus2SYMBOL))
ID.ALL = unlist(mget(probes.ALL, hgu133plus2UNIPROT))
ID.ALL = as.data.frame(ID.ALL)

test = as.data.frame(unlist(mget(probes.ALL, hgu133plus2REFSEQ)))


# construction du tableau contenant les Symboles des gènes + ID des sondes + les données brutes
table.ALL=cbind(PROBES = probes.ALL,
                #ID.ALL,
                SYMBOL = symbol.ALL,
                expr.val)

# On récupère la liste des échantillons
samples = samples[!samples %in% c("PROBES","SYMBOL")]

#X = table.ALL
#X=X[order(X$values,decreasing=T),]
#Y=X[which(!duplicated(X$SYMBOL)),]


# En parallèle on créer un DF contenant les ID des sondes et des gènes
table.ID = data.frame(row.names = row.names(table.ALL),
                      PROBES = table.ALL$PROBES,
                      SYMBOL = table.ALL$SYMBOL)

#table.ID = data.frame(row.names = ID.ALL$Probe.Set.ID, SYMBOL = ID.ALL$Gene.Symbol, PROBES = ID.ALL$Probe.Set.ID)

# On enlève tous les ID dupliqués
ID.unique = table.ID[which(!duplicated(table.ID$SYMBOL)),]
ID.unique = ID.unique[order(ID.unique$SYMBOL),]

test = merge(ID.unique,expr.val)

#ID.unique = merge(ID.unique, Y)

for (i in 1:length(samples)){
  SampleName = samples[i]
  # On récupère un DF avec un échantillon uniquement
  SampleValue = subset(x = table.ALL, select = c(SampleName,"SYMBOL") )
  
  # On réordonne les valeurs du plus grand au plus petit
  X=SampleValue[order(SampleValue[SampleName],decreasing=T),]
  # On enlève les lignes où le Symbole est dupliqué
  # On enlève donc les lignes dont les valeurs sont inférieurs au max
  Y=X[which(!duplicated(SampleValue$SYMBOL)),]
  Y=Y[order(Y$SYMBOL),]
  
  ID.unique = cbind(ID.unique,Y[SampleName])
  #ID.unique = merge(ID.unique, Y, by = "SYMBOL")
  
}



# Avec la médiane :
# on reprend la formation du tableau
expr.val = as.data.frame(expr.val)
expr.val$PROBES = row.names(expr.val)
table.ID = data.frame(row.names = row.names(table.ALL),
                      PROBES = table.ALL$PROBES,
                      SYMBOL = table.ALL$SYMBOL)


test = merge(table.ID,expr.val,by = "PROBES")

nhit(test$PROBES)

test2 = group_by(test,test$PROBES)








#############################################"""
#############################################"""
y = hgu133plus2UNIPROT[isNA(hgu133plus2UNIPROT)]
Unmapped = Lkeys(y)   
table(Unmapped)
ID.Unmapped = unlist(mget(Unmapped, hgu133plus2UNIPROT))

ID.mapped = ID.ALL[!ID.ALL %in% ID.Unmapped]
#############################################"""
#############################################"""








#############################################"""
#############################################"""
## How many probes?
dim(hgu133plus2UNIPROT)
##Make a mapping with multiple probes exposed
multi <- toggleProbes(hgu133plus2UNIPROT, "all")
singleOnly <- toggleProbes(hgu133plus2UNIPROT, "single")
## How many probes?
dim(multi) #on a donc environ 7000 sondes qui mappent le même gene
dim(singleOnly)

count.mappedLkeys(hgu133plus2UNIPROT)
count.mappedLkeys(multi)

mapIds(multi,hgu133plus2UNIPROT)

hasSingleProbes(singleOnly)
#############################################"""
#############################################"""