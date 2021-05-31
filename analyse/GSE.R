install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz ',repos=NULL)
BiocManager::install("oligoClasses", dependancies = TRUE)
BiocManager::install("affycoretools")

# à regarder : 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

celpath = "/home/antonin/Bureau/Cours M1 S7/Stage/GSE31684"
f <- list.files(path = celpath, pattern = "CEL.gz", full.names = TRUE)
txt.dir = paste0(celpath,"/GSE31684_table_of_clinical_details.txt")
tab = read.delim(txt.dir,check.names=FALSE,as.is=TRUE)

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
#att = grepl(pattern = "symbol", attributes, fixed = TRUE)
#test = attributes[grepl("symbol",attributes)]
#attributes = listAttributes(ensembl)
# Vérification que c'est bien ce qu'on veut :
#subset(attributes, attributes$name == test[1] | attributes$name == test[2])
# symbol human gene nomenclature comitee / uniprot symbol

affy_ensembl = c("affy_hg_u133_plus_2", "ensembl_gene_id", "hgnc_symbol") # regarder uniprot symbol
# si erreur getnodeset() : https://www.ensembl.org/info/website/archives/index.html le site est down
BM = getBM(attributes = affy_ensembl, mart = ensembl, uniqueRows = TRUE, useCache = FALSE)
# On récupère la liste des sondes réellement présentes dans la bdd
BM2 = subset(BM, BM$affy_hg_u133_plus_2 != "")

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
head(k)

Probes = select(hgu133plus2.db, keys=k, columns="SYMBOL", keytype="PROBEID")

#### Test de mapping : 

# On récupère la liste des sondes de notre échantillon
probeName = data.frame(PROBEID = row.names(expr.val))

library(dplyr)

test = rbind(probeName,Probes)





