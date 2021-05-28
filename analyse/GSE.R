install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz ',repos=NULL)
BiocManager::install("oligoClasses", dependancies = TRUE)
BiocManager::install("affycoretools")


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


expr.val = exprs(eset)
head(row.names(expr.val))

library(biomaRt)
mart = useMart(biomart = "ensembl"
                  #, dataset = 
                  )
# Toutes les annotations disponibles
Dataset = listDatasets(mart)

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

affy_ensembl = c("affy_hg_u133_plus_2", "ensembl_gene_id", "hgnc_symbol")
BM = getBM(attributes = affy_ensembl, mart = ensembl, uniqueRows = TRUE, useCache = FALSE)
# On récupère la liste des sondes réellement présentes dans la bdd
BM2 = subset(BM, BM$affy_hg_u133_plus_2 != "")


