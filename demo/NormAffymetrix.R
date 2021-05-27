#BiocManager::install("gcrma")
library(affy)
library(gcrma)
library(GEOquery)
#library(limma)
#library(affydata)
library(affyio)
library(stringr)


getGEOSuppFiles("GSE781", fetch_files = TRUE, baseDir = "./data")
untar(tarfile = "./data/GSE781/GSE781_RAW.tar", exdir = "./data/GSE781")
###############################
celpath = file.path("./data/GSE781/")
f <- list.files(path = celpath, pattern = "CEL.gz", full.names = TRUE)

# Pas indispensable
tab <- read.delim("./data/Affymetrix/filelist.txt",check.names=FALSE,as.is=TRUE)
rownames(tab) <- tab$Name
##

## Dépend de l'échantillon
celf <- lapply(f, function(x) ReadAffy(filenames = x))
table(sapply(celf, annotation))
ff <- split(f, sapply(f, function(x) read.celfile.header(x)$cdfName))

# En garder qu'un
abatch1 <- ReadAffy(filenames = ff$`HG-U133A`)
abatch2 <- ReadAffy(filenames = ff$`HG-U133B`)

abatch1$sample
pData(abatch1)
# Accéder à d'autres données : probes
feat = abatch1@featureData
# Tableau de design (ici vide)
exp = abatch1@experimentData
# Donner des renseignements sur les phénotypes (mutants vs sains) : 

#1 : Récupérer les données phénotypiques : 
gse <- getGEO('GSE781',GSEMatrix=TRUE)
Pheno = (pData(phenoData(gse[[1]])))
# Créer la matrice de design pour les analyses DEG
ph = abatch1@phenoData
ph@data[,1] = Pheno$title
ph@data
design = as.vector(ph@data[,1])
design = startsWith(design, "N")
design = ifelse(design == TRUE, 0, 1)
design
ph@data[,1] = design
design.matrix = as.data.frame(ph@data)
design.matrix

# MAS5 par défaut
eset.mas5 = mas5(abatch1)

# Pour plus de paramétrage : (on peut changer chaque méthode)
eset <- expresso(abatch1, bgcorrect.method="rma",
                 normalize.method="constant",
                 pmcorrect.method="pmonly"
                 ,summary.method="avgdiff"
                 )


abatchBG = bg.correct(abatch1, method = "rma")
abatchBG_constant = normalize.AffyBatch.constant(abatch = abatchBG)
abatchBG_constant_pmonly = pmcorrect.pmonly(abatchBG_constant)
abatchBG_constant_pmonly_avgdiff = express.summary.stat(abatchBG_constant, 
                                                        pmcorrect = "pmonly",
                                                        summary = "avgdiff")

# RMA par défaut
eset.rma = rma(abatch1)


eset.gcrma = gcrma(abatch1)



# Extraction pour analyse DEG
exprSet = exprs(eset)
exprSet.mas5 = exprs(eset.mas5)
exprSet.rma = exprs(eset.rma) #rma log les expressions
exprSet.gcrma = exprs(eset.gcrma) #rma log les expressions

# Si on veut log
exprSet = log(exprSet,2)


# Options personallisable allant avec la fonction expresso()
norm.methods = normalize.AffyBatch.methods()
#norm.methods = norm.methods[-c("methods","vsn")]
norm.methods = norm.methods [norm.methods != "methods" & norm.methods != "vsn" ]



bgcorrect.methods = bgcorrect.methods()
bgcorrect.methods = bgcorrect.methods[bgcorrect.methods == "rma" | bgcorrect.methods == "mas" ]
bgcorrect.methods

pmcorrect.methods = pmcorrect.methods()
pmcorrect.methods = pmcorrect.methods[pmcorrect.methods != "methods"]


sum.stat = express.summary.stat.methods()




# TRouver le nom des gènes (et non des sondes) : Avoir les packages correspondant aux bonnes sondes