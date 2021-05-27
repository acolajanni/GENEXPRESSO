library(GEOquery)


# GSM = données pour un individu
gsm <- getGEO("GSM11805")
# GDS = Tableau d'expression gènes pour plusieurs individus ? tableau d'expression = 
gds <- getGEO("GDS507")
gds@dataTable@columns
gds@dataTable@table

# Plus lourd que les autres, contient toutes les données
gse <- getGEO("GSE781",GSEMatrix=TRUE)
show(gse)

# récupérer un expression tableau d'expression à partir de gds
# Possible de faire une transfo log mais peut être pas intéressant ? (depénd de la normalisation)
eset <- GDS2eSet(gds
                 #,do.log2=TRUE
                 #
                 )

# DF d'expression
Expression = as.data.frame(eset@assayData[["exprs"]])
# Design : conditions normal vs cancer
Design = eset@phenoData@data

# Faire correspondre les "GSM" avec les individus pour avoir la condition exp
# Pareil pour les gènes

Genes = as.data.frame(eset@featureData@data)
Genes = data.frame(Genes$ID, Genes$`Gene symbol`)

