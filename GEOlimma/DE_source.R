

limma_DE <- function(X, Y, dataset = NULL){
    
    library(limma)
    #library(bigmemory)
    
    exprs <- t(X)
    design <- data.frame("Case" = as.numeric(as.character(Y)), "Control" = 1-as.numeric(as.character(Y)) )
    fit <- lmFit(exprs, design)
    cont.matrix <- makeContrasts(compare="Case-Control",levels=design)
    fit2  <- contrasts.fit(fit, cont.matrix)
    #print("Test if success")
    ###limma B -statistics
    fit21  <- eBayes(fit2)
    #print(fit21)
    #print(head(fit21$lods))
    #print(head(rownames(fit21$lods)))
    de <- topTable(fit21, number = nrow(exprs))
    #de <- data.frame(de, stringsAsFactors = F)
    #newde <- de[rownames(exprs), ]
    #de <- de[order(de$B, decreasing = T), ]
    de <- de[order(de$adj.P.Val, decreasing = F), ]
    #print(newde)
    #print(head(sort(newde$B, decreasing = T)))
    
    return(de)
    
}


GEOlimma_DE <- function(X, Y, dataset = NULL){
    
    library(limma)
    #library(bigmemory)
    #library(pheatmap)
    
    source("~/GIT/CPRD/GEOlimma/ebayesGEO.R")
    load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
    
    #count genes in prop and in X
    #if(!grepl("X", rownames(prop)[1])){
    #rownames(prop) = paste("X", rownames(prop), sep = "")
    #}
    
    print("DE frequency found in genes: ")
    print(length(intersect(colnames(X), rownames(prop))))
    
    exprs <- t(X)
    design <- data.frame("Case" = as.numeric(as.character(Y)), "Control" = 1-as.numeric(as.character(Y)) )
    
    #print(design)
    fit <- lmFit(exprs, design)
    cont.matrix <- makeContrasts(compare="Case-Control",levels=design)
    fit2  <- contrasts.fit(fit, cont.matrix)
    ###limma B -statistics
    
    #fit21  <- eBayes(fit2)
    #de <- topTable(fit21, number = nrow(exprs))
    #de <- de[order(de$B, decreasing = T), ]
    #nonchangeB = de
    ###geolimma test
    fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
    #head(fit22$p.value)
    #print(fit22)
    #print(head(fit22$lods))
    #print(head(rownames(fit22$lods)))
    #fit22 <- unclass(fit22)
    
    #toptable(fit22 = fit[c("coefficients", "stdev.unscaled")],
    #         coef = coef, number = number, genelist = genelist, A = fit$Amean,
    #    eb = fit[ebcols], adjust.method = adjust.method, sort.by = sort.by,
    #     resort.by = resort.by, p.value = p.value, lfc = lfc,
    #     confint = confint)
    
    de <- topTable(fit22, number = nrow(exprs))
    #de <- de[order(de$B, decreasing = T), ]
    #de <- data.frame(de, stringsAsFactors = F)
    #newde <- de[rownames(exprs), ]
    #print(newde)
    de <- de[order(de$adj.P.Val, decreasing = F), ]
    
    return(de)
    
}





