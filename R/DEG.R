###########################
########    DEG   #########
###########################

library(edgeR)
library(DESeq)
library(DESeq2)
library(RankProd)
library(limma)
source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))


#' Creates a Model matrix.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param cond1 Name of the first experimental condition.
#' @param cond2 Name of the second experimental condition.
#' @param ncond1 Number of sample in the first experimental condition.
#' @param ncond2 Number of sample in the the second experimental condition.
#'
#' @return Model matrix.
#' @export
#'
#' @examples
#' 
#' ajouter un exemple
#' 
make_designMatrix <- function(dataset,cond1 = "A", cond2 = "B",ncond1=(ncol(dataset)/2),ncond2=(ncol(dataset)/2)){
  status.control = rep(cond1,ncond1)
  status.test = rep(cond2,ncond2)
  status = c(status.control,status.test)
  design = model.matrix(~0+status)
  colnames(design) <- c(cond1,cond2)
  return(design)
}  

#' Differentially expressed genes analysis through the limma package.
#'
#' Uses the lm.fit(),eBayes() and topTable() from the limma package to automate the steps of the analysis.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param design vector of 0 and 1 values. 0 for the first experimental condition, 1 for the second one.
#' @param contrast.matrix design matrix like one produced by the make_designMatrix() function.
#'
#' @return
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL)
#' @export
#'
#' @examples
DEG_limma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  # To fit the model to the data, it needs the model matrix
  fit <- lmFit(dataset,design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  # Classifying genes through their pvalue being differentially expressed
  res.diff <- topTable(fit2, coef="diff",genelist=row.names(dataset), number=Inf)
  res.diff_limma <- data.frame(FoldChange = (res.diff$logFC) ,PValue=(res.diff$adj.P.Val),SYMBOL=res.diff$ID)
  
  colnames(res.diff_limma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_limma)
}

#' Differentially expressed genes analysis through the limma package, with the GEOlimma method.
#'
#' Uses the lm.fit(),eBayes() and topTable() from the limma package to automate the steps of the analysis.
#'
#' @param dataset dataframe of expression values with samples in columns and genes in row.
#' @param design vector of 0 and 1 values. 0 for the first experimental condition, 1 for the second one.
#' @param contrast.matrix design matrix like one produced by the make_designMatrix() function.
#' 
#' @return
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL) 
#' @export
#' @examples
#' 
#' 
DEG_GEOlimma <- function(dataset,design, contrast.matrix = cm){
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  # To fit the model to the data, it needs the model matrix
  fit <- lmFit(dataset,design)
  fit2  <- contrasts.fit(fit, cm)
  # GEOlimma_probabilities.rda is a file that contains real probabilities of a gene being differentially expressed
  load("~/GIT/CPRD/GEOlimma/GEOlimma_probabilities.rda")
  # Those data are used to help the model to fit
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  # Classifying genes through their pvalue being differentially expressed
  de <- topTable(fit22, number = nrow(dataset))
  res.diff_geolimma <- data.frame(FoldChange = (de$logFC), PValue=(de$adj.P.Val),genes=row.names(de))
  
  colnames(res.diff_geolimma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_geolimma)
}

#' Takes the results of DEG_GEOlimma and DEG_limma to extract the pvalues of alternative hypothesis.
#'
#' @param res.diff_limma 
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue),
#' the last column corresponds to the gene names (SYMBOL).
#'
#' @return
#' A dataframe with 3 columns is returned. The two columns corresponds to the pvalue
#' of a gene being respectively, more expressed and less expressed in the second condition
#' @export
#'
#' @examples
#' #get the results of limma analysis
#' res.diff = DEG_limma(dataset,design, contrast.matrix = cm)
#' #compute the results of alternative hypothesis
#' res.diff.alternative = DEG_alternative(res.diff)
DEG_alternative <- function(res.diff_limma){
  # two dataframes are created from the original given in argument
  res.up = copy(res.diff_limma)
  res.down = copy(res.diff_limma)
  # If the pvalue is significant, and depending on the sign of the logFoldChange, we can set the probability to one for one hypothesis
  # The idea is that, if the fist condition has significantly higher expression values than the second one, the gene can't be more expressed in the second condition, and vice versa.
  for (row in 1:nrow(res.diff_limma)){
    if (res.diff_limma[['logFC']][row] > 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.down[['PValue']][row] = 1 
      
    }else if (res.diff_limma[['logFC']][row] < 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.up[['PValue']][row] = 1
    }
  }
  # This way, two columns are created, one for the probability of a gene being more expressed in the second condition (Pvalue_Up), and less expressed (Pvalue_Down)
  res = data.frame(PValue_Up = res.up$PValue, Pvalue_Down = res.down$PValue, SYMBOL = res.diff_limma$SYMBOL)
  return(res)
}


#' Multiple Wilcoxons tests for each row of a dataframe
#'
#' Testing, for one row at the time, if the first series of values are different, greater or less than the values of the second condition.
#'
#' @param data dataframe of gene expression levels: Gene names in rows, samples in columns.
#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return
#' Dataframe with three columns, each corresponding to a hypothesis : 
#' weak hypothesis, and two for the strong hypothesis "alternative = "less" ", and "greater"
#' For the second condition has values that are less or greater than values of the first one.
#' @export
#'
#' @examples
wilcoxDEG <- function(data, n1, n2){
  wilcox <- data.frame()
  # for each row in the dataset, three Wilconxon test are made, depending on the hypothesis.
  for (i in 1:nrow(data)){
    x = data[i,1:n1]
    y = data[i,(n1+1):(n1+n2)]
    test1 <- wilcox.test(as.numeric(x),as.numeric(y))
    test2 <- wilcox.test(as.numeric(x),as.numeric(y),alternative = "less")
    test3 <- wilcox.test(as.numeric(x),as.numeric(y), alternative = "greater")
    pval <- c(test1$p.value,test2$p.value,test3$p.value)
    wilcox <- rbind(wilcox,pval)
  }
  # The dataframe named "wilcox" contains 3 columns  that contains the pvalue for the three possible hypothesis
  colnames(wilcox)[1] <- "wilcox two-sided"
  colnames(wilcox)[2] <- "wilcox less"
  colnames(wilcox)[3] <- "wilcox greater"
  row.names(wilcox) = row.names(data)
  return(wilcox)
}

#' Compute the pvalues of genes being up and downregulated for Microarray datasets.
#'
#' Through different functions contained in several packages, this function computes pvalues of differentially expressed genes
#'
#' @param data Dataframe of gene expression levels with sample names in columns and genes in rows
#' 
#' @param tool Different methods used to compute pvalues of differentially expressed genes for microarrays
#' "Wilcox" uses the wilcoxDEG() function implemented in this very same pacakge 
#' "limma" and "GEOlimma uses respectively the functions DEG_limma() and DEG_GEOlimma() that comes from the limma pacakge
#' "RankProduct","RankProduct.log" perform a Rank Product analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 
#' "RankSum","RankSum.log" perform a Rank Sum analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 
#' 
#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return
#' A dataframe with genes in rows and pvalues of a gene being upregulated and downregulated in columns
#' @export
#'
#' @examples
#' 
#' 
tools.DEG.Microarrays <- function(data,tool,n1,n2){
  if (tool == "GEOlimma" || tool == "limma"){
    design = make_designMatrix(dataset = data)
  } else if (tool%in%c("RankProduct","RankProduct.log","RankSum","RankSum.log")){
    # we need a vector with 0 or 1 corresponding to the experimental condition (0 for the first, 1 for the other)
    design = rep(c(0,1),c(n1,n2)) 
  }
  
  DEG_Microarrays_tools.fnc <- switch(tool,
                                      GEOlimma = {
                                        # Compute the Pvalue with GEOlimma
                                        res.diff = DEG_GEOlimma(data,design)
                                        res.diff = DEG_alternative(res.diff)
                                        # Columns names "up" for upregaluted (for a gene more expressed in the second condition)
                                        # Down for downregulated
                                        colnames(res.diff) = c("GEOlimma Up","GEOlimma Down","Gene.ID")
                                      },
                                      
                                      limma = {
                                        res.diff = DEG_limma(data,design)
                                        res.diff = DEG_alternative(res.diff)
                                        colnames(res.diff) = c("limma Up","limma Down","Gene.ID")
                                        
                                      },
                                      
                                      Wilcox = {
                                        res.diff = wilcoxDEG(data, n1, n2)
                                        res.diff$Gene.ID = row.names(res.diff)
                                        # the last columns corresponds to the weak hypothesis.
                                        # We don't need it since we have both alternative hypothesis in the first two columns
                                        res.diff = res.diff[,-1]
                                      },
                                      
                                      RankProduct = {
                                        # Computing the DEG analysis
                                        res.diff = RankProducts(data, # dataset
                                                                design, #vector with 0 and 1 corresponding the the experimental condition of each sample
                                                                rand = 123, # equivalent to set seed, to have reproducible results
                                                                logged = FALSE, # wether or not the data are in log scale
                                                                na.rm = TRUE , # removing missing data
                                                                calculateProduct = TRUE) # TRUE for Rankproducts, FALSE for RankSum
                                        # the wanted results are in the "pval" objects inside res.diff
                                        res.diff = res.diff[["pval"]]
                                        
                                      },
                                      
                                      RankProduct.log = {
                                        res.diff = RankProducts(data, design, rand = 123, logged = TRUE , na.rm = TRUE , calculateProduct = TRUE)
                                        res.diff = res.diff[["pval"]]
                                        
                                      },
                                      RankSum = {
                                        res.diff = RankProducts(data, design, rand = 123,logged = FALSE ,na.rm = TRUE ,calculateProduct = FALSE)
                                        res.diff = res.diff[["pval"]]
                                        
                                      },
                                      
                                      RankSum.log = {
                                        res.diff = RankProducts(data, design, rand = 123,logged = TRUE ,na.rm = TRUE ,calculateProduct = FALSE)
                                        res.diff = res.diff[["pval"]]
                                        
                                      },
                                      stop("Enter something that switches me!")
  )
  if (!tool%in%c("limma","GEOlimma","Wilcox")){
    # for the RanProd analysis, some changes are needed : adding row names (genes) and changing columns names
    res.diff = as.data.frame(res.diff)
    res.diff$Gene.ID = row.names(data)
    colnames(res.diff) = c(paste(tool,"Up"), paste(tool,"Down"),"Gene.ID")
    res.diff
  }
  return(res.diff)
}

#' Compute the pvalues of genes being up and downregulated for Nanostring datasets.
#'
#' Through different functions contained in several packages, this function computes pvalues of differentially expressed genes
#'
#' @param raw.data rcc type file. List of 4 elements: 
#' Samples.IDs is a dataframe with four columns: sample name, status (WildType, Mutated), filename. 
#' rcc.df contains expression levels with samples in column and genes in rows. 
#' annots.df is a dataframe with 3 columns: geneclass (endogenous, housekeeping, negative, positive), gene name, and the accession number. 
#' FOV is a dataframe that contains Field of views information.
#' 
#' @param tool Method to use to compute the pvalues of differentially expressed genes.
#' "Wilcox" uses the wilcoxDEG() function implemented in this very same pacakge 
#' "limma" uses the functions DEG_limma() that comes from the limma pacakge
#' "RankProduct" and "RankSum" perform respectively a Rank Product and a Rank Sum analysiswith the RankProducts() function from the RankProd package
#' "desq2" uses the DESeq() function from the DESeq2 package. This last one is particular, because if it is chosen it will ignore the argument "tool_norm"
#'   
#' @param data Dataframe of gene expression levels with sample names in columns and genes in rows
#' 
#' @param tool_norm Normalization tool used previously in the tools.norm.Nanostring() function
#'
#' @return 
#' A dataframe with genes in rows and pvalues of a gene being upregulated and downregulated in columns
#' 
#' @export
#'
#' @examples
tools.DEG.Nanostring <- function(raw.data, tool, data, tool_norm) {
  
  rcc.samples <- raw.data$rcc.df
  annots.df <- raw.data$annots.df
  samples.IDs <- raw.data$samples.IDs
  FOV <- raw.data$FOV
  
  method = paste0(tool_norm,"_",tool)
  
  if (!tool%in%c("desq2")){
    group = table(samples.IDs$tp53.status)
    n1 = as.integer(group[1])
    n2 = as.integer(group[2])
  }
  if (tool%in%c("RankProduct", "RankSum")){
    design = rep(c(0,1),c(n1,n2))
  }
  
  
  tools.fnc <- switch(tool,
                      desq2={
                        # removing the Positive and negetive controls in the dataset
                        rcc.samples <- rcc.samples[!grepl("^Pos|^Neg",annots.df$CodeClass),]
                        samples.IDs <- samples.IDs[match(colnames(rcc.samples),samples.IDs$title),]
                        # Transform the dataset into a DESeqDataSet
                        dds <- DESeqDataSetFromMatrix(countData = rcc.samples,
                                                      colData = samples.IDs,
                                                      design= ~ tp53.status)
                        # Normalization
                        dds <- DESeq(dds)
                        # Retrieving pvalues depending on the alternative hypothesis
                        resG <- results(dds, altHypothesis="greater", format = "DataFrame")
                        resL <- results(dds, altHypothesis="less", format = "DataFrame")
                        res.diff = data.frame("DESeq2_Up"=(resG$padj),"DESeq2_Down"=(resL$padj), "SYMBOL" = row.names(resG))
                      },
                      
                      limma = {
                        # Constructing the model matrix
                        design <- model.matrix(~0+samples.IDs$tp53.status)
                        colnames(design) <- c("A","B")
                        # Computing Pvalues
                        res.diff = DEG_limma(data,design)
                        res.diff = DEG_alternative(res.diff)
                        # creating names of the used method (normalization tool + Statistical method)
                        # Up and Down for upregulated and downregulated
                        method_Up = paste0(tool_norm,"_",tool,"_Up")
                        method_Down = paste0(tool_norm,"_",tool,"_Down")
                        
                        res.diff <- data.frame(A=(res.diff$PValue_Up), B=(res.diff$Pvalue_Down) ,SYMBOL=res.diff$SYMBOL)
                        colnames(res.diff) = c(method_Up,method_Down, "SYMBOL")
                        
                      },
                      
                      Wilcox = {
                        # Computing Pvalues
                        res.diff = wilcoxDEG(data, n1, n2)
                        # Discard the columns of weak hypothesis since we want the upregulated and downregulated genes
                        res.diff = res.diff[(-1)]
                        
                        res.diff$SYMBOL = row.names(res.diff)
                        method_less = paste0(tool_norm,"_",tool,"_less")
                        method_greater = paste0(tool_norm,"_",tool,"_greater")
                        colnames(res.diff) = c(method_less, method_greater, "SYMBOL")
                      },
                      
                      RankProduct = {
                        res.diff = RankProducts(data, design, rand = 123, logged = TRUE , na.rm = TRUE , calculateProduct = TRUE)
                        res.diff = res.diff[["pval"]]
                        
                      },
                      
                      RankSum = {
                        res.diff = RankProducts(data, design, rand = 123,logged = TRUE ,na.rm = TRUE ,calculateProduct = FALSE)
                        res.diff = res.diff[["pval"]]
                      },
                      
                      stop("Enter something that switches me!") 
                      
  )
  if (tool%in%c("RankProduct","RankSum")){
    # Rename columns with used methods
    res.diff = as.data.frame(res.diff)
    res.diff$SYMBOL = row.names(data)
    method_Up = paste0(method,"_Up")
    method_Down = paste0(method,"_Down")
    colnames(res.diff) = c(method_Up,method_Down,"SYMBOL")
  }
  return(res.diff)
}

#' Compute the pvalues of genes being differentially expressed for RNA-seq type data.
#'
#' Through different functions contained in several packages, this function computes pvalues of differentially expressed genes
#'
#' @param data dataframe of gene expression levels, with genes in rows and samples in columns. 
#' @param tool Method to Normalize datasets and statistically analysing them
#' "edgeR_RLE","edgeR_upperquartile","edgeR_TMMwsp", and "edgeR_TMM" are methods of normalization implemented in the edgeR package with the function calcNormFactors.DGEList()
#' If "tool" is one of these parameters, it will compute two analysis: the exact test from the exactTest() function 
#' and a Quasi-likelihood test with the glmQLFTest() function still in the edgeR package.
#' "deseq2.Wald" and "deseq2.LRT' uses the same normalization methods contained in the DESeq2 package with de DESeq() function.
#' The first method uses a Wald test and the second a Likelihood ratio test to determine Pvalues of differentially expressed genes.
#'
#' @return
#' Dataframe with genes in row, and methods used in columns. In contains the differentially expressed p-values for each gene. 
#' @export
#'
#' @examples
#' 
tools.DEG.RNAseq.inspect <- function(data,tool){
  data = as.matrix(data)
  # One column dataframe with genes inside
  res.diff = data.frame("genes" = row.names(data))
  
  tools_norm_RNAseq.fnc <- switch(tool,
                                  
                                  edgeR_TMM = {
                                    # Dataset into DEGlist (edgeR class)
                                    DGE = DGEList(count = data, 
                                                  group = rep(1:2,each=ncol(data)/2)) # functions for model with n1 = n2 
                                    # Normalization 
                                    res.norm <- calcNormFactors.DGEList(DGE, method = "TMM")
                                    tool_name = "TMM"
                                  },
                                  
                                  edgeR_RLE = {
                                    DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
                                    res.norm <- calcNormFactors.DGEList(DGE, method = "RLE")
                                    tool_name = "RLE"
                                  },
                                  
                                  edgeR_TMMwsp = {
                                    DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
                                    res.norm <- calcNormFactors.DGEList(DGE, method = "TMMwsp")
                                    tool_name = "TMMwsp"
                                  },
                                  
                                  edgeR_upperquartile = {
                                    DGE = DGEList(count = data, group = rep(1:2,each=ncol(data)/2))
                                    res.norm <- calcNormFactors.DGEList(DGE, method = "upperquartile")
                                    tool_name = "Upperquartile"
                                  },
                                  
                                  deseq2.Wald = {
                                    # data into DESeqDataSet (DESeq2 class) :
                                    condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
                                    type = factor(rep("single-read",ncol(data)))
                                    colData = data.frame(condition,type,row.names=colnames(data))
                                    dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)
                                    
                                    # Retrieve DEG analysis (test Wald)
                                    DEG <- DESeq(dds, test = "Wald")
                                    DEG <- results(DEG)

                                    # adding a column with the pvalues
                                    res.diff <- data.frame(deseq2_Wald=(DEG$padj),genes=row.names(DEG))
                                    
                                  },
                                  
                                  deseq2.LRT = {
                                    condition <-factor( rep(c("control","case"),each = (ncol(data)/2)) )
                                    type = factor(rep("single-read",ncol(data)))
                                    colData = data.frame(condition,type,row.names=colnames(data))
                                    dds<-DESeqDataSetFromMatrix(data,colData,design=~condition)
                                    
                                    DEG <- DESeq(dds, test = "LRT",reduced = ~1)
                                    DEG <- results(DEG)
                                    
                                    res.diff <- data.frame(deseq2_LRT=(DEG$padj),genes=row.names(DEG))
                                  },
                                  
                                  deseq = {
                                    condition <- rep(c("control","case"),each = (ncol(data)/2) )
                                    type = factor(rep("single-read",ncol(data)))
                                    design = data.frame(condition,type,row.names=colnames(data))
                                    singleSamples = design$type == 'single-read'
                                    countTable = data[,singleSamples]
                                    conds = design$condition[singleSamples]
                                    
                                    cds <- newCountDataSet(data, conds)
                                    cds = estimateSizeFactors(cds)
                                    cds = estimateDispersions(cds)
                                    DEG = nbinomTest(cds, condA = "control", condB = "case")
                                    res.diff <- data.frame(deseq = (DEG$padj),genes= (DEG$id))
                                    
                                  },
                                  stop("Enter something that switches me!") 
                                  
  )
  # DEG analysis for not yet treated data (only edgeR methods)
  if (!tool%in%c("deseq2.Wald","deseq2.LRT", "deseq")){
    # Exact Test (edgeR)
    colname1 = paste(tool_name,"ExactTest")
    # Dispersions 
    Disp = estimateCommonDisp(res.norm)
    Disp = estimateTagwiseDisp(Disp)
    # Retrieving pvalues
    pvalue = exactTest(Disp)$table[3]
    
    # 2 columns dataframe : Pvalues + genes
    res.diff1 = data.frame(genes = row.names(pvalue),pvalue = pvalue$PValue)
    
    # GLM (edgeR)
    colname2 = paste(tool_name,"GLM")
    # design matrix : 
    design = model.matrix(~0+group, data = res.norm$samples)
    colnames(design) <- c("Control","Test")
    # Dispersions
    y = estimateDisp(res.norm,design)
    # fitting model
    fit <- glmQLFit(y, design)
    BvsA <- makeContrasts(Control-Test, levels=design)
    # DEG analysis
    qlf <- glmQLFTest(fit, contrast=BvsA)
    # genes + pvalues are retrieved
    pvalue = qlf[["table"]][4]
    res.diff2 = data.frame(genes = row.names(pvalue),pvalue = pvalue$PValue)
    res.diff = merge(res.diff1,res.diff2,by = "genes",all=T)
    
    # On renaming the two last columns
    names(res.diff)[-1][-2] = colname1
    names(res.diff)[-1][-1] = colname2
  }
  return(res.diff)
}