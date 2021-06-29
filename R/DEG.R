###########################
########    DEG   #########
###########################

library(edgeR)
library(DESeq)
library(DESeq2)
library(RankProd)
library(limma)
library(data.table)
#source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
#source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))


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
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Creating the design matrix
#' design = make_designMatrix(dataset = Data, 
#'                             cond1 = "control", 
#'                             cond2 = "Case", 
#'                             ncond1 = 10, 
#'                             ncond2 = 10)
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
#' 
#' @import "limma"
#' @export
#'
#' @examples
#' # Import the dataset
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes  
#' 
#' #Build the design matrix
#' #design = make_designMatrix(dataset = Data,ncond1 = 10, ncond2 = 10)
#' #res.DEG = DEG_limma(Data,design)
DEG_limma <- function(dataset,design, contrast.matrix){
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
#' @param contrast.matrix design matrix like the one produced by the make_designMatrix() function.
#' 
#' @return
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue)
#' the last column corresponds to the gene names given in the dataset (SYMBOL) 
#' 
#' @docType data
#' @usage data(GEOlimma_probabilities)
#' @import "limma"
#' @export
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes  
#' 
#' #Build the design matrix
#' #design = make_designMatrix(dataset = Data,ncond1 = 10, ncond2 = 10)
#' #get the results of limma analysis
#' #res.diff = DEG_GEOlimma(dataset = Data,design)
DEG_GEOlimma <- function(dataset,design, contrast.matrix){
  source(file.path("~/GIT/CPRD/GEOlimma/","DE_source.R"))
  source(file.path("~/GIT/CPRD/GEOlimma/","ebayesGEO.R"))
  
  if(missing(contrast.matrix)){
    cm <- makeContrasts(diff = B-A, levels=design)
  }
  # To fit the model to the data, it needs the model matrix
  fit <- lmFit(dataset,design)
  fit2  <- contrasts.fit(fit, cm)
  # GEOlimma_probabilities.rda is a file that contains real probabilities of a gene being differentially expressed
  #data("GEOlimma_probabilities")
  prop = GENEXPRESSO::prop
  # These data are used to fit the model
  fit22  <- eBayesGEO(fit2, proportion_vector=prop[, 1, drop=F])
  # Classifying genes through their pvalue being differentially expressed
  de <- topTable(fit22, number = nrow(dataset))
  res.diff_geolimma <- data.frame(FoldChange = (de$logFC), PValue=(de$adj.P.Val),genes=row.names(de))
  
  colnames(res.diff_geolimma) <- c("logFC","PValue","SYMBOL")
  
  return(res.diff_geolimma)
}

#' Takes the results of DEG_GEOlimma and DEG_limma and extracts the pvalues of alternative hypothesis.
#'
#' @param res.diff_limma 
#' A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
#' The Pvalue associated with this Pvalue (PValue),
#' the last column corresponds to the gene names (SYMBOL).
#'
#' @return
#' A dataframe with 3 columns is returned. The two columns corresponds to the pvalue
#' of a gene being more expressed and less expressed, in the second condition
#' 
#' @import "data.table"
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' #Build the design matrix
#' #design = make_designMatrix(dataset = Data,ncond1 = 10, ncond2 = 10)
#' #get the results of limma analysis
#' #res.diff = DEG_limma(dataset = Data,design)
#' 
#' #compute the results of alternative hypothesis
#' #res.diff.alternative = DEG_alternative(res.diff)
DEG_alternative <- function(res.diff_limma){
  # two dataframes are created from the original given in argument
  res.up = copy(res.diff_limma)
  res.down = copy(res.diff_limma)
  # If the pvalue is significant, and depending on the sign of the logFoldChange, we can set the probability to one for one hypothesis
  # If the fist condition has significantly higher expression values than the second one, the gene can't be more expressed in the second condition, and vice versa.
  for (row in 1:nrow(res.diff_limma)){
    if (res.diff_limma[['logFC']][row] > 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.down[['PValue']][row] = 1 
      
    }else if (res.diff_limma[['logFC']][row] < 0 && res.diff_limma[['PValue']][row] <= 0.05){
      res.up[['PValue']][row] = 1
    }
  }
  # This way, two columns are created, one for the probability of a gene being more expressed (Pvalue_Up) and less expressed (Pvalue_Down) in the second condition
  res = data.frame(PValue_Up = res.up$PValue, Pvalue_Down = res.down$PValue, SYMBOL = res.diff_limma$SYMBOL)
  return(res)
}


#' Wilcoxon test for each row of a dataframe
#'
#' Testing, for one row at the time, if the first series of values are different, greater or less than the values of the second condition.
#'
#' @param data dataframe of gene expression levels: Gene names in rows, samples in columns.
#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return
#' dataframe with pvalues in one column and rownames of data as rownames.
#' @export
#' 
#' @import "dplyr"
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute Pvalues
#' res.DEG = wilcoxDEG(Data,10,10)
wilcoxDEG <- function(data, n1, n2){
  # for each row in the dataset, one wilcoxon test is made depending on the hypothesis (strong or weak) 
  pvals.less=apply(data,1,function(x) {
    wilcox.test(x[1:n1],x[n1+1:length(x)], 
                alternative="less")$`p.value`
  })
  
  pvals.greater=apply(data,1,function(x) {
    wilcox.test(x[1:n1],x[n1+1:length(x)], 
                alternative="greater")$`p.value`
  })
  
  # put it in a dataframe with gene SYMBOL as rownames
  wilcox=data.frame(Gene.ID = row.names(data) , 
                    less = pvals.less, 
                    greater = pvals.greater)
  
  # correcting the p-values with False Direcovery Rate
  wilcox$less = p.adjust(wilcox$less, method="BH")
  wilcox$greater = p.adjust(wilcox$greater, method="BH")
  row.names(wilcox) = NULL
  
  return(wilcox)
}

#' Compute the pvalues of genes being up- and down-regulated for Microarray datasets.
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
#' 
#' @import "RankProd" "limma" "dplyr"
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute Pvalues for the RankProduct method
#' res.DEG = tools.DEG.Microarrays(Data,"RankProduct",10,10)
tools.DEG.Microarrays <- function(data,tool,n1,n2, design){
  if (missing(design)){
    if (tool == "GEOlimma" || tool == "limma"){
      design = make_designMatrix(dataset = data)
    } else if (tool%in%c("RankProduct","RankProduct.log","RankSum","RankSum.log")){
      # we need a vector of 0 and 1 corresponding to the experimental condition (0 for the first, 1 for the other)
      design = rep(c(0,1),c(n1,n2)) 
    }
    
  } else if (tool%in%c("RankProduct","RankProduct.log","RankSum","RankSum.log", "Wilcox")){
    # we need a vector of 0 and 1 corresponding to the experimental condition (0 for the first, 1 for the other)
    design = as.data.frame(design)
    design = as.vector(ifelse(design[1] == 0, 
                    yes = 0,
                    no = 1))
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
                                        
                                        message("Reordering samples")
                                        
                                        names(design) = colnames(data)
                                        g1 = names(design)[design == 0]
                                        g2 = names(design)[design == 1]
                                        data = relocate(data, g1, g2)
                                        
                                        message("computing pvalues")
                                        res.diff = wilcoxDEG(data, length(g1), length(g2))
                                        colnames(res.diff) = c("Gene.ID","Wilcox Up","Wilcox Down")
                                        return(res.diff)


                                      },
                                      
                                      RankProduct = {
                                        # Computing the DEG analysis
                                        res.diff = RankProducts(data, # dataset
                                                                design, #vector of 0 and 1 corresponding the the experimental condition of each sample
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
  if (!tool%in%c("limma","GEOlimma")){
    # for the RanProd analysis, some changes are needed : adding row names (genes) and changing columns names
    res.diff = as.data.frame(res.diff)
    res.diff$Gene.ID = row.names(data)
    colnames(res.diff) = c(paste(tool,"Up"), paste(tool,"Down"),"Gene.ID")
    
    if(tool %in%c("RankProduct.log", "RankProduct", "RankSum", "RankSum.log")){
      # Compute adjusted pvalues (first for Up then for Down) for RankProducts pvalues
      res.diff[[ colnames(res.diff)[1] ]]=p.adjust(res.diff[[ colnames(res.diff)[1] ]],method = "BH")
      res.diff[[ colnames(res.diff)[2] ]]=p.adjust(res.diff[[ colnames(res.diff)[2] ]],method = "BH")
    }
  }
  return(res.diff)
}

#' Compute the pvalues of genes being up and downregulated for Nanostring datasets.
#'
#' Through different functions contained in several packages, this function computes pvalues of diferentially expressed genes
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
#' @import "DESeq2" "RankProd" "limma"
#' @export
#'
#' @examples 
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Normalizing data using one method
#' # Norm.data = tools.norm.Nanostring(raw.data = Data, tool = "nappa.NS")
#' # Analyze normalized data with one DEG method
#' #res.DEG = tools.DEG.Nanostring(raw.data = Data, 
#' #                                 data = Norm.data, 
#' #                                 tool_norm = "nappa.NS", 
#' #                                 tool = "RankProduct")
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


#' Compute pvalues for a gene to be differentially expressed.
#'
#' This function calls various function to compute the pvalues for each rows of the dataframe (gene).
#' It has been conceived to compute pvalues from raw count matrix with pre calculated norm factors
#'
#' @param count.matrix.raw 
#' Dataframe of count with samples in columns and genes SYMBOL in rows.
#' 
#' @param nf 
#' vector of normalisation / size factors to analyse with "ExactTest", "GLM", "nbinom", "nbinom.Wad", "nbinom.LRT"
#' Matrix of normalized count or "Elist" for the parameter "limma"
#' 
#' @param tool.norm 
#' "TMM","TMMwsp", "RLE", "Upperquartile", "voom", "vst", "vst2"
#' Note that those character strings are just a way to rename columns for the returned pvalues 
#' 
#' 
#' @param tool.DEG 
#' Character string among : "ExactTest", "GLM", "nbinom", "nbinom.Wad", "nbinom.LRT"
#' "ExactTest calls the \link{edgeR}{exactTest} function.
#' "GLM" uses a linear model with the \link{edgeR}{glmQLFit} and \link{edgeR}{glmQLFTest} functions.
#' "nbinom" is the DESeq equivalent with the \link{DESeq}{nbinomTest} function.
#' "nbinom.Wald" and "nbinom.LRT" calls the same function with different parameters : \link{DESeq2}{DESeq}.
#' 
#' @param design 
#' Vector of 1 and 2 of the same length of colnames(count.matrix).
#' 1 for the first group and 2 for the second.
#'
#' @import "DEFormats" "edgeR" "DESeq2" "DESeq" "limma"
#' 
#' @return
#' @export
#'
#' @examples
#' # load a count matrix (example with a random dataset)
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute design vector
#' design = c(rep(1,10), rep(2,10)) # 10 from group 1, 10 from group 2
#' 
#' # Compute normalization factors
#' nf = runif(min = 0.75, max = 1.25, n = 20)
#' 
#' DEG = tools.DEG.RNAseq(Data, nf, "Random.norm", "ExactTest", design )
#' 
tools.DEG.RNAseq <- function(count.matrix.raw, nf, tool.norm, tool.DEG, design){
  storage.mode(count.matrix.raw) = "integer"
  # limma only apply for voom normalization
  if (class(nf) == "EList"){
    # Elist are class object for the voom normalization
    tool.norm = "voom"
    tool.DEG = "limma"
    data = nf 
    method = "voom.limma"
    
  }
  
  else if(class(nf) == "matrix"){
    # Matrix are the output for vst normalization
    tool.norm = "vst"
    tool.DEG = "limma"
    data = nf
    method = "vst.limma"
  }
  
  else{
    method = paste0(tool.norm,".",tool.DEG)
  }
  
  # Specific class is needed for edgeR analysis
  if (tool.DEG %in% c("ExactTest","GLM")){
    edgeR.dgelist = DGEList(counts = count.matrix.raw, group = factor(design))
    
    # Estimating dispersion, retrieve and apply normalization factors
    edgeR.dgelist[["samples"]][["norm.factors"]] = nf
    edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
    edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist,
                                        trend = "movingave")
    
  }
  # Same for DESeq2
  else if (tool.DEG%in%c("nbinom.Wald","nbinom.LRT")){
    design = data.frame(design,row.names=colnames(count.matrix.raw))
    design$design = as.factor(design$design)
    dds<-DESeqDataSetFromMatrix(count.matrix.raw,
                                colData = design,
                                design= ~design)
    sizeFactors(dds) = nf
    
    
    dds = estimateDispersions(dds, 
                              fitType = "local")
    
  }
  # Same again for DESeq
  else if (tool.DEG == "nbinom"){
    DESeq.cds = newCountDataSet(countData = count.matrix.raw,
                                conditions = factor(design))
    # Initializing size factors
    sizeFactors(DESeq.cds) = nf
  }
  
  tools_norm_RNAseq.fnc <- switch(tool.DEG,
                                  # EdgeR exact test
                                  ExactTest = {
                                    
                                    DEG = exactTest(edgeR.dgelist)
                                    DEG.pval = DEG$table$PValue
                                    
                                    
                                  },
                                  
                                  # EdgeR General linear model
                                  GLM = {
                                    design = model.matrix(~0+group, data = edgeR.dgelist$samples)
                                    colnames(design) <- levels(edgeR.dgelist$samples$group)
                                    
                                    # Fitting the GLM model
                                    fit <- glmQLFit(edgeR.dgelist, design)
                                    # Testing outliers
                                    qlf <- glmQLFTest(fit, contrast=c(-1,1))
                                    # Retrieving pvalues
                                    DEG.pval = qlf[["table"]]$PValue
                                    
                                  },
                                  
                                  # DESeq analysis
                                  nbinom = {
                                    # Dispersions
                                    DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum",     
                                                                    method = "pooled", 
                                                                    fitType = "local")
                                    
                                    # Searching for DE genes
                                    DESeq.test = nbinomTest(DESeq.cds, "1", "2")
                                    DEG.pval = DESeq.test$pval
                                    
                                  },
                                  
                                  # DESeq2 Wald Test
                                  nbinom.Wald = {
                                    DEG <- DESeq(dds, test = "Wald")
                                    
                                  },
                                  
                                  # DESeq2 Likelihood ratio test for GLMs
                                  nbinom.LRT = {
                                    DEG <- DESeq(dds, test = "LRT",reduced = ~1)
                                    
                                  },
                                  
                                  
                                  limma = {
                                    # Fitting model and computing pvalues
                                    fitlimma = lmFit(data, design = model.matrix(~factor(design)))
                                    fitbayes = eBayes(fitlimma)
                                    DEG.pval = fitbayes$p.value[, 2]
                                    
                                    # Dataframe with adjusted pvalues
                                    DEG = data.frame( pval = p.adjust(DEG.pval, method = "BH"),
                                                      Gene.ID = row.names(count.matrix.raw)
                                    )
                                    colnames(DEG) = c(method,"SYMBOL")
                                    return(DEG)
                                    
                                  },
                                  stop("Enter something that switches me !")
  )
  
  if (tool.DEG%in% c("nbinom.Wald","nbinom.LRT")){
    DEG.pval = results(DEG)$pvalue
  }
  # Recompute pvalues with the false discovery method
  DEG.padj = p.adjust(DEG.pval, method = "BH")
  
  DEG = data.frame(SYMBOL = row.names(count.matrix.raw), DEG.padj)
  colnames(DEG) = c("SYMBOL", method)
  
  return(DEG)
}



#' Merge the DEG Pvalues for each tool used by tools.DEG.RNAseq in one dataframe
#'
#' @param data 
#' Dataframe with genes in row, and methods used in columns. In contains the differentially expressed p-values for each gene.
#' @param tools list character string.
#' By default all the methods present in tools.DEG.RNAseq() are used.
#' Any tools given in this list could be passed as argument : 
#' "edgeR_RLE","edgeR_upperquartile","edgeR_TMMwsp","deseq2.Wald","deseq2.LRT", "deseq"
#' @return 
#' Dataframe of DEG pvalues with genes in columns and tools in rows.
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute all the possible pvalues
#' # res.DEG = tools.DEG.RNAseq.merge(Data)
#' 
#' # Compute the pvalues, only with edgeR methods
#' tools = c("edgeR_RLE","edgeR_upperquartile","edgeR_TMM")
#' res.DEG = tools.DEG.RNAseq.merge(data = Data, tools = tools)
tools.DEG.RNAseq.merge <- function(data,tools){
  if(missing(tools)){
    tools = c("edgeR_TMM","edgeR_RLE","edgeR_upperquartile","edgeR_TMMwsp","deseq2.Wald","deseq2.LRT", "deseq")
  }

  data_to_comp = data.frame(genes = row.names(data))
  for (tool in tools){
    # for each tool, the corresponding analysis is computed
    print(tool)
    tmp = tools.DEG.RNAseq(data,tool)
    # the pvalue column associated to the tool is added
    data_to_comp = merge(data_to_comp,tmp,by = "genes",all=T)  
  }
  # The gene column is placed as row names
  row.names(data_to_comp) <- data_to_comp$genes
  data_to_comp <- data_to_comp[,-1]
  # we obtain a dataframe with genes in columns and methods in rows
  data_to_comp = as.data.frame(t(data_to_comp))
  return(data_to_comp)
}

#' Merge the DEG Pvalues for each tool used by tools.DEG.Nanostring in one dataframe
#'
#' @param tools_DEG Method to use to compute the pvalues of differentially expressed genes.
#' "Wilcox" uses the wilcoxDEG() function implemented in this very same pacakge 
#' "limma" uses the functions DEG_limma() that comes from the limma package
#' "RankProduct" and "RankSum" perform respectively a Rank Product and a Rank Sum analysis with the RankProducts() function from the RankProd package
#'   
#' @param tools_norm Normalization tool. "nappa.NS", "nappa.param1","nappa.param2","nappa.param3" are different parameters used with the NAPPA() function from the NAPPA package.
#'  "nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2" use the NanoStringNorm() normalization function from the package NanoStringNorm
#'  "nanostringR" uses the HKnorm() function from the package nanostringr.
#'  "nanoR.top100","nanoR.total" uses the nsNormalize() function from the nanoR package.
#'  For the nanoR package, it is needed to give the file path to rcc files
#'  
#' @param DESeq logical values.
#' TRUE to use DESeq2 DEG analysis with its own normalization
#' 
#' 
#' @param raw.data rcc type file. List of 4 elements: 
#' Samples.IDs is a dataframe with four columns: sample name, status (WildType, Mutated), filename. 
#' rcc.df contains expression levels with samples in column and genes in rows. 
#' annots.df is a dataframe with 3 columns: geneclass (endogenous, housekeeping, negative, positive), gene name, and the accession number. 
#' FOV is a dataframe that contains Field of views information.
#' 
#' @param dir directory of rcc files. 
#' This parameter is only necessary if the nanoR normalizations are wanted.
#'
#' @return 
#' Dataframe with genes in columns and as many rows as there is possible combination of methods
#' @export
#'
#' @examples
#' # Retrieve gene expression data from Nanostring
#' # raw.data = Simul.data(type = "Nanostring")
#' 
#' # Compute the DEG pvalues with some of the possible methods
#' tools_DEG = c("Wilcox","limma")
#' tools_norm = c("nappa.NS","nanoR.total")
#' # Give the location of your rcc files
#' #RCC.dir <- file.path("./data/NANOSTRING","GSE146204_RAW")
#' #res.DEG = tools.DEG.Nanostring.merge(raw.data =raw.data,
#' #                                      tools_DEG = tools_DEG,
#' #                                      tools_norm = tools_norm,
#' #                                      DESeq = FALSE,
#' #                                      dir = RCC.dir)
#' 
#' # Compute the DEG pvalues with only "nappa.default" normalization 
#' #and Wilcoxon's test only
#' #res.DEG = tools.DEG.Nanostring.merge(raw.data =raw.data,
#' #                                      "Wilcox",
#' #                                      "nappa.default",
#' #                                      DESeq = FALSE)
tools.DEG.Nanostring.merge <- function(raw.data,tools_DEG,tools_norm,DESeq=T,dir = NULL){
  # By default, all the possible normalization methods are computed
  if (missing(tools_norm)){
    tools_norm <- c("nappa.NS","nappa.param1", "nappa.param2","nappa.param3","nanostringnorm.default","nanostringnorm.param1","nanostringnorm.param2","nanoR.top100","nanoR.total","nanostringR")
  }
  # Same for the DEG analysis methods
  if (missing(tools_DEG)){
    tools_DEG = c("limma","Wilcox","RankProduct","RankSum" )
  }
  # putting the directory into another variable
  RCC.dir = dir
  # To create the dataframe, it is easier to start with DESeq2 since this analysis uses an integrated normalization method
  # Otherwise, a dataframe with one column with genes is conserved to merge further methods with it
  if(DESeq){
    print("DESeq2")
    data.to.comp <- tools.DEG.Nanostring(raw.data = raw.data,data =  raw.data, tool = "desq2", tool_norm = NULL)
  }
  else{
    data.to.comp = data.frame(SYMBOL = row.names(raw.data$rcc.df))
  }
  # For each normalization tool, the data are modified in order to be analyzed with a DEG method.
  for (tool_norm in tools_norm){
    nanoR=F
    raw <- raw.data
    dir = NULL
    # nanoR methods needs a particular setup. It needs the directory parameters passed in parameters
    if (tool_norm%in%c("nanoR.top100","nanoR.total")){
      nanoR <- T
      dir = RCC.dir
    }
    # the normalized dataset is computed
    tmp <- tools.norm.Nanostring(raw,tool_norm,nanoR = nanoR, dir = dir)
    # Now it can be analyzed with each DEG tools
    for (tool_diff in tools_DEG){
      print(paste(c(tool_norm, tool_diff)))
      tmp2 = tools.DEG.Nanostring(raw.data = raw.data, data =  tmp, tool = tool_diff, tool_norm = tool_norm)
      # adding the pvalues by merging with the SYMBOL columns that contains the genes
      data.to.comp <- merge(data.to.comp,tmp2,by="SYMBOL",all=T) 
    }
  }
  # The gene column is placed as row names
  row.names(data.to.comp) <- data.to.comp$SYMBOL
  data.to.comp <- data.to.comp[,-1]
  # if you need to remove genes not present in all analysis: NA (missing) data
  data.to.comp <- na.omit(data.to.comp)
  # we obtain a dataframe with genes in columns and methods in rows
  data.to.comp <- as.data.frame(t(data.to.comp))
  return(data.to.comp)
}

#' Merge the DEG Pvalues for each tool used by tools.DEG.Microarrays in one dataframe
#'
#' @param data Dataframe with genes in row, and methods used in columns. It contains the differentially expressed p-values for each gene.
#' @param tools Different methods are used to compute pvalues of differentially expressed genes for microarrays
#' "Wilcox" uses the wilcoxDEG() function implemented in this very same package 
#' "limma" and "GEOlimma uses respectively the functions DEG_limma() and DEG_GEOlimma() that come from the limma package
#' "RankProduct","RankProduct.log" perform a Rank Product analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 
#' "RankSum","RankSum.log" perform a Rank Sum analysis with the RankProducts() function from the RankProd package for normal and logged values respectively 

#' @param n1 Number of samples for the first experimental condition
#' @param n2 Number of samples for the second experimental condition
#'
#' @return Dataframe of pvalues of genes being differentially expressed with genes in columns and methods in rows
#' @export
#'
#' @examples
#' # Import the dataset
#' Data = matrix(runif(5000, 10, 100), ncol=20)
#' group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
#' genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Data) = group
#' row.names(Data) = genes 
#' 
#' # Compute Pvalues for all the methods 
#' tools = c("limma", "Wilcox","RankProduct","RankProduct.log","RankSum","RankSum.log")
#' res.DEG = tools.DEG.Microarrays.merge(Data,tools,10,10)
tools.DEG.Microarrays.merge <-function(data,tools,n1,n2){
  # By default, all the implemented methods are used
  if (missing(tools)){
    tools = c("limma", "GEOlimma", "Wilcox","RankProduct","RankProduct.log","RankSum","RankSum.log")
  }
  if (missing(n1)){
    n1=ncol(data)/2
  }
  if (missing(n2)){
    n2=ncol(data)/2
  }
  # Creating the dataframe with only genes in it to merge further results with it
  data_to_comp = data.frame(Gene.ID = row.names(data))
  
  for (tool in tools){
    print(tool)
    tmp = tools.DEG.Microarrays(data,tool,n1, n2)
    data_to_comp = merge(data_to_comp,tmp,by = "Gene.ID",all=T)  
  }
  # adding the genes as row names
  row.names(data_to_comp) <- data_to_comp$Gene.ID
  data_to_comp <- data_to_comp[,-1]
  # Genes in columns, methods in rows
  data_to_comp = as.data.frame(t(data_to_comp))
  return(data_to_comp)
}

