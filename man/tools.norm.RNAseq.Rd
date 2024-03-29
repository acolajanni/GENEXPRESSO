% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalization.R
\name{tools.norm.RNAseq}
\alias{tools.norm.RNAseq}
\title{Compute normalization / size factors for RNAseq count matrix}
\usage{
tools.norm.RNAseq(count.matrix, tool, design)
}
\arguments{
\item{count.matrix}{Dataframe of count with samples in columns and genes SYMBOL in rows.}

\item{tool}{Character string among "TMM","TMMwsp", "RLE", "Upperquartile", "voom", "vst", "vst2".

"TMM","TMMwsp", "RLE", "Upperquartile" calls the \link[edgeR]{calcNormFactors} function.
"voom" calls the \link[limma]{voom} function.
"vst" calls the \link[DESeq]{estimateSizeFactors} function on a CountDataSet.
"vst2" does the same but also calls the \link[DESeq2]{varianceStabilizingTransformation} function.}

\item{design}{Vector of 1 and 2 of the same length of colnames(count.matrix).
1 for the first group and 2 for the second.}
}
\value{
"TMM","TMMwsp", "RLE", "Upperquartile" and "vst" returns a vector of the same size as colnames(count.matrix)
"voom" returns an "Elist" class object.
"vst2" returns the normalized count.matrix with the variance stabilizing transformation
}
\description{
This function calls various function to compute normalization, size factors, or normalized dataset 
to call another function that analyses differentially expressed genes
}
\examples{
# load a count matrix (example with a random dataset)
Data = matrix(runif(5000, 10, 100), ncol=20)
group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
colnames(Data) = group
row.names(Data) = genes 

# Compute design vector
design = c(rep(1,10), rep(2,10)) # 10 from group 1, 10 from group 2

Norm = tools.norm.RNAseq(Data, "TMM", design)

}
