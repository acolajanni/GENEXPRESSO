% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEG.R
\docType{data}
\name{DEG_GEOlimma}
\alias{DEG_GEOlimma}
\title{Differentially expressed genes analysis through the limma package, with the GEOlimma method.}
\usage{
data(GEOlimma_probabilities)
}
\arguments{
\item{dataset}{dataframe of expression values with samples in columns and genes in row.}

\item{design}{vector of 0 and 1 values. 0 for the first experimental condition, 1 for the second one.}

\item{contrast.matrix}{design matrix like the one produced by the make_designMatrix() function.}
}
\value{
A dataframe with 3 columns is returned. It contains the Log Fold change value (logFC),
The Pvalue associated with this Pvalue (PValue)
the last column corresponds to the gene names given in the dataset (SYMBOL)
}
\description{
Uses the \link[limma]{lmFit} and \link[limma]{eBayes} from the limma package to automate the steps of the analysis.
}
\examples{
# Import the dataset
Data = matrix(runif(5000, 10, 100), ncol=20)
group = paste0(rep(c("control", "case"), each = 10),rep(c(1:10),each = 1))
genes <- paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
colnames(Data) = group
row.names(Data) = genes  

#Build the design matrix
#design = make_designMatrix(dataset = Data,ncond1 = 10, ncond2 = 10)
#get the results of limma analysis
#res.diff = DEG_GEOlimma(dataset = Data,design)
}
