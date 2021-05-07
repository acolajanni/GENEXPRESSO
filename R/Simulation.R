###########################
####### Simulation ########
###########################

library(DEFormats)
library(madsim)

#' Create a dataset
#'
#' Create a dataset depending on the wanted data type 
#'
#' @param type Data type wanted. "RNAseq" simulates data with SimulateRnaSeqData() from the DEFormats package. "Microarrays" simulates data with the function madsim() from madsim package. "Nanostring" imports an already existing dataset, from real data.
#' @param n.cond1 Sample number in the first group 
#' @param n.cond2 Sample number in the second group
#' @param nb.genes Gene numbers
#'
#' @return 
#' Dataframe of gene expression values
#' @docType data
#' @usage data(Data_Nanostring)
#' @import "madsim" "DEFormats"
#' @export 
#' 
#' @examples
#' # To get a RNAseq type data with 1000 genes and 2 groups of 15 samples
#' Data = Simul.data(type = "RNAseq", n.cond1 = 15, n.cond2 = 15, nb.genes = 1000)
#' # To get Microarray type data with 1000 genes and 2 groups of 15 samples
#' Data = Simul.data(type = "Microarrays", n.cond1 = 15, n.cond2 = 15, nb.genes = 1000)
#' # To get the non simulated Nanostring data
#' Data = Simul.data(type = "Nanostring")
Simul.data <-function(type,n.cond1,n.cond2,nb.genes){
  
  Choose.type<-switch(type,
                      
                      RNAseq = {
                        # Creating two different datasets
                        counts1 <- as.data.frame(simulateRnaSeqData (n=nb.genes,m=n.cond1,seed = 222))
                        counts2 <- as.data.frame(simulateRnaSeqData (n=nb.genes,m=n.cond2,seed = 200))
                        
                        # Merging both 
                        counts=merge(counts1, counts2, by = "row.names", all = TRUE)
                        row.names(counts)=counts$Row.names
                        counts=counts[,-1]
                        
                        # Rename the samples
                        group = paste0(rep(c("control", "case"), each = 15),rep(c(1:15),each = 1))
                        colnames(counts) = group
                        data = as.data.frame(counts)
                        
                        
                      },
                      
                      Nanostring = {
                        data = GENEXPRESSO::Data_Nanotring
     
                      },
                      
                      Microarrays = {
                        
                        if(nb.genes < 100){
                          stop("The genes number is too low. \n The minimal value is 100 genes")
                        }
                        
                        # Simulations parameters
                        fparams <- data.frame(m1 = n.cond1, m2 = n.cond2, shape2 = 4, lb = 4, ub = 14,pde=0.06,sym=0.5)
                        dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
                        sdn <- 0.4
                        rseed <- 50
                        data <- madsim(mdata = NULL, n = nb.genes, ratio = 0, fparams, dparams, sdn, rseed)
                        data =data$xdata
                        
                        # Creating samples names
                        listecol1 = paste0(rep("Control",each = n.cond1),rep(1:n.cond1, each = 1) )
                        listecol2 = paste0(rep("Case",each = n.cond2),rep(1:n.cond2, each = 1) )
                        listecol = c(listecol1,listecol2)
                        colnames(data)<-listecol
                        
                        # Creating Gene names
                        listenom <- paste0(rep(LETTERS[1:26], each=400), rep(1:400, 26))
                        listenom = listenom[1:nb.genes]
                        row.names(data) <- listenom
                        data = data.frame(data)
                        
                      },
                      stop("Enter something that switches me!") 
  )
  return(data)
                      
}