###########################
####### Simulation ########
###########################

#' Creates dataset
#'
#' Creates a dataset depending on the wanted data type 
#'
#' @param type Data type wanted. "RNAseq" Simulates data with SimulateRnaSeqData() from the DEFormats package. "Microarrays" simulates data with the function madsim() from madsim package. "Nanostring" imports an already existing dataset, from real data.
#' @param n.cond1 Sample number in the first group 
#' @param n.cond2 Sample number in the second group
#' @param nb.genes Gene numbers
#'
#' @return 
#' @export 
#'
#' @examples
#' # To get a RNAseq type data with 1000 genes and 2 groups of 15 samples
#' Data = Simul.data(type = "RNAseq", n.cond1 = 15, n.cond2 = 15, nb.genes = 1000)
Simul.data <-function(type,n.cond1,n.cond2,nb.genes){

  library(DEFormats)
  library(madsim)
  
  Choose.type<-switch(type,
                      
                      RNAseq = {
                        # Création de deux jeux de données RNAseq
                        counts1 <- as.data.frame(simulateRnaSeqData (n=nb.genes,m=n.cond1,seed = 222))
                        counts2 <- as.data.frame(simulateRnaSeqData (n=nb.genes,m=n.cond2,seed = 200))
                        
                        # Fusion des deux jeux de données
                        counts=merge(counts1, counts2, by = "row.names", all = TRUE)
                        row.names(counts)=counts$Row.names
                        counts=counts[,-1]
                        
                        # Renommer les échantillons des deux groupes
                        group = paste0(rep(c("control", "case"), each = 15),rep(c(1:15),each = 1))
                        colnames(counts) = group
                        data = as.data.frame(counts)
                        
                        
                      },
                      
                      Nanostring = {
                        
                        
                      },
                      
                      Microarrays = {
                        
                        if(nb.genes < 100){
                          stop("The genes number is too low. \n The minimal value is 100 genes")
                        }
                        
                        # Paramètres de simulation
                        fparams <- data.frame(m1 = n.cond1, m2 = n.cond2, shape2 = 4, lb = 4, ub = 14,pde=0.06,sym=0.5)
                        dparams <- data.frame(lambda1 = 0.13, lambda2 = 2, muminde = 1, sdde = 0.5)
                        sdn <- 0.4
                        rseed <- 50
                        data <- madsim(mdata = NULL, n = nb.genes, ratio = 0, fparams, dparams, sdn, rseed)
                        data =data$xdata
                        
                        # Création des noms des échantillons selon les deux conditions
                        listecol1 = paste0(rep("Control",each = n.cond1),rep(1:n.cond1, each = 1) )
                        listecol2 = paste0(rep("Case",each = n.cond2),rep(1:n.cond2, each = 1) )
                        listecol = c(listecol1,listecol2)
                        colnames(data)<-listecol
                        
                        # Création des noms des gènes
                        listenom <- paste0(rep(LETTERS[1:26], each=400), rep(1:400, 26))
                        listenom = listenom[1:nb.genes]
                        row.names(data) <- listenom
                        
                      },
                      stop("Enter something that switches me!") 
  )
  return(data)
                      
}
