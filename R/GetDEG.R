#' Extract DEG from a binary matrix.
#'
#' Analyze the output of Upset.Binary.Dataframe() to obtain DEG for one method 
#'
#' @param binary.matrix Dataframe or matrix containing only 0 and 1 values
#' @param method Character string that corresponds to any columns of the binary matrix.
#' @param alternative Logical value. 
#' If TRUE, Alternative hypothesis will be searched : for both up and downregulated. 
#' In this case, the method argument should be the name of the analysis (DESeq2 to extract up and downregulated genes).
#' If FALSE, only the method given in argument is analyzed.
#' 
#' @return A list of one vector if alternative = FALSE, 2 vectors if alternative = TRUE
#' @import "stringr" "dplyr" 
#' @export
#'
#' @examples
#' # Getting a binary matrix, DEG are designed by "1" 
#' # 250 genes for 10 methods (up and downregulated)
#' Data = matrix(runif(5000, 0, 1), ncol=20)
#' Binary.matrix = round(Data)
#' methods = paste0(rep("method", each = 20),rep(c(1:10),each = 2),rep(c("_Up","_Down"),each=1))
#' genes = paste0(rep(LETTERS[1:25], each=10), rep(c(1:10),each = 1))
#' colnames(Binary.matrix) = methods
#' row.names(Binary.matrix) = genes
#' 
#' # Retriving both up and down regulated genes for the "method5"
#' DEG = Get.DEG(Binary.matrix,alternative = TRUE, method = "method5")
Get.DEG <- function(binary.matrix, method, alternative){
  # By default argument : 
  # We get the results for upregulated and down regulate genes
  if (missing(alternative)){
    alternative = TRUE
  }
  
  if (alternative){
    # Series of logical value wether or not the "method" argument is contained in the names of columns 
    pattern = str_detect(colnames(binary.matrix),method)
    if (!TRUE%in% pattern){
      stop("Unknown method")
    }
    # actual colnames corresponding to the pattern of "method" is extracted
    # example : 
    # method = "DESeq2"
    # method.name = "DESeq2_Up", "DESeq2_Down"
    method.name = colnames(binary.matrix)[pattern]
  }
  else{
    if (!method%in% colnames(binary.matrix)){
      stop("Unknown method")
    }
    method.name = method
  }
  binary.matrix = as.data.frame(binary.matrix)
  # Removing all the other methods
  result = select(binary.matrix, all_of(method.name))
  result$genes = row.names(result)
  
  # if alternative hypothesis does not interest us
  if (!alternative){
    # Only one column (given in argument) is retrieved
    DEG = subset(result, result[[method.name]] == 1)
    DEG = list(DEG = DEG[["genes"]])
  }
  else{
    # If we want alternative hypothesis, two dataframe are needed (up and downregulated)
    DEG_up = data.frame(Upregulated = result[[ method.name[1] ]], genes = result$genes)
    # only genes with "1" are differentially expressed, other does not interest us
    DEG_up = filter(DEG_up, Upregulated == 1 )
    
    DEG_down = data.frame(Downregulated = result[[ method.name[2] ]], genes = result$genes)
    DEG_down = filter(DEG_down, Downregulated == 1 )
    
    # A list of 2 vectors is returned : for both up and downregulated hypothesis
    DEG = list(Upregulated = DEG_up$genes , Downregulated = DEG_down$genes)
  }
  return(DEG)
}