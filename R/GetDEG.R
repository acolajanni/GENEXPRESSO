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
  # We get the results for both upregulated and down regulate genes
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







#' Extract DEG from a binary matrix.
#'
#' Analyze the output of Upset.Binary.Dataframe() to obtain Dthe union or the intersection of predicted DEG of each method.
#'
#' @param binary.matrix Dataframe or matrix containing only 0 and 1 values, with genes in rows and methods in columns
#' @param method Character string 
#' "union" or "intersect", "union" to get all the predicted DEG and "intersect" for the intersection of predicted DEG among the set of methods
#' @param alternative Logical value. 
#' TRUE to compute the union or intersection of predicted DEG if the binary matrix is composed of strong hypothesis prediction (for Upregulated and Downregulated genes)
#' FALSE if the dataframe contains no such information, just weak hypothesis prediction 
#'
#' @import "stringr" "dplyr"
#' @return A list of two vectors of character string is returned for alternative = TRUE : for upregulated and downregulated genes. Otherwise, only a vector of DEG is returned 
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
#' # Retrieving all the predicted down and upregulated genes  
#' DEG = Get.DEG.2(Binary.matrix,alternative = TRUE, method = "union")
Get.DEG.2 = function(binary.matrix, method, alternative){
  # By default the binary matrix is composed of Upregulated and downregulated predicted DEG
  # In this case, we have the alternative hypothesis
  if(missing(alternative)){
    alternative = TRUE
  }
  # Union, all the genes predicted
  if (method == "union"){ 
    threshold = 1
    # Intersect = genes found by all the methods  
  } else if (method == "intersect") {
    threshold = ncol(binary.matrix)
  } else {
    stop("Give me a method that switches me !")
  }

  # the dataframe is composed of 0 and 1 values, if the row sums is greater than the threshold,
  # it means that it is a gene that is predicted as DEG
  if (!alternative){
    
    if (TRUE %in% str_detect(colnames(binary.matrix), "Up") || TRUE %in% str_detect(colnames(binary.matrix),"Down") ){
      warning("Upregulated and Downregulated predicted genes are present in your binary matrix. You may want to pass alternative = TRUE in parameters.")
    }
    
    res = binary.matrix[rowSums(binary.matrix) >= threshold , ]
    # In this case, only DEG in weak hypothesis is given in argument 
    # So only a vector of DEG is returned
    genes = row.names(res)
  }  else{
    # Otherwise, the binary matrix is subsetted in Up and Down method
    Up = binary.matrix[,str_detect(colnames(binary.matrix),"Up")]
    Down = binary.matrix[,str_detect(colnames(binary.matrix),"Down")]
    
    # threshold/2 because Up and Down should be the same size
    res.up = Up[rowSums(Up) >= (threshold/2) , ]
    res.down = Down[rowSums(Down) >= (threshold/2) , ]
    
    # a list of 2 vectors is returned : one for downregulated and the second for upregulated genes
    genes = list(Upregulated = row.names(res.up),Downregulated = row.names(res.down) )
  }
  
  if (length(genes) == 0){
    message("The intesection of predicted DEG is null")
  }
  
  return(genes)
}
