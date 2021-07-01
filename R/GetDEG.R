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
#' @param sets 
#' Specific sets to look at (Include as combinations. Ex: c("Name1", "Name2") 
#' where Name1 and Name2 are two columns of the binary.matrix .
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
#' 
#' # Retrieving all the non redundant predicted down and upregulated genes  
#' DEG = Get.DEG(Binary.matrix,alternative = TRUE, method = "union")
#' 
#' # Retrieving common DE genes for a set of method 
#' subset.DEG = Get.DEG(Binary.matrix,sets = methods[1:4],alternative = TRUE, method = "intersect")
#' 
Get.DEG = function(binary.matrix, sets, method, alternative){
  # By default the binary matrix is composed of Upregulated and downregulated predicted DEG
  # In this case, we have the alternative hypothesis
  
  if (!missing(sets)){
    if (any(!sets %in%  colnames(binary.matrix)) ){
      message(sets[(!sets %in% colnames(binary.matrix))], " is not an known column")
      stop("Enter a valid column name")
    }
    binary.matrix = binary.matrix[,sets]
  }
  
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

