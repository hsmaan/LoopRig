#' Export a LoopRanges or \emph{italics}ElementRanges object to a BED file 
#'
#' Uses ElementRanges objects to output to BED format, and LoopRanges objects to output to BEDPE format
#' @param obj An object of \emph{italics}LoopRanges or \emph{italics}ElementRanges class 
#' @param file A string indicating the name and save location of BED file 
#' @export 


ExportBED <- function(obj, file) {

  if (class(obj) == "LoopRanges") {
    
    
  }
  
  else if (class(obj) == "ElementRanges") {
    
    
  }
  
  else {
    stop("Incorrect object type. Please enter either a 'LoopRanges' or 'ElementRanges' object")
  }
  
}