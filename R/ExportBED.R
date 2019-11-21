#' Export a LoopRanges or \emph{italics}ElementRanges object to a BED/BEDPE file 
#'
#' Uses ElementRanges objects to output to BED format, and LoopRanges objects to output to BEDPE format
#' @param obj An object of \emph{italics}LoopRanges or \emph{italics}ElementRanges class 
#' @param index List index of LoopRanges or ElementRanges object to output
#' @param mcol A boolean specifying whether the first mcol of the object are to be output (default=FALSE)
#' @param file_name A string indicating the name and save location of the BED/BEDPE file 
#' @importFrom utils write.table
#' @export 


ExportBED <- function(obj, index = NULL, mcol = FALSE, file_name) {
  
  if (is.null(index) == TRUE || index <= 0) {
    stop("Please enter a positive integer indicating the index of LoopRanges or ElementRanges object to export")
  }
  
  if(file.exists(file_name) == TRUE) {
    stop("File already exists, cannot overwrite")
  }
  
  ### Dir.exist conditional ###
  
  if (class(obj) == "LoopRanges") {
    
    if (mcol == TRUE) {
      
      if(ncol(mcols(obj[[index]][[1]])) == 0){
        stop("Object has no mcols, please choose mcol = FALSE")
      }
      
      ranges <- obj[[index]]
      range_table <- data.frame("chr1" = as.vector(seqnames(ranges[[1]])), "str1" = start(ranges[[1]]), "end1" = end(ranges[[1]]), "chr2" = as.vector(seqnames(ranges[[2]])), "str2" = start(ranges[[2]]), "end2" = end(ranges[[2]]), "mcol" = mcols(ranges[[1]])[,1])
      write.table(range_table, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("BEDPE file exported to ", file_name, sep = ""))
      
    } else {
      
      ranges <- obj[[index]]
      range_table <- data.frame("chr1" = as.vector(seqnames(ranges[[1]])), "str1" = start(ranges[[1]]), "end1" = end(ranges[[1]]), "chr2" = as.vector(seqnames(ranges[[2]])), "str2" = start(ranges[[2]]), "end2" = end(ranges[[2]]))
      write.table(range_table, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("BEDPE file exported to ", file_name, sep = ""))
      
    }
  }
  
  else if (class(obj) == "ElementRanges") {
    
    if (mcol == TRUE) {
      
      if(ncol(mcols(obj[[index]])) == 0){
        stop("Object has no mcols, please choose mcol = FALSE")
      }
      
      ranges <- obj[[index]]
      range_table <- data.frame("chr" = as.vector(seqnames(ranges)), "str" = start(ranges), "end" = end(ranges), "mcol" = mcols(ranges)[,1])
      write.table(range_table, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("BED file exported to ", file_name, sep = ""))
      
      
    } else {
      
      ranges <- obj[[index]]
      range_table <- data.frame("chr" = as.vector(seqnames(ranges)), "str" = start(ranges), "end" = end(ranges))
      write.table(range_table, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      print(paste("BED file exported to ", file_name, sep = ""))
      
    }
  }
  
  else {
    stop("Incorrect object type. Please enter either a 'LoopRanges' or 'ElementRanges' object")
  }
}

