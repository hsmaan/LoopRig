#' Create a list of ranges objects from looping data
#' 
#' Uses tab delimited loop data in the form of BedPE format
#' @param ... Any number of tab delimited loop data files in BedPE format
#' @param loop_names A character vector of names for the looping datasets (optional)
#' @return A 'LoopRanges' class object: list of GRanges loop data objects 
#' @export 

LoopsToRangesBed <- function(..., loop_names = NULL) {
  
  pre_list <- list(...)
  
  meta_classes <- c("character", "numeric", "numeric", "character", "numeric", "numeric")

  loop_list <- structure(lapply(pre_list, read.table, colClasses = c(meta_classes), header = FALSE), class = "Looplist")
  
  loop_to_range <- function(loop_df) {
    
    anchor_1_ranges <- GRanges(seqnames = loop_df[[1]], ranges = IRanges(start = loop_df[[2]], end = loop_df[[3]]))
    anchor_2_ranges <- GRanges(seqnames = loop_df[[4]], ranges = IRanges(start = loop_df[[5]], end = loop_df[[6]]))
    anchor_list <-  GRangesList(anchor_1_ranges, anchor_2_ranges)
    
  }
  
  l_ranges_list <- structure(lapply(loop_list, loop_to_range), class = "LoopRanges")
  
  if (is.null(loop_names) == FALSE) {
    names(e_ranges_list) <- loop_names
    l_ranges_list
  }
  
  else { 
    l_ranges_list
  }
}