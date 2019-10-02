#' Create a list of ranges objects from looping data
#' 
#' Uses tab delimited looping data in the form of BED12 files to create GRanges loop objects
#' @param ... Up to 10 tab delimited looping data files in HiCCUPS output format 
#' @param loop_names A character vector of names for the loop datasets (optional)
#' @return A 'LoopRanges' class object: list of GRanges looping data objects 
#' @export 

LoopsToRanges <- function(..., loop_names = NULL) {
  
  pre_list <- list(...)

  loop_classes <- c("character", "numeric", "numeric", "character", "numeric", "numeric")
  
  l_list <- structure(lapply(pre_list, read.table, colClasses = c(loop_classes, "NULL", "numeric", rep("NULL", 12)), header = TRUE), class = "LoopList")

  loop_to_range <- function(loop_df) {
    anchor_1_ranges <- GRanges(seqnames = paste("chr", loop_df[[1]], sep =""), ranges = IRanges(start = loop_df[[2]], end = loop_df[[3]]), mcols = loop_df[[7]])
    anchor_2_ranges <- GRanges(seqnames = paste("chr", loop_df[[4]], sep =""), ranges = IRanges(start = loop_df[[5]], end = loop_df[[6]]), mcols = loop_df[[7]])
    
    anchor_list <-  GRangesList(anchor_1_ranges, anchor_2_ranges)


  }
    
    ranges_list <- structure(lapply(l_list, loop_to_range), class = "LoopRanges")
    
    if (is.null(loop_names) == FALSE) {
      names(ranges_list) <- loop_names
      ranges_list
    } 
    
    else { 
      ranges_list
    }
}
