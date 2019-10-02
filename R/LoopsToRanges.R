#' Create a list of ranges objects from looping data
#' 
#' Uses tab delimited looping data in the form of BEDPE files to create LoopList objects
#' @param ... Any number of tab delimited loop data files in BEDPE format 
#' @param loop_names A character vector of names for the loop datasets (optional)
#' @param custom_cols An integer indicating the number of extra columns in the BEDPE file 
#' @return A 'LoopRanges' class object: list of GRanges looping data objects 
#' @export 

LoopsToRanges <- function(..., loop_names = NULL, custom_cols = 0) {
  
  pre_list <- list(...)
  
  if (length(pre_list) == 0) {
    stop("Please enter at least one BEDPE looping data file")
  }
  
  test_table <- tryCatch(
    {
    read.table(pre_list[[1]], header = FALSE)
    },
    error=function(error_text) {
      return(NA)
    }
  )
  
  if (class(test_table) != "data.frame") {
    stop("Error in reading BEDPE file(s), please ensure file specification and format accuracy")
  }
  
  if (length(colnames(test_table)) != (6 + custom_cols)) {
    stop("Incorrect number of columns, please check your files again")
  }
  
  loop_classes <- c("character", "numeric", "numeric", "character", "numeric", "numeric")
  
  l_list <- structure(lapply(pre_list, read.table, colClasses = c(loop_classes, rep("NULL", custom_cols)), header = FALSE), class = "LoopList")

  loop_to_range <- function(loop_df) {
    anchor_1_ranges <- GRanges(seqnames = loop_df[[1]], ranges = IRanges(start = loop_df[[2]], end = loop_df[[3]]))
    anchor_2_ranges <- GRanges(seqnames = loop_df[[4]], ranges = IRanges(start = loop_df[[5]], end = loop_df[[6]]))
    anchor_list <-  GRangesList("Anchor 1" = anchor_1_ranges, "Anchor 2" = anchor_2_ranges)
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
