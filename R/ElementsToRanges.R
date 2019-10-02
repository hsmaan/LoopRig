#' Create a list of ranges objects from element data
#' 
#' Uses tab delimited element data in the form of BED4 or BED12 files to create GRanges element objects.
#' @param ... Any number of tab delimited element data files in BED4 or BED12 format  
#' @param file_format The type of input element file - either 'BED4' or 'BED12' format
#' @param element_names A character vector of names for the element datasets (optional)
#' @return A 'ElementRanges' class object: list of GRanges element data objects 
#' @export 

ElementsToRanges <- function(..., file_format, element_names = NULL) {
  
  pre_list <- list(...)
  
  if (length(pre_list) == 0) {
    stop("Please enter at least one BED4/12 element data file")
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
    stop("Error in reading BED file(s), please ensure file specification and format accuracy")
  }
  
  meta_classes <- c("character", "numeric", "numeric", "character")
  
  if (file_format == "BED12") {
    el_list <- structure(lapply(pre_list, read.table, colClasses = c(meta_classes,rep("NULL", 8)), header = FALSE), class = "ElementsList")
  } 
  
  else if (file_format == "BED4") {
    el_list <- structure(lapply(pre_list, read.table, colClasses = c(meta_classes), header = FALSE), class = "ElementsList")
  } 
  
  else {
    stop("Please enter an appropriate format for the element files (BED4/BED12)")
  }
  
  
  element_to_range <- function(element_df) {
    element_ranges <- GRanges(seqnames = element_df[[1]], ranges = IRanges(start = element_df[[2]], end = element_df[[3]]), mcols = element_df[[4]])
    element_ranges
  }
    
  e_ranges_list <- structure(lapply(el_list, element_to_range), class = "ElementRanges")
  
  if (is.null(element_names) == FALSE) {
    names(e_ranges_list) <- element_names
    e_ranges_list
  } 
  
  else { 
    e_ranges_list
  }
}