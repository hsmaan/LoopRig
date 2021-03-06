#' Create a list of ranges objects from element data
#' 
#' Uses tab delimited element data in the form of BED4 or BED12 files to create GRanges element objects
#' @param ... Any number of tab delimited element data files in BED4 or BED12 format  
#' @param element_names A character vector of names for the element datasets (optional)
#' @param custom_cols An integer indicating the number of extra columns in the BED file.
#' @param custom_mcols An integer or vector of integers indicating which columns are used for metadata (optional) 
#' @return An \emph{ElementRanges} class object: list of GRanges element data objects 
#' @examples 
#' # Load enhancer and promoter elements into an ElementRanges object
#' enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
#' promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)
#' element_ranges <- ElementsToRanges(enhancers, promoters, 
#' element_names = c("enhancers", "promoters"), 
#' custom_cols = 1, custom_mcols = 4)
#' element_ranges
#' @import GenomicRanges
#' @import IRanges
#' @importFrom utils read.table
#' @export 

ElementsToRanges <- function(..., element_names = NULL, custom_cols = NULL, custom_mcols = NULL) {
  
  pre_list <- tryCatch(
    {
      pre_list <- list(...)
    },
    
    warning=function(warning_text) {
      return(NA)
    }
  )
  
  if (length(pre_list) == 0) {
    stop("Please enter at least one BED4/12 element data file")
  }
  
  if(length(unique(pre_list)) != length(pre_list)) {
    stop("Duplicate BED files entered, only unique entries allowed")
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
    stop("Error in reading BED files, please ensure file specification and format accuracy")
  }
  
  if (length(which(sapply((test_table[1, c(2, 3)]), is.numeric))) < 2) {
    stop("Chromosome coordinates are not numeric. Please ensure file format accuracy and remove any header columns")
  }
  
  meta_classes <- c("character", "numeric", "numeric")
  
  if (is.null(custom_cols) == TRUE) {
    stop("Please enter an appropriate value for the custom_cols parameter")
  } 
  
  else {
    
    if (length(colnames(test_table)) != (3 + custom_cols)) {
      stop("Incorrect number of custom columns, please check your files again")
    }
    
    el_list <- structure(lapply(pre_list, read.table, colClasses = c(meta_classes,rep("character", custom_cols)), header = FALSE), class = "ElementsList")
    
  }
  
  
  element_to_range <- function(element_df) {
    
    if (is.null(custom_mcols) == TRUE) {
      
      element_ranges <- GRanges(seqnames = element_df[, 1], ranges = IRanges(start = element_df[, 2], end = element_df[, 3]))
      element_ranges
    } 
    
    else {
      
      if (!(custom_mcols %in% c(4:(4 + custom_cols))) == TRUE) {
        stop("Incorrect custom mcols location, please check your files again")
      }
      
      element_ranges <- GRanges(seqnames = element_df[, 1], ranges = IRanges(start = element_df[, 2], end = element_df[, 3]), mcols = element_df[, custom_mcols])
      element_ranges
    }
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