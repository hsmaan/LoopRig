#' Finds elements stacked atop loop anchors
#' 
#' Using a finalized \emph{LoopRanges} loops object and two sets of elements, determines which pairs of elements from the two sets are found overlapping with each other and a loop anchor.
#' @param loop_ranges A single \emph{LoopRanges} object that has been subset to one range through ConsensusLoops().  
#' @param element_ranges_x The first element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param element_ranges_y The second element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param range_out_x A boolean indicating if the output should be a subset element ranges object for element_ranges_x instead of the default dataframe. (default = FALSE) 
#' @param range_out_y A boolean indicating if the output should be a subset element ranges object for element_ranges_y instead of the default dataframe. (default = FALSE) 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 base-pair.
#' @return Returns a dataframe indicating stacked elements in the form of their first metadata columns. The element_ranges_x are always under column 1, and element_ranges_y are always under column 2. 
#' @import GenomicRanges
#' @import IRanges
#' @export


StackedElements <- function(loop_ranges, element_ranges_x, element_ranges_y, range_out_x = FALSE, range_out_y = FALSE, overlap_threshold = 1) {
  
  if(class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the loop_ranges parameter")
  }
  
  if(length(loop_ranges) != 1) {
    stop("Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  }
  
  if((range_out_x == TRUE) & (range_out_y == TRUE)) {
    stop("Can only output either element_ranges_x or element_ranges_y as ElementRanges object")
  }
  
  loops_collapsed <- unlist(loop_ranges[[1]])
  
  element1_overlaps <- as.data.frame(findOverlaps(loops_collapsed, element_ranges_x, minoverlap = overlap_threshold))
  element2_overlaps <- as.data.frame(findOverlaps(loops_collapsed, element_ranges_y, minoverlap = overlap_threshold))
  hits_merge <- merge(element1_overlaps, element2_overlaps, by = "queryHits")
  
  if(range_out_x == TRUE) {
    
    range_x_subset <- element_ranges_x[hits_merge[[2]]]
    return(structure(list("el_x_stacked" = range_x_subset), class = "ElementRanges"))
  }
  
  if (range_out_y == TRUE) {
    
    range_y_subset <- element_ranges_y[hits_merge[[3]]]
    return(structure(list("el_y_stacked" = range_y_subset), class = "ElementRanges"))
  }
  
  out_df <- data.frame("Element_1" = mcols(element_ranges_x[hits_merge[[2]]])[[1]], "Element_2" = mcols(element_ranges_y[hits_merge[[3]]])[[1]])
  return(unique(out_df))
    
}