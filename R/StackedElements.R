#' Finds elements stacked atop loop anchors
#' 
#' Using a finalized loops object and two sets of elements, determines which pairs of elements from the two sets are found overlapping with each other and a loop anchor.
#' @param loop_ranges A single GRanges loop object in the form of GRangesList class. Should ideally be the consensus loops object created through ConsensusLoops(). 
#' @param element_ranges_1 The first element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param element_ranges_2 The second element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 base-pair.
#' @return Returns a dataframe indicating stacked elements in the form of their first metadata columns. The first element ranges nodes are always under one column, and the second under the other. 
#' @export


StackedElements <- function(loop_ranges, element_ranges_1, element_ranges_2, overlap_threshold = 1) {
  
  loops_collapsed <- unlist(loop_ranges)
  
  element1_overlaps <- as.data.frame(findOverlaps(loops_collapsed, element_ranges_1, minoverlap = overlap_threshold))
  element2_overlaps <- as.data.frame(findOverlaps(loops_collapsed, element_ranges_2, minoverlap = overlap_threshold))
  hits_merge <- merge(element1_overlaps, element2_overlaps, by = "queryHits")
  
  out_df <- data.frame("Element 1" = mcols(element_ranges_1[hits_merge[[2]]])[[1]], "Element 2" = mcols(element_ranges_2[hits_merge[[3]]])[[1]])
    
  return(unique(out_df))
    
}