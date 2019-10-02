#' Determines elements linked by loop anchors 
#' 
#' Using a finalized loops object and two sets of elements, determines which elements are linked by loops through overlap of anchor regions
#' @param loop_ranges A single GRanges loop object in the form of GRangesList class. Should ideally be the consensus loops object created through ConsensusLoops(). 
#' @param element_ranges_1 The first element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param element_ranges_2 The second element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 base-pair.
#' @return Returns a dataframe indicating the unique element links present in the form of their first metadata columns. The first element ranges nodes are always under one column, and the second under the other. 
#' @export

LinkedElements <- function(loop_ranges, element_ranges_1, element_ranges_2, overlap_threshold = 1) {
  
  anchor1_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_1, minoverlap = overlap_threshold))
  anchor2_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_2, minoverlap = overlap_threshold))
  hits_merge_1 <- merge(anchor1_overlaps_1, anchor2_overlaps_1, by = "queryHits")
  
  anchor1_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_1, minoverlap = overlap_threshold))
  anchor2_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_2, minoverlap = overlap_threshold))
  hits_merge_2 <- merge(anchor1_overlaps_2, anchor2_overlaps_2, by = "queryHits")
  
    
  out_df <- data.frame("Anchor 1" = c(mcols(element_ranges_1[hits_merge_1[[2]]])[[1]], mcols(element_ranges_1[hits_merge_2[[2]]])[[1]]), "Anchor 2" = c(mcols(element_ranges_2[hits_merge_1[[3]]])[[1]], mcols(element_ranges_2[hits_merge_2[[3]]])[[1]]))
  
  return(unique(out_df))
    
}