#' Determines elements linked by loop anchors 
#' 
#' Using a finalized \emph{LoopRanges} loops object and two sets of elements, determines which elements are linked by loops through overlap of anchor regions
#' @param loop_ranges A single \emph{LoopRanges} object that has been subset to one range through ConsensusLoops().  
#' @param element_ranges_1 A subset \emph{ElementRanges} class object. Subset appropriate ranges by using list indexing. 
#' @param element_ranges_2 A subset \emph{ElementRanges} class object. Subset appropriate ranges by using list indexing. 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 bp. 
#' @return Returns a dataframe indicating the unique element links present in the form of their first metadata columns. The first element ranges nodes are always under one column 1, and the second under column 2.
#' @import GenomicRanges
#' @import IRanges
#' @export

LinkedElements <- function(loop_ranges, element_ranges_1, element_ranges_2, overlap_threshold = 1) {
  
  if(class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the loop_ranges parameter")
  }
  
  if(length(loop_ranges) != 1) {
    stop("Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  }
  
  loop_ranges <- loop_ranges[[1]]
  
  anchor1_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_1, minoverlap = overlap_threshold))
  anchor2_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_2, minoverlap = overlap_threshold))
  hits_merge_1 <- merge(anchor1_overlaps_1, anchor2_overlaps_1, by = "queryHits")
  
  anchor1_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_1, minoverlap = overlap_threshold))
  anchor2_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_2, minoverlap = overlap_threshold))
  hits_merge_2 <- merge(anchor1_overlaps_2, anchor2_overlaps_2, by = "queryHits")
    
  out_df <- data.frame("Anchor_1" = c(mcols(element_ranges_1[hits_merge_1[[2]]])[[1]], mcols(element_ranges_1[hits_merge_2[[2]]])[[1]]), "Anchor_2" = c(mcols(element_ranges_2[hits_merge_1[[3]]])[[1]], mcols(element_ranges_2[hits_merge_2[[3]]])[[1]]))
  
  return(unique(out_df))
    
}
