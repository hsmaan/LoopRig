#' Finds elements within loops scaffolded by another set of elements
#' 
#' Using a finalized \emph{LoopRanges} loops object and two sets of elements, determines which loops are 'anchored' by the first set of elements (x) and finds which elements from the second set (y) are within these loops. 
#' @param loop_ranges A single \emph{LoopRanges} object that has been subset to one range through ConsensusLoops().  
#' @param element_ranges_x The first element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param element_ranges_y The second element ranges to be considered. If using an 'ElementRanges' class object, subset by using list indexing. 
#' @param range_out_x A boolean indicating if the output should be a subset element ranges object for element_ranges_x instead of the default dataframe. (default = FALSE) 
#' @param range_out_y A boolean indicating if the output should be a subset element ranges object for element_ranges_y instead of the default dataframe. (default = FALSE) 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 base-pair.
#' @return Returns a dataframe indicating scaffold linked elements in the form of their first metadata columns and the loop IDs for their scaffolds. The element_ranges_x (scaffold) are always under column 1, and element_ranges_y (links) are always under column 2.
#' @import GenomicRanges
#' @import IRanges
#' @export 

ScaffoldElements <- function(loop_ranges, element_ranges_x, element_ranges_y, range_out_x = FALSE, range_out_y = FALSE, overlap_threshold = 1) {
  
  if(class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the loop_ranges parameter")
  }
  
  if(length(loop_ranges) != 1) {
    stop("Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  }
  
  if((range_out_x == TRUE) & (range_out_y == TRUE)) {
    stop("Can only output either element_ranges_x or element_ranges_y as ElementRanges object")
  }
  
  loops_collapsed <- loop_ranges[[1]]
  
  loops_collapsed[[1]]$id <- seq(1, length(loops_collapsed[[1]]))
  loops_collapsed[[2]]$id <- seq(1, length(loops_collapsed[[1]]))
  
  scaff_anchors <- mergeByOverlaps(element_ranges_x, unlist(loops_collapsed), minoverlap = overlap_threshold)
  
  loop_a1_subset <- loops_collapsed[[1]][loops_collapsed[[1]]$id %in% scaff_anchors$id]
  loop_a2_subset <- loops_collapsed[[2]][loops_collapsed[[2]]$id %in% scaff_anchors$id]
  
  loops_subset <- GRangesList(loop_a1_subset, loop_a2_subset)
  names(loops_subset) <- c("Anchor 1", "Anchor 2")
  
  loops_full_circ <- GRanges(seqnames = as.vector(seqnames(loops_subset[[1]])), ranges = IRanges(start = start(loops_subset[[1]]), end = end(loops_subset[[2]])))
  
  loops_full_circ$id <- loops_subset[[1]]$id
  
  element_loops <- mergeByOverlaps(element_ranges_y, loops_full_circ, minoverlap = overlap_threshold)
  
  joined_df <- merge(scaff_anchors[,c(2,4)], element_loops[,c(2,4)], by = "id", all = TRUE)
  colnames(joined_df) <- c("Loop ID", "Element_1_Scaffold", "Element_2_Link")
  
  if(range_out_x == TRUE) {
    
    range_x_subset <- unique(element_ranges_x[which(mcols(element_ranges_x)[,1] %in% joined_df[,2])])
    return(structure(list("el_x_scaffold" = range_x_subset), class = "ElementRanges"))
  }
  
  if (range_out_y == TRUE) {
    
    range_y_subset <- unique(element_ranges_y[which(mcols(element_ranges_y)[,1] %in% joined_df[,3])])
    return(structure(list("el_y_scaffold_linked" = range_y_subset), class = "ElementRanges"))
  }
  
  return(unique(as.data.frame(joined_df[,c(2,3)])))
  
}

