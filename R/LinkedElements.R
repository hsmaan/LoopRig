#' Determines elements linked by loop anchors 
#' 
#' Using a finalized \emph{LoopRanges} loops object and two sets of elements, determines which elements are linked by loops through overlap of anchor regions
#' @param loop_ranges A single \emph{LoopRanges} object that has been subset to one range through ConsensusLoops().  
#' @param element_ranges_x A subset \emph{ElementRanges} class object. Subset appropriate ranges by using list indexing. 
#' @param element_ranges_y A subset \emph{ElementRanges} class object. Subset appropriate ranges by using list indexing. 
#' @param range_out_x A boolean indicating if the output should be a subset element ranges object for element_ranges_x instead of the default dataframe. (default = FALSE) 
#' @param range_out_y A boolean indicating if the output should be a subset element ranges object for element_ranges_y instead of the default dataframe. (default = FALSE) 
#' @param overlap_threshold Single numerical input for significant base-pair overlap to be considered 'overlapping'. Default is 1 bp. 
#' @return Returns a dataframe indicating the unique element links present in the form of their first metadata columns. The element_ranges_x are always under column 1, and element_ranges_y are always under column 2. 
#' @examples 
#' # Load enhancer and promoter elements into an ElementRanges object
#' enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
#' promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)
#' element_ranges <- ElementsToRanges(enhancers, promoters, 
#' element_names = c("enhancers", "promoters"), 
#' custom_cols = 1, custom_mcols = 4)
#'
#' # Load loops into LoopRanges object and determine consensus loops (keep all)
#' ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", 
#' package = "LoopRig", mustWork = TRUE)
#' pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", 
#' package = "LoopRig", mustWork = TRUE)
#' loops <- LoopsToRanges(ovary_loops, pancreas_loops, custom_cols = 0)
#' consensus_loops <- ConsensusLoops(loops, keep_all = TRUE)
#' 
#' # Based on consensus loops, determine which enhancers and promoters 
#' # are connected based on loop anchors 
#' LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]])  
#' @import GenomicRanges
#' @import IRanges
#' @export

LinkedElements <- function(loop_ranges, element_ranges_x, element_ranges_y, range_out_x = FALSE, range_out_y = FALSE, overlap_threshold = 1) {
  
  if(class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the loop_ranges parameter")
  }
  
  if(length(loop_ranges) != 1) {
    stop("Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  }
  
  if((range_out_x == TRUE) & (range_out_y == TRUE)) {
    stop("Can only output either element_ranges_x or element_ranges_y as ElementRanges object")
  }
  
  loop_ranges <- loop_ranges[[1]]
  
  anchor1_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_x, minoverlap = overlap_threshold))
  anchor2_overlaps_1 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_y, minoverlap = overlap_threshold))
  hits_merge_1 <- merge(anchor1_overlaps_1, anchor2_overlaps_1, by = "queryHits")
  
  anchor1_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[2]], element_ranges_x, minoverlap = overlap_threshold))
  anchor2_overlaps_2 <- as.data.frame(findOverlaps(loop_ranges[[1]], element_ranges_y, minoverlap = overlap_threshold))
  hits_merge_2 <- merge(anchor1_overlaps_2, anchor2_overlaps_2, by = "queryHits")
  
  if(range_out_x == TRUE) {

    range_x_hits <- unique(c(hits_merge_1[[2]], hits_merge_2[[2]]))
    range_x_subset <- unique(element_ranges_x[range_x_hits])
    return(structure(list("el_x_linked" = range_x_subset), class = "ElementRanges"))
  }
  
  if (range_out_y == TRUE) {
    
    range_y_hits <- unique(c(hits_merge_1[[3]], hits_merge_2[[3]]))
    range_y_subset <- unique(element_ranges_y[range_y_hits])
    return(structure(list("el_y_linked" = range_y_subset), class = "ElementRanges"))
  }
  
  out_df <- data.frame("Anchor_1" = c(mcols(element_ranges_x[hits_merge_1[[2]]])[[1]], mcols(element_ranges_x[hits_merge_2[[2]]])[[1]]), "Anchor_2" = c(mcols(element_ranges_y[hits_merge_1[[3]]])[[1]], mcols(element_ranges_y[hits_merge_2[[3]]])[[1]]))
  return(unique(out_df))
    
}
