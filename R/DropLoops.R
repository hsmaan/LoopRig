#' Drop loops from \emph{LoopRanges} objects using anchor and loop sizes
#'
#' Subset loops based on loop and anchor size filters. Can be used to filter the loops in the \emph{LoopRanges} object before or after calling ConsensusLoops
#' @param loop_ranges A \emph{LoopRanges} class object 
#' @param type A string indicating the type of filtering when determining which loops to drop:
#' \itemize{
#'  \item "anchor_size" - Drops loops based on given anchor sizes 
#'  \item "loop_size" - Drops loops based on end-to-end loop size
#'}
#' @param size A numerical vector indicating size range to keep (e.g. c(start, end)) 
#' @return A subsetted \emph{LoopRanges} class object 
#' @examples 
#' # Load loops into LoopRanges object 
#' ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", 
#' package = "LoopRig", mustWork = TRUE) 
#' loops_ovary <- LoopsToRanges(ovary_loops, custom_cols = 0)
#' 
#' # Subset loops based on total length between 100 to 100000 bp 
#' DropLoops(loops_ovary, type = "loop_size", size = c(100, 100000))
#' 
#' # Subset loops based on anchor size between 1000 to 25000 bp 
#' DropLoops(loops_ovary, type = "anchor_size", size = c(1000, 25000))
#' @import GenomicRanges
#' @import IRanges
#' @export  

DropLoops <- function(loop_ranges, type = NULL, size = NULL) {
  
  if (class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the 'loop_ranges' parameter")
  }
  
  if (length(size) != 2) {
    stop("Please enter a numerical vector of two integers for the 'size' parameter")
  }

  if (type == "anchor_size") {
    
    anchor_dropper <- function(loop_range, anchors) {
      
      anchor_1_hits <- which(IRanges::width(loop_range[[1]]) %in% seq(anchors[1], anchors[2]))
      anchor_2_hits <- which(IRanges::width(loop_range[[2]]) %in% seq(anchors[1], anchors[2])) 
      anchor_intersect <- intersect(anchor_1_hits, anchor_2_hits)
      
      anchor_1_subset <- loop_range[[1]][anchor_intersect]
      anchor_2_subset <- loop_range[[2]][anchor_intersect]
      
      subset_loop <- GRangesList("Anchor 1" = anchor_1_subset, "Anchor 2" = anchor_2_subset)
    }
    
    loop_ranges_sub <- lapply(loop_ranges, function(x) anchor_dropper(loop_range = x, anchors = size))
    return(structure(loop_ranges_sub, class = "LoopRanges"))
  }
  
  else if (type == "loop_size") {
    
    loop_dropper <- function(loop_range, loops) {
      
      loop_subset <- which((IRanges::end(loop_range[[2]]) - IRanges::start(loop_range[[1]])) %in% seq(loops[1], loops[2]))
      anchor_1_subset <- loop_range[[1]][loop_subset]
      anchor_2_subset <- loop_range[[2]][loop_subset]
      subset_loop <- GRangesList("Anchor 1" = anchor_1_subset, "Anchor 2" = anchor_2_subset)
    }
    
    loop_ranges_sub <- lapply(loop_ranges, function(x) loop_dropper(loop_range = x, loops = size))
    return(structure(loop_ranges_sub, class = "LoopRanges"))
  }
  
  else {
    stop("Please enter either 'loop_size' or 'anchor_size' for the 'type' parameter")
  }
}