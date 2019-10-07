#' A variety of options to drop loops from \emph{LoopRanges} objects
#'
#' Drops loops based on a type of filter and given options. Can be used to filter the loops in the \emph{LoopRanges} object before or after determining consensus
#' @param loop_ranges A \emph{LoopRanges} class object 
#' @param type A string indicating the type of filtering when determining which loops to drop:
#' \itemize{
#'  \item "Anchor" - Drops loops based on given anchor sizes 
#'  \item "Loop" - Drops loops based on end-to-end loop size
#'}
#' @param size A numerical vector indicating size range to keep (e.g. c(start, end)) 
#' @return A subsetted \emph{LoopRanges} class object 
#' @export  

DropLoops <- function(loop_ranges, type, size) {
  
  if (class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of LoopRanges class for the 'loop_ranges' parameter")
  }

  if (type == "Anchor Size") {
    
    anchor_dropper <- function(loop_range, anchors) {
      
      anchor_1_subset <- loop_range[[1]][width(loop_range[[1]]) %in% anchors] 
      anchor_2_subset <- loop_range[[2]][width(loop_range[[2]]) %in% anchors] 
      subset_loop <- GRangesList(anchor_1_subset, anchor_2_subset)
      names(subset_loop) <- c("Anchor 1", "Anchor 2")
      subset_loop
    }
    
    loop_ranges_sub <- lapply(loop_ranges, function(x) anchor_dropper(loop_range = x, anchors = size))
    return(structure(loop_ranges_sub, class = "LoopRanges"))
  }
  
  else if (type == "Loop Size") {
    
    loop_dropper <- function(loop_range, loops) {
      
      loop_subset <- which(start((loop_range[[1]]) - end(loop_range[[2]])) %in% loops)
      anchor_1_subset <- loop_range[[1]][loop_subset]
      anchor_2_subset <- loop_range[[2]][loop_subset]
      subset_loop <- GRangesList(anchor_1_subset, anchor_2_subset)
      names(subset_loop) <- c("Anchor 1", "Anchor 2")
      subset_loop
    }
    
    loop_ranges_sub <- lapply(loop_ranges, function(x) loop_dropper(loop_range = x, loops = size))
    return(structure(loop_ranges_sub, class = "LoopRanges"))
  }
}