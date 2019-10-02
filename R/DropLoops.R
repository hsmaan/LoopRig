#' A variety of options to drop loops from the 'LoopRanges' object
#'
#' Drops loops based on a type of filter and given options. Can be used to filter the loops in the 'LoopRanges' object before or after determining consensus.
#' @param type Type of filtering when determining which loops to drop:
#' \itemize{
#'  \item "Anchor Size" - Drops loops based on given anchor sizes as the option.
#'  \item 
#'}
#' @param ... Options for filtering based on method type selected:
#' \itemize{
#'  \item For method "Anchor Size" - a single integer or numerical vector of anchor sizes. 
#'  \item For method "Loop Size" 
#'}
#' @export  

DropLoops <- function(type, ...) {
  
  f_options <- list(...)

  if (type == "Anchor Size") {
    
    anchor_dropper <- function(loop_range, anchors) {
      
      '%ni%' <- Negate('%in%')
      
      anchor_1_subset <- loop_range[[1]][width(loop_range[[1]]) %ni% (anchors+1)] 
      anchor_2_subset <- loop_range[[2]][width(loop_range[[2]]) %ni% (anchors+1)] 
      subset_loop <- GRangesList(anchor_1_subset, anchor_2_subset)
      names(subset_loop) <- c("Anchor 1", "Anchor 2")
      subset_loop
      
    }
    
    loop_ranges <- lapply(loop_ranges, function(x) anchor_dropper(loop_range = x, anchors = f_options[[1]]))
    
  }
  
  else if (type == "Loop Size") {
    
    
    
  }
}