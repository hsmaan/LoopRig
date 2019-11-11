#' Subset an object of class LoopRanges using consensus options
#'
#' Performs filtering on looping data in LoopRanges objects based on custom parameters and returns a single GRangesList object indicating one looping dataset with two anchors
#' @param loop_ranges An object of 'LoopRanges' class created from the LoopsToRanges() function
#' @param stringency Integer (n>=0) indicating the number of looping datasets a loop from a given dataset must overlap with to be considered a consensus loop
#' @param overlap_threshold Single numerical input in either percentage (0<=n<=1) overlap format if split_anchors = TRUE, or in base pair number format (n>=0) in split_anchors=FALSE (default=1)
#' @param split_anchors A boolean (TRUE/FALSE) that determines if the different loop anchor sizes are considered together (default=TRUE) or seperately (FALSE)
#' @param resolutions An optional numerical vector of anchor sizes - to be used only when split_anchors=TRUE
#' @param keep_all If TRUE, keeps all of the loops (concatenation of looping datasets)
#' @return A 'LoopRanges' class object for the consensus loops 
#' @import GenomicRanges
#' @import IRanges
#' @importFrom S4Vectors queryHits subjectHits intersect.Vector pc 
#' @export

ConsensusLoops <- function(loop_ranges, stringency = 1, overlap_threshold = 1, split_anchors = FALSE, resolutions = NULL, keep_all = FALSE) {
  
  if (class(loop_ranges) != "LoopRanges") {
    stop("Please enter an object of class 'LoopRanges' for the loop_ranges parameter")
  }
  
  if (length(loop_ranges) < 2) {
    stop("Please enter a LoopRanges object with at least 2 loopsets")
  }
  
  if (keep_all == TRUE) {
    
    consensus_loops <- Reduce(S4Vectors::pc, loop_ranges)
    names(consensus_loops) <- c("Anchor 1", "Anchor 2")
    return(structure(list("Consensus" = consensus_loops), class = "LoopRanges"))
  }

  if (split_anchors == FALSE) {
    
    hit_counter_1 <- function(looprange1, looprange2, threshold = overlap_threshold, type = "other") {
      anchor_1_hits <- findOverlaps(looprange1[[1]], looprange2[[1]], minoverlap = threshold)
      anchor_2_hits <- findOverlaps(looprange1[[2]], looprange2[[2]], minoverlap = threshold)
      consensus_hits <- intersect.Vector(anchor_1_hits, anchor_2_hits)
      
      if (type == "self") {
        return(consensus_hits[(which(queryHits(consensus_hits) != subjectHits(consensus_hits)))])
      } 
  
      else if (type == "other") {
        consensus_hits2 <- consensus_hits[which(!duplicated(queryHits(consensus_hits)))]
        return(consensus_hits2)
      }
    }
    
    total_counts <- lapply(loop_ranges, function(x) lapply(loop_ranges, function(y) hit_counter_1(x, y)))
    for (i in 1:length(loop_ranges)) { 
      total_counts[[i]][[i]] <- NULL 
    }
    
    hit_subset <- function(count_list) {
      table(as.vector(unlist(sapply(count_list, queryHits))))
    }
  }
  
  else if (split_anchors == TRUE) {
    
    if (is.null(resolutions) == TRUE) {
      stop("Please enter the loop anchor resolutions in numerical vector format")
    }
    
    if ((overlap_threshold < 0) || (overlap_threshold > 1)) {
      stop("For split anchors, please enter a percentage overlap threshold between 0 and 1")
    }
    
    resolutions <- resolutions + 1
    
    hit_counter_2 <- function(loop_range1, loop_range2, threshold = overlap_threshold, type = "other", divisions = resolutions) {
      
      overlapper <- function(looprange1 = loop_range1, looprange2 = loop_range2, anchor, threshold = overlap_threshold, division) {
        range1_split <- looprange1[[anchor]][width(looprange1[[anchor]]) == division]
        range2_split <- looprange2[[anchor]][width(looprange2[[anchor]]) == division]
        overlaps <- as.matrix(findOverlaps(range1_split, range2_split, minoverlap = division*threshold))
      }
      
      anchor_1_hits <- base::Reduce(base::rbind, mapply(overlapper,  division = divisions, anchor = 1))
      anchor_2_hits <- base::Reduce(base::rbind, mapply(overlapper,  division = divisions, anchor = 2))
      
      consensus_hits <- merge(anchor_1_hits, anchor_2_hits)
      
      if (type == "self") {
        return(consensus_hits[(which(queryHits(consensus_hits) != subjectHits(consensus_hits)))])
      } 
      else if (type == "other") {
        return(unique(consensus_hits[[1]]))
      }
      
    }
    
    total_counts <- lapply(loop_ranges, function(x) lapply(loop_ranges, function(y) hit_counter_2(x, y)))
    for (i in 1:length(loop_ranges)) { 
      total_counts[[i]][[i]] <- NULL 
    }
    
    hit_subset <- function(count_list) {
      table(as.vector(unlist(count_list)))
    }
  } 
  
  else {
    stop("Please enter either 'TRUE' or 'FALSE' for the split_anchors parameter")
  }
  
  count_table <- lapply(total_counts, hit_subset)
  
  subset_ranges <- function(loop_range, count_list) {
    anchor1 <- loop_range[[1]][as.numeric(names(count_list[count_list >= (stringency)]))]
    anchor2 <- loop_range[[2]][as.numeric(names(count_list[count_list >= (stringency)]))]
    return(GRangesList(anchor1, anchor2))
  }
  
  subset_loops <- mapply(subset_ranges, loop_range = loop_ranges, count_list = count_table)
  combined_subset <- base::Reduce(pc, subset_loops)
  dup_indices <- which(duplicated(combined_subset[[1]]) & duplicated(combined_subset[[2]]))
  anchor_1_consensus <- combined_subset[[1]][dup_indices]
  anchor_2_consensus <- combined_subset[[2]][dup_indices]
  consensus_loops <- GRangesList(anchor_1_consensus, anchor_2_consensus)
  names(consensus_loops) <- c("Anchor 1", "Anchor 2")
  return(structure(list("Consensus" = consensus_loops), class = "LoopRanges"))
  
}
