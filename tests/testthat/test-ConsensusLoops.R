library(LoopRig)

# Load files

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

# Create LoopRanges objects

loops <- LoopsToRanges(ovary_loops, pancreas_loops, spleen_loops, custom_cols = 0)

loops_single <- LoopsToRanges(ovary_loops, custom_cols = 0)

# Error handling

test_that("error handling", {
  
  expect_error(ConsensusLoops(loop_ranges = c("a","b","c")), "Please enter an object of class 'LoopRanges' for the loop_ranges parameter")
  
  expect_error(ConsensusLoops(loop_ranges = loops_single), "Please enter a LoopRanges object with at least 2 loopsets")
  
  expect_error(ConsensusLoops(loop_ranges = loops, split_anchors = TRUE), "Please enter the loop anchor resolutions in numerical vector format")
  
  expect_error(ConsensusLoops(loop_ranges = loops, split_anchors = TRUE, resolutions = c(200, 2000), overlap_threshold = 22), "For split anchors, please enter a percentage overlap threshold between 0 and 1")
  
  expect_error(ConsensusLoops(loop_ranges = loops, split_anchors = "ABC"), "Please enter either 'TRUE' or 'FALSE' for the split_anchors parameter")
  
  
  
})

