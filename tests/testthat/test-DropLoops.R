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
  
  expect_error(DropLoops(ovary_loops, type = "loop_size", size = c(1000, 10000)), "Please enter an object of LoopRanges class for the 'loop_ranges' parameter")
  
  expect_error(DropLoops(loops, type = "type", size = c(1000, 20000)), "Please enter either 'loop_size' or 'anchor_size' for the 'type' parameter")
  
  expect_error(DropLoops(loops, type = "loop_size", size = 1000), "Please enter a numerical vector of two integers for the 'size' parameter")
  
})

# Expected class output

test_that("class output", {
  
  expect_is(DropLoops(loops, type = "loop_size", size = c(100, 10000)), "LoopRanges")
  
  expect_is(DropLoops(loops, type = "anchor_size", size = c(1000, 25000)), "LoopRanges")
  
})