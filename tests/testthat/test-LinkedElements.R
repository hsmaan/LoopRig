library(LoopRig)

# Load files

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)

# Get LoopRanges and ElementRanges objects

element_ranges <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1, custom_mcols = 4)

loops <- LoopsToRanges(ovary_loops, pancreas_loops, spleen_loops, custom_cols = 0)

consensus_loops <- ConsensusLoops(loops)

# Error handling

test_that("error handling", {
  
  expect_error(LinkedElements(ovary_loops, element_ranges[[1]], element_ranges[[2]]), "Please enter an object of LoopRanges class for the loop_ranges parameter")
  
  expect_error(LinkedElements(loops, element_ranges[[1]], element_ranges[[2]]), "Please enter a conseus LoopRanges object with only one range for the loop_ranges parameter")
  
  expect_error(LinkedElements(consensus_loops, element_ranges[[1]], elements_ranges[[2]], range_out_x = TRUE, range_out_y = TRUE), "Can only output either element_ranges_x or element_ranges_y as ElementRanges object")
  
})

# Expected class output 

test_that("class output", {
  
  expect_is(LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]]), "data.frame")
  
  expect_is(LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]], range_out_x = TRUE), "ElementRanges")
  
  expect_is(LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]], range_out_y = TRUE), "ElementRanges")
  
  expect_is(LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]], overlap_threshold = 10000000), "data.frame")
  
})