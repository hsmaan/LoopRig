library(LoopRig)

# Load files

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

# Error handling 

test_that("error handling", {
  
  expect_error(LoopsToRanges(loop_names = NULL, custom_cols = 0, custom_mcols = NULL), "Please enter at least one BEDPE looping data file")
  
  expect_error(LoopsToRanges(ovary_loops, ovary_loops), "Duplicate BEDPE data-files entered")
  
  expect_error(LoopsToRanges(ovary_loops, pancreas_loops, custom_cols = 20), "Incorrect number of custom columns, please check your files again")
  
  expect_error(LoopsToRanges(c("100", "200", "300"), ovary_loops), "Error in reading BEDPE files, please ensure file specification, format, and location accuracy")
  
  expect_error(LoopsToRanges(pancreas_loops, ovary_loops, custom_mcols = 20), "Incorrect custom mcols location, please check your files again")
  
  ### Include chromosome coordinate check ### 
  
})

# Expected class output

test_that("class output", {
  
  expect_is(LoopsToRanges(ovary_loops, spleen_loops, custom_cols = 0), "LoopRanges")
  
})

