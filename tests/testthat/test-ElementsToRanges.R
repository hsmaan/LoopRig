library(LoopRig)

# Load files

enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)

# Error handling

test_that("error handling", {
  
  expect_error(ElementsToRanges(element_names = NULL), "Please enter at least one BED4/12 element data file")
  
  expect_error(ElementsToRanges(enhancers, enhancers), "Duplicate BED files entered, only unique entries allowed")
  
  expect_error(ElementsToRanges(c(2,4,5,6)), "Error in reading BED files, please ensure file specification and format accuracy")
  
  expect_error(ElementsToRanges(enhancers, promoters), "Please enter an appropriate value for the custom_cols parameter")
  
  expect_error(ElementsToRanges(enhancers, promoters, custom_cols = 25), "Incorrect number of custom columns, please check your files again")
  
  expect_error(ElementsToRanges(enhancers, promoters, custom_cols = 1, custom_mcols = 22), "Incorrect custom mcols location, please check your files again")
  
  ### Include chromosome coordinate check ### 
  
})

test_that("class output", {
  
  expect_is(ElementsToRanges(enhancers, promoters, custom_cols = 1, custom_mcols = 4), "ElementRanges")
  
  expect_is(ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1, custom_mcols = 4), "ElementRanges")
  
})