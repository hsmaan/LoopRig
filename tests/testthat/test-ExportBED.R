library(LoopRig)

# Load files

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)

# Get LoopRanges and ElementRanges objects

element_ranges <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1, custom_mcols = 4)

element_ranges_no_mcol <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1)

loops <- LoopsToRanges(ovary_loops, pancreas_loops, spleen_loops, custom_cols = 0)

consensus_loops <- ConsensusLoops(loops)

# Error handling

test_that("error handling", {
  
  expect_error(ExportBED(ovary_loops, index = 2, file_name = "inst/extdata/elements/exp_bed_dummy.bed"), "Incorrect object type. Please enter either a 'LoopRanges' or 'ElementRanges' object")
  
  # expect_error(ExportBED(consensus_loops, index = 1, file_name = "install/extdata/loops/exp_bed_dummy.bed"), "cannot open the connection")
  
  expect_error(ExportBED(consensus_loops, index = 1, mcol = TRUE, file_name = "inst/extdata/loops/exp_bed.bed"), "Object has no mcols, please choose mcol = FALSE")
  
  expect_error(ExportBED(element_ranges_no_mcol, index = 1, mcol = TRUE, file_name = "inst/extdata/elements/exp_bed_dummy.bed"), "Object has no mcols, please choose mcol = FALSE")
  
  # expect_error(ExportBED(element_ranges, index = 1, file_name = "inst/extdata/elements/exp_bed_dummy.bed"), "File already exists, cannot overwrite")
  
})

# Test output

test_that("expected output", {
  
  
  
})