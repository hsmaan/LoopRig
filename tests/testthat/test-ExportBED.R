library(LoopRig)

# Load files

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
ovary_loops_mcols <- system.file("extdata/loops", "ovary_hg19_mcol.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)


enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
promoters <- system.file("extdata/elements", "promoters.bed", package = "LoopRig", mustWork = TRUE)

loops_dir <- system.file("extdata/loops", package = "LoopRig", mustWork = TRUE)
elements_dir <- system.file("extdata/elements", package = "LoopRig", mustWork = TRUE)

  
# Get LoopRanges and ElementRanges objects

element_ranges <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1, custom_mcols = 4)

element_ranges_no_mcol <- ElementsToRanges(enhancers, promoters, element_names = c("enhancers", "promoters"), custom_cols = 1)

loops <- LoopsToRanges(ovary_loops, pancreas_loops, spleen_loops, custom_cols = 0)

loops_mcols <- LoopsToRanges(ovary_loops_mcols, custom_cols = 1, custom_mcols = 7, loop_names = "ovary_mcols")

consensus_loops <- ConsensusLoops(loops)

# Error handling

test_that("error handling", {
  
  expect_error(ExportBED(ovary_loops, index = -1, file_name = paste(loops_dir, "exp_bed_dummy2.bed", sep = "/")), "Please enter a positive integer indicating the index of LoopRanges or ElementRanges object to export")
  
  expect_error(ExportBED(ovary_loops, index = 2, file_name = paste(loops_dir, "exp_bed_dummy2.bed", sep = "/")), "Incorrect object type. Please enter either a 'LoopRanges' or 'ElementRanges' object")
  
  expect_error(ExportBED(consensus_loops, index = 1, mcol = TRUE, file_name = paste(loops_dir, "exp_bed_dummy2.bed", sep = "/")), "Object has no mcols, please choose mcol = FALSE")
  
  expect_error(ExportBED(element_ranges_no_mcol, index = 1, mcol = TRUE, file_name =  paste(elements_dir, "exp_bed_dummy.bed2", sep = "/")), "Object has no mcols, please choose mcol = FALSE")
  
  expect_error(ExportBED(element_ranges, index = 1, file_name = paste(elements_dir, "exp_bed_dummy.bed", sep = "/")), "File already exists, cannot overwrite")
  
})

# Test output

test_that("expected output", {
  
  expect_equal(ExportBED(element_ranges, index = 1, file_name = paste(elements_dir, "exp_bed_dummy2.bed", sep = "/")), paste("BED file exported to ", paste(elements_dir, "exp_bed_dummy2.bed", sep = "/"), sep = ""))
  
  expect_equal(ExportBED(element_ranges, index = 1, mcol = TRUE, file_name = paste(elements_dir, "exp_bed_dummy3.bed", sep = "/")), paste("BED file exported to ", paste(elements_dir, "exp_bed_dummy3.bed", sep = "/"), sep = ""))
  
  expect_equal(ExportBED(consensus_loops, index = 1, file_name = paste(loops_dir, "exp_bed_dummy2.bed", sep = "/")), paste("BEDPE file exported to ", paste(loops_dir, "exp_bed_dummy2.bed", sep = "/"), sep = ""))
  
  expect_equal(ExportBED(loops_mcols, index = 1, mcol = TRUE, file_name = paste(loops_dir, "exp_bed_dummy3.bed", sep = "/")), paste("BEDPE file exported to ", paste(loops_dir, "exp_bed_dummy3.bed", sep = "/"), sep = ""))
  
})

file.remove(paste(elements_dir, "exp_bed_dummy2.bed", sep = "/"))
file.remove(paste(elements_dir, "exp_bed_dummy3.bed", sep = "/"))
file.remove(paste(loops_dir, "exp_bed_dummy2.bed", sep = "/"))
file.remove(paste(loops_dir, "exp_bed_dummy3.bed", sep = "/"))