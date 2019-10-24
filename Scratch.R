library(LoopRig)

### LoopsToRanges Tests ###

ovary_loops <- system.file("extdata/loops", "ovary_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
pancreas_loops <- system.file("extdata/loops", "pancreas_hg19.bedpe", package = "LoopRig", mustWork = TRUE)
spleen_loops <- system.file("extdata/loops", "spleen_hg19.bedpe", package = "LoopRig", mustWork = TRUE)

loop_ranges <- LoopsToRanges(ovary_loops, spleen_loops, custom_cols = 0, loop_names = c("test1", "test2"))

loop_ranges

loop_ranges <- LoopsToRanges(ovary_loops, ovary_loops, custom_cols = 14, loop_names = "test")

### ElementsToRanges Tests ### 

enhancers <- system.file("extdata/elements", "enhancers.bed", package = "LoopRig", mustWork = TRUE)
lncrnas <- system.file("extdata/elements", "lncrnas.bed", package = "LoopRig", mustWork = TRUE)

element_ranges <- ElementsToRanges(enhancers, lncrnas, element_names = c("enhancers", "lncrnas"), custom_cols = 8, custom_mcols = 4)

element_ranges

### DropLoops Tests ###

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(300000, 1000000))

drop_loop_ranges

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(0, 0))

drop_loop_ranges

### ConsensusLoops Tests ###

loops <- LoopsToRanges(ovary_loops, spleen_loops, pancreas_loops)

consensus_loops <- ConsensusLoops(loop_ranges, keep_all = TRUE)

consensus_loops

### LinkedElements Tests ###

linked_elements <- LinkedElements(loop_ranges, element_ranges[[1]], element_ranges[[2]])

linked_elements <- LinkedElements(loop_ranges[[1]], element_ranges[[1]], element_ranges[[2]])

linked_elements <- LinkedElements(consensus_loops, element_ranges[[1]], element_ranges[[2]])

linked_elements

View(linked_elements)
