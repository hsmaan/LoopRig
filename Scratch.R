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

element_ranges <- ElementsToRanges(enhancers, element_names = "enhancers", custom_cols = 8, custom_mcols = 10)

element_ranges

### DropLoops Tests ###

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(300000, 1000000))

drop_loop_ranges

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(0, 0))

drop_loop_ranges

### ConsensusLoops Tests ###

consensus_loops <- ConsensusLoops(loop_ranges)

consensus_loops


