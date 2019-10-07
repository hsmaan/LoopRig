library(LoopRig)

### LoopsToRanges Tests ###

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/yue_looplists")

loop_ranges <- LoopsToRanges("Rubin_2017.Epidermal_Keratinocyte_day0.hg19.peakachu-merged.loops", custom_cols = 0, loop_names = "test")

loop_ranges

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/rao_looplists/")

loop_ranges <- LoopsToRanges("GSE63525_HeLa_HiCCUPS_looplist.txt", custom_cols = 14, loop_names = "test")

### ElementsToRanges Tests ### 

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/elements")

element_ranges <- ElementsToRanges("enhancers.bed", element_names = "enhancers", custom_cols = 8, custom_mcols = 10)

element_ranges

### DropLoops Tests ###

drop_loop_ranges <- DropLoops(loop_ranges, "Anchor Size", c(0, 10000))
