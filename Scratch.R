library(LoopRig)

### LoopsToRanges Tests ###

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/yue_looplists")

loop_ranges <- LoopsToRanges("Rubin_2017.Epidermal_Keratinocyte_day0.hg19.peakachu-merged.loops", "ENCODE3.LNCAP.hg19.peakachu-merged.loops", custom_cols = 0, loop_names = c("test1", "test2"))

loop_ranges

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/rao_looplists/")

loop_ranges <- LoopsToRanges("GSE63525_HeLa_HiCCUPS_looplist.txt", custom_cols = 14, loop_names = "test")

### ElementsToRanges Tests ### 

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/elements")

element_ranges <- ElementsToRanges("enhancers.bed", element_names = "enhancers", custom_cols = 8, custom_mcols = 10)

element_ranges

### DropLoops Tests ###

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(300000, 1000000))

drop_loop_ranges

drop_loop_ranges <- DropLoops(loop_ranges, "loop_size", c(0, 0))

drop_loop_ranges

### ConsensusLoops Tests ###

setwd("~/Documents/lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/yue_looplists")

loop_ranges <- LoopsToRanges("Rubin_2017.Epidermal_Keratinocyte_day0.hg19.peakachu-merged.loops", "Dixon_2015.H1-TRO.hg19.peakachu-merged.loops", "Schmitt_2016.Ovary.hg19.peakachu-merged.loops", "Schmitt_2016.Ovary.hg19.peakachu-merged.loops", custom_cols = 0, loop_names = c("test1", "test2", "test3"))

consensus_loops <- ConsensusLoops(loop_ranges)



