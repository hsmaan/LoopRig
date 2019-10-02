library(LoopRig)

setwd("../lncrna_loopnet_looptest/lncrna_loopnet_looptest/data/")

setwd("yue_looplists/")

### LoopsToRanges Test ###

loop_ranges <- LoopsToRanges("Rubin_2017.Epidermal_Keratinocyte_day0.hg19.peakachu-merged.loops", custom_cols = 5, loop_names = "test")

loop_ranges

