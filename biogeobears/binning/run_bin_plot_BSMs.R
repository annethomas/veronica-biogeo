bgb_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/biogeobears/"
source("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography/biogeobears/bin_plot_BSMs.R")

## run with default A/S/L
run_dir=file.path(bgb_dir,"sortadate50_time_expertBands/")
run_dir=file.path(bgb_dir,"test/test_incl_lineage")
run_bin_plot_BSM_multi(run_dir,"test_sortadate50",plotonly=TRUE)

run_dir=file.path(bgb_dir,"sortadate50_time_expertBands_merged")
run_bin_plot_BSM_multi(run_dir,prefix="sortadate50_merged",areas=c('L','M'),
                       all.codes=c('O','M','L','ML'),bin.width=0.5,colors=c('blue','green'))

run_dir=file.path(bgb_dir,"sortadate50_time_expertBands_merged1Myr")
dir.create(run_dir)
run_bin_plot_BSM_multi(run_dir,prefix="sortadate50_merged1Myr",areas=c('L','M'),
                       all.codes=c('O','M','L','ML'),bin.width=1,colors=c('blue','green'))

run_dir=file.path(bgb_dir,"sortadate50_time_expertBands_mergedNS2")
run_bin_plot_BSM_multi(run_dir,prefix="sortadate50_mergedNS",areas=c('A','B','C','D'),
                       all.codes=c("O","A","B","C","D","AB","AC","AD","BC","BD","CD",  
                               "ABC","ABD","ACD","BCD","ABCD"),bin.width=1,colors=c('#0000FF','#00CD00','#FFFF00','#FF0000'))

## with 100 trees
run_dir=file.path(bgb_dir,"biogeobears_output/sortadate50_merged")
run_bin_plot_BSM_multi(run_dir,prefix="sortadate50_merged",areas=c('L','M'),
                       all.codes=c('O','M','L','ML'),bin.width=0.5,colors=c('blue','green'))


