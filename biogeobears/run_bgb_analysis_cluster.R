## run BioGeoBEARS and BSMs and plot binned cladogenetic/anagenetic events and rates for tree set

bgb_dir="~/biogeobears"
tree_dir=file.path(bgb_dir,"trees")
out_dir=file.path(bgb_dir,"biogeobears_output")

prefix= "sortadate50_ASL"
#prefix= "sortadate50_merged"

bgbout_dir=file.path(out_dir,prefix) ## manually create directory and copy/create necessary input files
print(bgbout_dir)

if(prefix=="sortadate50_ASL") geogfn = file.path(bgbout_dir,"band_summary_expert_reduced.txt")
else if(prefix=="sortadate50_merged") geogfn = file.path(bgbout_dir,"band_summary_expert_merged.txt")


PREP_TREES=FALSE
BGB=FALSE
BSM=FALSE
BIN=TRUE
PLOTONLY=FALSE

###########################################
## randomly select and prepare 100 trees ##
###########################################

if(PREP_TREES){
  print("prepping posterior trees")
  source("~/biogeobears/scripts/prep_trees.R")
  
  ## random tree selection (only do once)
  tree_dir=file.path(hyb_dir,"time_calibration/StarBeast2/sortadate50/rank50_strictEstimate/phase3/run2/")
  treefile=file.path(tree_dir,"rank50_strictEstimate.species.trees")
  trees100file=sample_trees(treefile=treefile,burnin=.55,final.number=100,format="nex",prefix="sortadate50",outdir=bgb_dir)
  
  ## subset and format trees and save separately
  trees100file=file.path(bgb_dir,"sortadate50_100trees.nex")
  outdir=file.path(bgb_dir,"trees")
  prep_trees_BioGeoBEARS(treefile=trees100file,geogfn=geogfn,prefix="sortadate50",treeoutdir=outdir)
}

################################################
## run BioGeoBEARS and 100 BSMs for 100 trees ##
################################################
if(BGB){
  source("~/biogeobears/scripts/biogeobears_functions.R")

  print("running BioGeoBEARS")
  print(Sys.time())

  ## parallel function version
  run_bgb_parallel(batch_size=20,treedir=tree_dir,geogfn=geogfn, bgboutdir=bgbout_dir,BSM=BSM)
  print(Sys.time())

#  ## parallel script version (untested)
#  batch_size=20
#  source("~/biogeobears/scripts/run_bgb_parallel.R")
#  ## non-parallel version
#  run_BioGeoBEARS_plusBSM_multitrees(treedir=tree_dir,geogfn=geogfn, bgboutdir=bgbout_dir,BSM=BSM)
}

#########################################################################
## Bin and plot anagenetic/cladogenetic events by area across all runs ##
#########################################################################
if(BIN){
  source("~/biogeobears/scripts/process_BSM_events.R")
  print("binning/plotting")
  if(prefix=="sortadate50_ASL"){
    ## run with default A/S/L
    run_bin_plot_BSM_multi(bgbout_dir,prefix,plotonly=PLOTONLY)
  } else if(prefix=="sortadate50_merged"){
    run_bin_plot_BSM_multi(bgbout_dir,prefix=prefix,areas=c('L','M'),
                       all.codes=c('O','M','L','ML'),bin.width=0.5,colors=c("red","#288C8C"))
  } else if(prefix=="sortadate50_merged1Myr"){
    ## might need to reconfigure prefix system or copy input files for this one (same input, different binning)
    bgboutdir=file.path(bgb_dir,"sortadate50_time_expertBands_merged1Myr")
    dir.create(bgboutdir)
    run_bin_plot_BSM_multi(bgboutdir,prefix="sortadate50_merged1Myr",areas=c('L','M'),
                       all.codes=c('O','M','L','ML'),bin.width=1,colors=c("red","#288C8C"))

  } else if(prefix=="sortadate50_mergedNS"){
    run_bin_plot_BSM_multi(bgbout_dir,prefix="sortadate50_mergedNS",areas=c('A','B','C','D'),
                       all.codes=c("O","A","B","C","D","AB","AC","AD","BC","BD","CD",  
                               "ABC","ABD","ACD","BCD","ABCD"),bin.width=1,colors=c('#0000FF','#00CD00','#FFFF00','#FF0000'))

  }
}



