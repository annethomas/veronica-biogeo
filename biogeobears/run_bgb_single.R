#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS")
library(BioGeoBEARS)
library(ape)

## note: before running this, check that there are no graphical devices open as the pdfs often glitch when saving.

## paths
data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
BGB_dir=file.path(data_dir,"biogeobears") # parent bgb dir
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/biogeobears/"
source(file.path(script_dir,"biogeobears_functions.R"))



## basic setup
genus = "Veronica"
trfn = file.path(BGB_dir,"sortadate50_tree.tre")
trfn_TTR = file.path(BGB_dir,"sortadate50_TTRtips.tre")
trfn_V = file.path(BGB_dir,"sortadate50_tree_V.tre") ## labelled with V_ instead of Veronica_

##
prefix="sortadate50_merged"
prefix="sortadate50_ASL"
prefix="sortadate50_merged_TTRtips"
prefix="sortadate50_ASL_TTRtips"

## if first time: copy time_periods.txt, manual_dispersal_multipliers.txt and areas_allowed.txt to new run_dir

if(prefix=="sortadate50_merged"){
  ## two bands (merged mountain band) with all BGB tips
  run_dir=file.path(BGB_dir,"sortadate50_time_expertBands_merged")
  geogfn = file.path(run_dir,"band_summary_expert_merged.txt")
  trfn=trfn
  time_strat=TRUE
} else if(prefix=="sortadate50_ASL"){
  ## three bands with all BGB tips
  run_dir=file.path(BGB_dir,"sortadate50_time_expertBands")
  geogfn = file.path(run_dir,"band_summary_expert_reduced.txt")
  trfn=trfn
  time_strat=TRUE
} else if(prefix=="sortadate50_merged_TTRtips"){
  ## merged areas with TTR tips
  run_dir=file.path(BGB_dir,"sortadate50_merged_TTRtips")
  geogfn = file.path(run_dir,"band_summary_expert_merged.txt")
  trfn=trfn_TTR
  time_strat=TRUE
} else if(prefix=="sortadate50_ASL_TTRtips"){
  ## merged areas with TTR tips
  run_dir=file.path(BGB_dir,"sortadate50_merged_TTRtips")
  geogfn = file.path(run_dir,"band_summary_expert_reduced.txt")
  trfn=trfn_TTR
  time_strat=TRUE
}


#runs the models
run_biogeobears(trfn=trfn,geogfn=geogfn,run_dir=run_dir,name=genus,max_range_size = 6,time_strat=time_strat)
#source("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography/biogeobears/BioGeoBEARS_generic.R")
restable = read.table(file.path(run_dir,"Veronica_restable_AICc_rellike.txt"))


## see deprecated/biogeobears_prep_stream.R for infomap5 


