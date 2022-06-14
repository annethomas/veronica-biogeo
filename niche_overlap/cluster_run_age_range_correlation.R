#####################
##### phyloclim #####
#####################
library(raster)
library(ape)
library(rgdal)
library(phyloclim)
source("~/TTR/niche_overlap/scripts/age_range_correlation_warren.R")

set.seed(123) #could also make an commandline input

TTR_dir="~/TTR"
data_dir="/production/TTR_env_data"
pred_dir=file.path(data_dir,"predict/predict_grids_current")
resolution=0.01
input_file=commandArgs(trailingOnly=TRUE)[1]

input_list=readRDS(input_file)
phy=input_list$phy
ovlap=input_list$overlap

print("starting age-range correlation")
print(Sys.time())
arc.list = age.range.correlation.SES(phy,ovlap,n=100)
print(Sys.time())


saveRDS(arc.list,file=file.path(TTR_dir, "niche_overlap",paste0("arc.seed123_",basename(input_file))))

