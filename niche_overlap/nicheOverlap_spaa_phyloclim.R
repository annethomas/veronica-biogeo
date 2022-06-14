### Calculate niche overlap from TTR predicted extent (reflecting "fundamental niche") with two packages, phyloclim and spaa

##### Notes:
### phyloclim requires SpatialGrid input and outputs a "niolap" matrix with Schoener's D and Hellinger's I; the format is necessary for other phyloclim functions but not flexible with others
### spaa requires a matrix from the raster cell pres/abs values, can choose from a range of metrics including D, and outputs a dist object that can easily be converted to a matrix (but this can't be used for phyloclim downstream)
### Schoener's D is identical for both; it's just the format and the other metric options that differ
### TTR prediction files are from scripts run on the cluster (08-predict_species_dist.R)

library(raster)
library(ape)
library(rgdal)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/"
dist_dir=file.path(data_dir,"distribution_data")
bgb_dir=file.path(data_dir,"biogeobears")
TTR_dir=file.path(data_dir,"TTR")
niche_dir=file.path(data_dir,"niche_overlap")

phy=read.tree(file.path(bgb_dir,"sortadate50_TTRtips.tre"))

## load functions (prep_species_raster, plot_species_raster)
source(file.path(script_dir,"niche_overlap/overlap_util_functions.R"))

#########################
## prep common to both ##
#########################
## load regional environmental grid for the predictions to create background raster
extent_fn=file.path(TTR_dir,"region_extent_NZ.Rdata")
load(extent_fn)
r <- raster(region.extent)

resolution=.01 ## default
#resolution=.05

res(r)<-resolution

if(resolution==.01){
  pred_dir=file.path(TTR_dir,"predict_grids_current")
  load(file.path(TTR_dir,"env_grid_current_NZ.Rdata"))
} else if(resolution==.05){
  pred_dir=file.path(TTR_dir,"predict_grids_current05")
  load(file.path(TTR_dir,"env_grid_current_NZ_05.Rdata"))
}


predfiles=list.files(pred_dir)
output_suffix = ".Rdata"

#####################
##### phyloclim #####
#####################
library(phyloclim)

spatial_grids=list()

for(i in 1:length(predfiles)) {
  print(i)
  #load the presence/absence grid
  f <-predfiles[i]
  species = substr(f,nchar("predict_grd_ft__")+1,nchar(f)-nchar(output_suffix))
  print(species)
  
  pred.list=readRDS(file.path(pred_dir,predfiles[i]))
  pred.grd = pred.list[[3]]
  print(head(pred.grd))
  
  pred.r=prep_species_raster(r,env.grid,pred.grd) 
  #plot_species_raster(pred.r,title=species)
  
  pred.spgr = as(pred.r,"SpatialGridDataFrame")
  
  spatial_grids[[species]] = pred.spgr
  
}

TTR.no.full = phyloclim::niche.overlap(spatial_grids) #The upper triangle contains pairwise comparisons of niche overlap in terms of D, whereas the lower triangle contains values of I.
TTR.no.phy = clean_matrix(TTR.no.full,phy$tip.label)


##########
## spaa ##
##########
library(spaa)
## get length of individual pred.grd
pred.list=readRDS(file.path(pred_dir,"predict_grd_ft__Veronica_amplexicaulis.Rdata"))
pred.grd = pred.list[[3]]
mat.pa=matrix(nrow=length(pred.grd),ncol=length(predfiles))
SPNAMES=c()

for(i in 1:length(predfiles)) {
  
  print(i)
  f <-predfiles[i]
  species = substr(f,nchar("predict_grd_ft__")+1,nchar(f)-nchar(output_suffix))
  print(species)
  
  pred.list=readRDS(file.path(pred_dir,predfiles[i]))
  pred.grd = pred.list[[3]]
  
  #put the data in the matrices, save species names for column labels later
  mat.pa[,i]<-pred.grd
  SPNAMES<-c(SPNAMES,species)
  
}

sp.site=mat.pa
colnames(sp.site)<-SPNAMES 

schoener.predicted<-spaa::niche.overlap((sp.site),method="schoener")
sch.mat.full=as.matrix(schoener.predicted)

sch.mat.phy = clean_matrix(sch.mat.full,phy$tip.label)

###############################
## compare phyloclim and spaa #
###############################
plot(sch.mat[upper.tri(sch.mat)],TTR.no.full[upper.tri(TTR.no.full)])

## save phyloclim niolap objects and spaa Schoener's D (default to 0.01 degree resolution)
saveRDS(list(TTR.no.full=TTR.no.full,TTR.no.phy=TTR.no.phy,sch.mat.full=sch.mat.full,sch.mat.phy=sch.mat.phy),
        file=file.path(niche_dir,"TTRoverlap_phyloclim_spaa_01.Rdata"))
## option if both resolutions were run (edit script to reflect naming)
saveRDS(list(TTR.no.full01=TTR.no.full01,TTR.no.full05=TTR.no.full05,sch.mat01=sch.mat01,sch.mat05=sch.mat05),
        file=file.path(niche_dir,"TTRoverlap_phyloclim_spaa_01_05.Rdata"))

# TTR.no.list=readRDS(file.path(niche_dir,"TTRoverlap_phyloclim_spaa_01_05.Rdata"))

####################################################
## run/save other spaa metrics (incl schoener's D) #
####################################################
print("pianka")
pianka.predicted<-niche.overlap((sp.site),method="pianka")
print("levins")
levins.predicted<-niche.overlap((sp.site),method="levins")
print("morisita")
morisita.predicted<-niche.overlap((sp.site),method="morisita")

## save the spaa niche overlap calcs
saveRDS(list(schoener=schoener.predicted,pianka=pianka.predicted,levins=levins.predicted,morisita=morisita.predicted),
        file=file.path(TTR_dir,"niche_overlap_TTR_spaa.Rdata"))    ## was sdm_overlap_predicted_1.Rdata
niche_idxs=readRDS(file.path(TTR_dir,"niche_overlap_TTR_spaa.Rdata"))

