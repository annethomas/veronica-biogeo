library(rgdal)
library(raster)
library(rgeos)
library(sf)
library(rnaturalearth)
library(dplyr)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
dist_dir = file.path(data_dir,"distribution_data/")
TTR_dir=file.path(data_dir,"TTR")
niche_dir=file.path(data_dir,"biogeography_paper/niche_overlap")
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github"

## load functions (incl. clean_matrix)
source(file.path(script_dir,"niche_overlap/overlap_util_functions.R"))

###############################
### tree and occurrence data ##
###############################
phy=read.tree(file.path(TTR_dir,"sortadate50_TTRtips.tre"))
# occ=read.csv(file.path(dist_dir,"Veronica_eflora_strictClean_norepeats_20210226.csv"))
occ=read.csv(file.path(dist_dir,"Veronica_gbif_eflora_lucas_stricterClean.csv"))

####################
### spatial data ###
####################
wgs_proj4str = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
nztm_proj4str = "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
nztm.epsg=2193
## rnaturalearth
nz <- ne_countries(scale = 10, country = "New Zealand", returnclass = "sf") #might prompt some additional installations
nz=st_transform(nz[,-c(1:94)],crs=nztm.epsg) # remove fields
e=extent(162,179.99,-49.01080, -31.22159)
e.nztm=extent(796319.105,	2266528.473, 4512772.651, 6524698.951)
nz.sf=st_crop(nz,e.nztm)
plot(nz.sf)

## raster version (similar to tnc -> r.nz in other scripts, but different input NZ polygon)
nz.rast=raster(as_Spatial(nz.sf))
res(nz.rast) = 5000
nz.rast=rasterize(as_Spatial(nz.sf),nz.rast,1)
plot(nz.rast)

png_dir=file.path(dist_dir,"species_maps/eflora_buffers_wgbif")
#png_dir=file.path(dist_dir,"species_maps/eflora_buffers")
dir.create(png_dir)

species=unique(occ$species)
#species=species[c(2,13)] ## 2 species good for testing

#####################################
## create buffer polygons and save ##
## rasterize and save spatial grid ##
#####################################

plot=TRUE

buffer.polygons=list()
spatial_grids=list()

for(sp in species){
  print(sp)
  ## load and reproject points using sf
  sp.sfpoints=st_as_sf(occ[occ$species==sp,c("longitude","latitude")],coords=c("longitude","latitude"),crs=wgs_proj4str)
  #st_is_longlat(sp.sfpoints)
  sp.sfpoints=st_transform(sp.sfpoints,crs=nztm.epsg)
  
  ## buffer: 50km around points, subtract 40 km from outside perimeter
  sp.buffer=st_buffer(sp.sfpoints,dist=50000) %>% st_union() %>% st_buffer(dist=-40000)

  # intersect and crop with NZ boundaries
  nz.buffer=st_intersection(nz.sf,sp.buffer)
  
  buffer.polygons[[sp]]=nz.buffer
  
  if(plot){
    png(file.path(png_dir,paste0(gsub(" ","_",sp),".png")),3000,2018)
    plot(nz.sf)
    plot(nz.buffer,col="grey",add=TRUE)
    plot(sp.sfpoints,pch=16,cex=.1,col="red",add=TRUE)
    dev.off()
  }

  ##rasterize buffer
  nz.buff.rast=raster(as_Spatial(nz.buffer),ext=extent(nz.rast))
  res(nz.buff.rast) = 5000
  nz.buff.rast=rasterize(as_Spatial(nz.buffer),nz.buff.rast,1)

  ## sum layers to combine species range with NZ background
  nz.buff.sum=sum(nz.rast, nz.buff.rast,na.rm=TRUE)
  #plot(nz.buff.sum)
  nz.buff.sum@data@values[nz.buff.sum@data@values == 0] = NA
  nz.buff.sum=nz.buff.sum - 1
  range.grid = as(nz.buff.sum,"SpatialGridDataFrame")
  
  spatial_grids[[sp]] = range.grid
}

#####################
##### phyloclim #####
#####################
ro.full = phyloclim::niche.overlap(spatial_grids) #The upper triangle contains pairwise comparisons of niche overlap in terms of D, whereas the lower triangle contains values of I.
row.names(ro.full)=gsub(" ","_",row.names(ro.full))
colnames(ro.full)=gsub(" ","_",colnames(ro.full))
ro.phy = clean_matrix(ro.full,phy$tip.label)

## currently the preferred dataset for downstream analyses, see process_ARC_range_niche.R
#saveRDS(list(ro.full=ro.full,ro.phy=ro.phy),file=file.path(niche_dir,"range_overlap_phyloclim.Rdata"))
saveRDS(list(ro.full=ro.full,ro.phy=ro.phy),file=file.path(niche_dir,"range_overlap_phyloclim_wgbif.Rdata"))


########################
## jaccard similarity ##
########################

jacc.scores=matrix(nrow=length(species),ncol=length(species))
for(i in 1:length(species)){
  for(j in 1: length(species)){
    print(paste(species[i],species[j]))
    buffer_int=st_intersection(buffer.polygons[[i]],buffer.polygons[[j]])
    if(length(st_dimension(buffer_int))==0) {jacc=0; print("0 overlap")} else{
      
      buffer_union=st_union(buffer.polygons[[i]],buffer.polygons[[j]])
      int.area=st_area(buffer_int) 
      union.area=st_area(buffer_union)
      jacc=int.area/union.area
      print(paste(int.area,"/",union.area,"=",jacc))
    }
    jacc.scores[i,j]=jacc
  }
}
colnames(jacc.scores)=gsub(" ","_",species)
rownames(jacc.scores)=gsub(" ","_",species)
diag(jacc.scores)=NA

write.csv(jacc.scores,file.path(niche_dir,"range_overlap_jaccard.csv"))
# jacc.scores=as.matrix(read.csv(file.path(niche_dir,"range_overlap_jaccard.csv"),row.names=1))



##########################################
## Pairwise geographic distance matrix ###
######### (from Matt Larcombe) ###########
##########################################
library(raster)
library(ape)

phy=read.tree(file.path(TTR_dir,"sortadate50_TTRtips.tre"))


occ=read.csv(file.path(dist_dir,"Veronica_eflora_strictClean_norepeats_20210226.csv"))
occ$species=gsub(" ","_",occ$species)
occ[which(occ$species=="Veronica_plano-petiolata"),"species"] = "Veronica_planopetiolata"
occ[which(occ$species=="Veronica_cockayniana"),"species"] = "Veronica_cockayneana"

species=unique(occ$species)
n.species<-length(species) 

## make a matrix to populate with the data
veronica.distances <- matrix(nrow = n.species, ncol = n.species)
rownames(veronica.distances) <- species
colnames(veronica.distances) <- species

## use the pointDistance function from raster to generate pairwise matrices for all species 
for(i in 1:n.species){
  print(i)
  for(j in 1:n.species){
    pr1 <- as.matrix(occ[occ$species==species[i],c("longitude","latitude")])
    pr2 <- as.matrix(occ[occ$species==species[j],c("longitude","latitude")])
    disM <- pointDistance(pr1, pr2, lonlat=T,allpairs=T)
    veronica.distances[i,j] <- mean(disM) #take the mean as the point distance
    
  }
}

diag(veronica.distances)<-0
saveRDS(veronica.distances, file = file.path(niche_dir,"Veronica_geo_distances.Rdata"))
veronica.distances=readRDS(file.path(niche_dir,"Veronica_geo_distances.Rdata"))
