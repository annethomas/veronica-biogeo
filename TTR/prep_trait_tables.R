## generates several modified/expanded trait/species tables from existing ones and saves as csvs (only run once)
require(ape)
require(nlme)
require(dplyr)
require(ggplot2)
require(phytools)

#script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/"
#install.packages(file.path(script_dir,"TTR.sdm_0.4.tar.gz"),repos=NULL,type="source")
library(TTR.sdm)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github"
TTR_dir = file.path(data_dir,"TTR")
bgb_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"comparative_phylo")
dist_dir=file.path(data_dir,"distribution_data")


#############################################
## add habitat categories to alt_lat table ##
#############################################
assign_band=function(data.row,merged=FALSE){
  band=c()
  if(merged){
    if(data.row$alpine | data.row$subalpine) band=paste0(band,"M")
    if(data.row$lowland) band=paste0(band,"L")
  } else{
    if(data.row$alpine) band=paste0(band,"A")
    if(data.row$subalpine) band=paste0(band,"S")
    if(data.row$lowland) band=paste0(band,"L")
  }
  
  return(band)
}

alt_lat=read.csv(file.path(dist_dir,"alt_lat_bands.csv"))
for(i in 1:nrow(alt_lat)){
  alt_lat[i,"category"]=assign_band(alt_lat[i,])
  alt_lat[i,"merged_category"]=assign_band(alt_lat[i,],merged=TRUE)
  
}
#write.csv(alt_lat,file.path(dist_dir,"alt_lat_bands.csv"),row.names=FALSE)

#################################
## join traits and habitat info #
#################################

## from phylo_PCA_TTRtraits.R, for selecting traits and matching to phylogeny
traitset="TTR9"
trait_input=readRDS(file=file.path(TTR_dir,paste0("mvmorph_input_",traitset,".Rdata")))

traits=as.data.frame(trait_input$traits_matrix)
traits$species=row.names(traits)

traits_bands=left_join(traits,alt_lat,by="species")
traits_bands=dplyr::relocate(traits_bands,"species")

#write.csv(traits_bands,file.path(TTR_dir,"TTR9_traits_bands.csv"),row.names=FALSE)

#############################################
## correct unconstrained vars in sdm stats ## 
#############################################

sdm.stats=read.csv(file.path(TTR_dir,"Veronica_sdm_stats.csv"))

## change unconstrained params 

## thorn_range: the 0-1 normalization scales of the env data on which TTR was run (this version requires TTR package)
thorn_range = list(Tmax=fTMAX,Q=c(1e-5,20),W1=c(0,100),Nshoot=c(0,0.05),Tmean=fTMEAN,
                   W2=c(0,100),Nsoil=c(10,1500),Tmin=fTMIN,Tmean2=fTMEAN)

colnames(sdm.stats)[8:31]
unconstr=sdm.stats[,8:31]>1
apply(unconstr,2,any)
sdm.stats.corrected=sdm.stats

which(sdm.stats$w24>1)
sdm.stats.corrected$tmax3.back[sdm.stats$tmax3>1] = max(thorn_range$Tmax)/10
sdm.stats.corrected$tmax4.back[sdm.stats$tmax4>1] = max(thorn_range$Tmax)/10
sdm.stats.corrected$w12.back[sdm.stats$w12>1] = max(thorn_range$W1)
sdm.stats.corrected$w22.back[sdm.stats$w22>1] = max(thorn_range$W2)
sdm.stats.corrected$w23.back[sdm.stats$w23>1] = max(thorn_range$W2)
sdm.stats.corrected$w24.back[sdm.stats$w24>1] = max(thorn_range$W2)
sdm.stats.corrected$tmin3.back[sdm.stats$tmin3>1] = max(thorn_range$Tmin)/10
sdm.stats.corrected$tmin4.back[sdm.stats$tmin4>1] = max(thorn_range$Tmin)/10
sdm.stats.corrected$tmean21.back[sdm.stats$tmean21>1] = max(thorn_range$Tmean)/10
sdm.stats.corrected$tmean22.back[sdm.stats$tmean22>1] = max(thorn_range$Tmean)/10

write.csv(sdm.stats.corrected,file=file.path(TTR_dir,"Veronica_sdm_stats_corrected.csv"),row.names=FALSE)



###############################################
## add back corrected traits to traits_bands ##
###############################################

traits_bands=read.csv(file.path(TTR_dir,"TTR9_traits_bands.csv")) 
table(traits_bands$category)
# A  AS ASL   L   S  SL 
# 12  43   7  12   5  19
table(traits_bands$merged_category)
# L  M ML 
# 12 60 26 

traits_bands$species=gsub("Veronica","V.",traits_bands$species)
## add corrected back-transformed
traits_bands_back=dplyr::left_join(traits_bands,sdm.stats.corrected[,c(1,32:55)],by=c("species"="sp"))
#write.csv(traits_bands_back,file=file.path(TTR_dir,"TTR9_traits_bands_back.csv"))
