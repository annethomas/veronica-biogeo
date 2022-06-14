require(ape)
require(dispRity)
require(Claddis)
require(mvMORPH)
require(dplyr)
require(ggplot2)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography"
TTR_dir = file.path(data_dir,"TTR")
bgb_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"biogeography_paper/comparative_phylo")
dist_dir=file.path(data_dir,"distribution_data")

source(file.path(script_dir,"ancestral_slice_functions.R"))

## choose BGB dataset
dataset="merged"
#dataset="ASL"

if(dataset=="merged"){
  out_dir=file.path(bgb_dir,"sortadate50_merged_TTRtips")
  areas=c('L','M')
  area_codes=c('O','M','L','ML')
  focal_areas='M'
} else if(dataset=="ASL"){
  out_dir=file.path(bgb_dir,"sortadate50_ASL_TTRtips")
  areas<-c('L','S','A')
  area_codes<-c('O','A','S','L','AS','AL','SL','ASL')
  focal_areas=c('A','S')
}

##########################
# interactive phylo PCA #
##########################
## only need to do this once for a given set of TTR runs/species

## pPCA for choosing variables
require(phytools)

## TTR params that were unconstrained
traits_to_drop=c(15,16)

## check phylogenetic signal of individual parameters
lambdas=apply(ttr_traits_matrix,2,function(x){phylosig(ttr_tree,x,method='lambda')})
## check correlations between raw variables
pairs(ttr_traits_matrix[,-traits_to_drop])

## run pPCA
ttr.pca.lambda=phyl.pca(ttr_tree, ttr_traits_matrix[,-traits_to_drop],method="lambda")
plot(ttr.pca.lambda)
## choose PCs with most variance
cumsum(diag(ttr.pca.lambda$Eval)/sum(diag(ttr.pca.lambda$Eval))) #PC1-7 ~ 90%
num_PCs=7
## examine loadings (eigenvector)
row.names(ttr.pca.lambda$S)=gsub("Veronica_","",row.names(ttr.pca.lambda$S))
row.names(ttr.pca.lambda$Evec)=full_traits[-traits_to_drop]
# check biplots for correlated variables after filtering by loading below, choose highest loading
biplot(ttr.pca.lambda)
biplot(ttr.pca.lambda,choices=c(2,3))

loading_thresh=sqrt(1/length(full_traits[-traits_to_drop]))
evec_loadings_filtered=ttr.pca.lambda$Evec[,c(1:num_PCs)]
evec_loadings_filtered[which(abs(evec_loadings_filtered) < loading_thresh)]=0
write.csv(evec_loadings_vals,file.path(TTR_dir,"TTR22_eigenvector_loadings_filtered.csv"))
## see hand-annotated excel spreadsheet
eigen_ranks=read.csv(file.path(TTR_dir,"TTR22_eigenvector_loadings_filtered_ranked.csv"))
## top variables:
## strong: tmax3, q1, w11, w21, tmin3
## weak (single instance in PCs 6/7): tmax1, tmean2,tmin2, tmean22



## archived: check unconstrained TTR param proportions
sdm.stats=read.csv(file.path(TTR_dir,"Veronica_sdm_stats.csv"))
sdm.stats$sp=gsub("plano-petiolata","planopetiolata",sdm.stats$sp)
row.names(sdm.stats) = sdm.stats$sp

sdm.stats.lim=sdm.stats[,c(8:31)]
sdm.stats.lim=apply(sdm.stats.lim,2,function(x){return(ifelse(x>1,1,x))})


apply(sdm.stats.lim,2,function(x){length(which(x==1))/length(x)})
# tmax1      tmax2      tmax3      tmax4         q1         q2        w11        w12        ns1        ns2     tmean1     tmean2        w21        w22 
# 0.00000000 0.00000000 0.44444444 0.53703704 0.00000000 0.00000000 0.00000000 0.03703704 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.08333333 
# w23        w24     nsoil1     nsoil2      tmin1      tmin2      tmin3      tmin4    tmean21    tmean22 
# 0.87962963 0.95370370 0.00000000 0.00000000 0.00000000 0.00000000 0.30555556 0.38888889 0.03703704 0.16666667 



