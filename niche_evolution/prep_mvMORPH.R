## prepare input for running mvMORPH on the cluster
require(ape)
require(dispRity)
require(Claddis)
require(mvMORPH)
require(dplyr)
require(ggplot2)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github"
TTR_dir = file.path(data_dir,"TTR")
BGB_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"biogeography_paper/comparative_phylo")
mv_dir=file.path(data_dir,"biogeography_paper/mvMORPH_time_slice")

source(file.path(script_dir,"time_slice/ancestral_slice_functions.R"))

########################
## original code for traits nexus file below, however it involved hand editing and hasn't been retested
########################
sdm.stats=read.csv(file.path(TTR_dir,"Veronica_sdm_stats.csv"))
sdm.stats$sp=gsub("plano-petiolata","planopetiolata",sdm.stats$sp)
row.names(sdm.stats) = sdm.stats$sp

sdm.stats.lim=sdm.stats[,c(8:31)]
sdm.stats.lim=apply(sdm.stats.lim,2,function(x){return(ifelse(x>1,1,x))})


### writing nexus file for continuous data 
sdm.stats.list=list()
for(i in 1:nrow(sdm.stats.lim)){
  sp=row.names(sdm.stats.lim)[i]
  sdm.stats.list[[sp]]=unlist(sdm.stats.lim)
}

#write.nexus.data(sdm.stats.list,file=file.path(TTR_dir,"test.TTR.nex"),format="continuous",datablock=FALSE)

## hand edit: remove "interleave=no", breaks read_nexus_matrix for continuous data; make sure there's no ; directly after CONTINUOUS
## write.csv(data.frame(colnames(sdm.stats.lim)),file=file.path(trait_dir,"TTR_traits.csv")) ## for manually adding trait names to nexus file, not strictly necessary

###############################
# read in and prep trait data #
###############################
trait_file<-file.path(TTR_dir,"TTR.traits.nex")
tree_file<-file.path(BGB_dir,"sortadate50_TTRtips.tre")

trait_input <- prep_traits_tree(trait_file,tree_file)
ttr_traits_matrix <- trait_input$traits_matrix$matrix_1$matrix
ttr_tree <- trait_input$traits_phy$tree

## if you don't have below file, go to PCA_select_TTRtraits.R and run phyloPCA pipeline and hand annotate
eigen_file <- file.path(TTR_dir,"TTR22_eigenvector_loadings_filtered_ranked.csv")
traitset="TTR9" ## number of traits
#traitset="TTR1"
if(traitset=="TTR9") {
  trait_input_sub <- subset_traits(trait_input=trait_input,eigen_file=eigen_file) ## default is top 9 traits
} else if(traitset=="TTR1"){
  trait_input_sub <- subset_traits(trait_input=trait_input,eigen_file=eigen_file,num_to_keep = 1) ## default is top 9 traits
}

########################
## create input files ##
########################

## choose evolution model type
#model_type="single" ## prepares input for run_mvMORPH_single.R
model_type="multi"   ## prepares input for run_mvMORPH_multiBSM.R

## choose BGB dataset (only applicable for multi-regime)
dataset="merged"
#dataset="ASL"

bsmdir="BSMs_run2" ## currently only for merged dataset (reruns); BSMs_run1 output is identical to what's in the main run_dir
bsmdir="BSMs_run1"

if(model_type=="single"){
  
  input=list()
  input$traits_matrix=trait_input_sub$traits_matrix$matrix_1$matrix
  input$traits_phy=trait_input_sub$traits_phy
  saveRDS(input,file=file.path(mv_dir,paste0("mvmorph_input_",traitset,".Rdata")))
  
} else if(model_type=="multi"){
  
  if(dataset=="merged"){
    run_dir=file.path(BGB_dir,"sortadate50_merged_TTRtips",bsmdir)
    out_dir=file.path(mv_dir,"sortadate50_merged_TTRtips")
    areas=c('L','M')
    area_codes=c('O','M','L','ML')
    focal_areas='M'
  } else if(dataset=="ASL"){
    run_dir=file.path(BGB_dir,"sortadate50_ASL_TTRtips")
    out_dir=file.path(mv_dir,"sortadate50_ASL_TTRtips")
    areas<-c('L','S','A')
    area_codes<-c('O','A','S','L','AS','AL','SL','ASL')
    focal_areas=c('A','S')
  }
  
  #### make simmap trees from BSMs
  source(file.path(script_dir,'niche_evolution/make.simmap.BGB.original.R'))
  source(file.path(script_dir,'niche_evolution/rpanda_stratified_BGB_to_tables.R'))
  
  load(file.path(run_dir,"RES_clado_events_tables.Rdata"))
  
  bsm_idx=1:length(RES_clado_events_tables)
  ##bsm_idx=2:5 ## for testing OUM on slow dataset (ASL)
  
  for(i in bsm_idx){
    print(i)
    smap<-stratified_BGB_to_tables(ttr_tree,RES_clado_events_tables,i)
    simmap=make.simmap.BGB(anc.phylo=ttr_tree,subclade.phylo=ttr_tree,ana.events=smap$ana.int,clado.events=smap$clado.int,return.mat=FALSE)
    plot(simmap$geo.simmap)
    
    input=list()
    input$traits_matrix=trait_input_sub$traits_matrix$matrix_1$matrix
    input$simmap=simmap$geo.simmap
  
    saveRDS(input,file=file.path(out_dir,"mvmorph_simmap_input",bsmdir,paste0("mvmorph_siminput_",traitset,"_BSM",i,".Rdata")))
  }
}

