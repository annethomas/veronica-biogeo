require(ape)
require(dispRity)
require(Claddis)
require(mvMORPH)
require(dplyr)
require(ggplot2)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/"
TTR_dir = file.path(data_dir,"TTR")
BGB_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"biogeography_paper/comparative_phylo")
mv_dir=file.path(data_dir,"biogeography_paper/mvMORPH_time_slice")

source(file.path(script_dir,"time_slice/ancestral_slice_functions.R"))

## choose BGB dataset
dataset="merged"
#dataset="ASL"

## note: the median distance calculations and summary stats/models in run_slice() etc are currently calculated 
## from the whole dist_df, and thus for all the focal areas combined.
## If only interested in node-by-node distances/disparity in dist_df, set distdf_only to TRUE 
## If interested in plotting both mountain and lowland disparity but not combined stats, set ML_distdf to TRUE
## For a full run focused on mountains, set ML_distdf to FALSE.

ML_distdf = FALSE ## only affects merged dataset

if(dataset=="merged"){
  run_dir=file.path(BGB_dir,"sortadate50_merged_TTRtips")
  out_dir=file.path(mv_dir,"sortadate50_merged_TTRtips")
  areas=c('L','M')
  area_codes=c('O','M','L','ML')
  if(ML_distdf) {
    focal_areas=c('M','L')
    distdf_only = TRUE
  }  else{
    focal_areas='M'
    distdf_only=FALSE
  }  
} else if(dataset=="ASL"){
  run_dir=file.path(BGB_dir,"sortadate50_ASL_TTRtips")
  out_dir=file.path(mv_dir,"sortadate50_ASL_TTRtips")
  areas<-c('L','S','A')
  area_codes<-c('O','A','S','L','AS','AL','SL','ASL')
  focal_areas=c('A','S')
  distdf_only=FALSE
}

###############################
# read in and prep trait data #
###############################
trait_file<-file.path(trait_dir,"TTR.traits.nex")
tree_file<-file.path(BGB_dir,"sortadate50_TTRtips.tre")

trait_input <- prep_traits_tree(trait_file,tree_file)
ttr_traits_matrix <- trait_input$traits_matrix$matrix_1$matrix
ttr_tree <- trait_input$traits_phy$tree

## if you don't have below file, go to PCA_select_TTRtraits.R and run phyloPCA pipeline and hand annotate
eigen_file <- file.path(TTR_dir,"TTR22_eigenvector_loadings_filtered_ranked.csv")
trait_input_sub <- subset_traits(trait_input=trait_input,eigen_file=eigen_file)


################################################
## load ancestral matrix and make morphospace ##
##          slice and plot disparity          ##
################################################

traitset="TTR9"
model="OUM"
#model = "OU"
model_name=paste0(traitset,model)
bsm_dir="BSMs_run2"

if(model=="OUM") {
  ## make list of BSM-specific model fit files to be read in loop
  bsmOU_files=list.files(file.path(out_dir,"mvmorph_output",bsm_dir),pattern="BSM*")
  prefix=paste0("mvmorph_output_",model_name,"_BSM")
} else{
  ## read in the overall model fit file and prepare morphospace
  mv=readRDS(file.path(mv_dir,paste0("mvmorph_output_",model_name),".Rdata"))
  row.names(mv$mvanc$estimates) <- trait_input_sub$traits_phy$tree$node.label
  morphospace <- get_morphospace(trait_input_sub$traits_matrix,mv)
}

## calculate slices and stats for each BSM (for MCC tree)
load(file.path(run_dir,bsm_dir,"RES_clado_events_tables.Rdata"))
BSM_out <- list()
for(i in 1:length(bsmOU_files)){
  if(model=="OUM"){
    f=bsmOU_files[i]
    bsm=substr(f,nchar(prefix)+1,nchar(f)-nchar(".Rdata"))
    
    mv <- readRDS(file.path(out_dir,"mvmorph_output",bsm_dir,f))
    ## to account for different list structures
    if(is.null(mv$mvanc))  mv$mvanc <- mv$OU$mvanc 
        
    row.names(mv$mvanc$estimates) <- trait_input_sub$traits_phy$tree$node.label
    morphospace <- get_morphospace(trait_input_sub$traits_matrix,mv)
  } else{ bsm = i }
  
  print(paste("BSM",bsm))

  clado_table=RES_clado_events_tables[[as.numeric(bsm)]]
  lm.output <- run_slice(clado_table=clado_table, areas=areas, area_codes=area_codes,
                         focal_areas=focal_areas,tree=ttr_tree,morphospace=morphospace,
                         title_base=paste(model_name,dataset,paste0("BSM",bsm),sep="_"),out_dir=out_dir,
                         distdf_only=distdf_only)
  BSM_out[[bsm]] <- lm.output
}

## save list of BSM slice data
if(ML_distdf) {
  saveRDS(BSM_out,file.path(out_dir,paste0("slice_output/slice_disparity_MLdistdf_output_",length(bsmOU_files),"BSMs.Rdata")))
} else {
  saveRDS(BSM_out,file.path(out_dir,paste0("slice_output/slice_disparity_output_",length(bsmOU_files),"BSMs.Rdata")))
}

## latest non-ML version saved as "goodBSMs" instead of number of BSMs; see slice_BSMs_OUM_202201.R for record of how this was compiled

## 3/12/22 adding new BSMs 

BSM_out1=readRDS(file.path(out_dir,paste0("slice_output/slice_disparity_output_goodBSMs.Rdata")))
names(BSM_out)=paste0(names(BSM_out),"b")
BSM_out_combine=c(BSM_out1,BSM_out)
saveRDS(BSM_out_combine,file.path(out_dir,paste0("slice_output/slice_disparity_output_50BSMs.Rdata")))


BSM_out1=readRDS(file.path(out_dir,paste0("slice_output/slice_disparity_ML_output_44BSMs.Rdata")))
names(BSM_out)=paste0(names(BSM_out),"b")
BSM_out_combine2=c(BSM_out1,BSM_out)
saveRDS(BSM_out_combine2,file.path(out_dir,paste0("slice_output/slice_disparity_MLdistdf_output_50BSMs.Rdata")))
