## compare different niche evolution models with TTR traits

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github"
TTR_dir = file.path(data_dir,"TTR")
BGB_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"biogeography_paper/comparative_phylo")
mv_dir=file.path(data_dir,"biogeography_paper/mvMORPH_time_slice")

source(file.path(script_dir,'util/nodiv_plotting_functions_edited.R'))

####################
## compare models ##
####################
## note: the list structure changed when writing cluster scripts for individual and BSM runs, so will need to be careful
TTR9_OU <- readRDS(file.path(mv_dir,paste0("mvmorph_output_TTR9OU.Rdata")))
TTR9_BM <- readRDS(file.path(mv_dir,paste0("mvmorph_output_TTR9BM.Rdata")))
TTR9_EB <- readRDS(file.path(mv_dir,paste0("mvmorph_output_TTR9EB.Rdata")))
#TTR9_BMM <- readRDS(file.path(mv_dir,paste0("mvmorph_output_TTR9BMM.Rdata"))) # unreliable
TTR9_OUM <- readRDS(file.path(mv_dir,paste0("mvmorph_output_TTR9OUM.Rdata")))

model_list=list()
model_list[["OU"]] <- TTR9_OU[["OU"]]
model_list[["BM"]] <- TTR9_BM[["BM"]]
model_list[["EB"]] <- TTR9_EB[["EB"]]
model_list[["OUM"]][["mvfit"]] <- TTR9_OUM$mvfit
model_list[["OUM"]][["mvanc"]] <- TTR9_OUM$OU$mvanc 

for(model in names(model_list)){
  print(model)
  print(model_list[[model]]$mvfit$AICc)
}

library(nodiv)
mvanc_mat=TTR9_OUM$OU$mvanc$estimates
mvanc_mat=TTR9_OU$OU$mvanc$estimates
plot_nodes_phylo(variable=mvanc_mat[,"tmin3"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"q1"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"tmax3"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"w21"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"w11"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"tmean22"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"tmean2"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"tmin2"], tree=ttr_tree)
plot_nodes_phylo(variable=mvanc_mat[,"tmax1"], tree=ttr_tree)

plot_nodes_tips_phylo(node.variable=mvanc_mat[,"tmin3"], tip.variable= ttr_traits_matrix[,"tmin3"],tree=ttr_tree,show.tip.label = TRUE)
colnames(mvanc_mat)


mvanc_mat1=TTR9_OU$OU$mvanc$estimates
mvanc_mat2=TTR9_OUM$OU$mvanc$estimates

pdf(file = file.path(mv_dir,"TTR9_OUvsOUM_traits.pdf"))
for(tr in colnames(mvanc_mat1)){
  plot_nodes_phylo(variable=mvanc_mat1[,tr], tree=ttr_tree,main=paste("TTR9 OU",tr))
  plot_nodes_phylo(variable=mvanc_mat2[,tr], tree=ttr_tree,main=paste("TTR9 OU_ML",tr))
}
dev.off()



pdf(file = file.path(mv_dir,"TTR9_OUvsOUM_traits2.pdf"),height=12,width=12)
for(tr in colnames(mvanc_mat1)){
  plot_nodes_tips_phylo(node.variable=mvanc_mat1[,tr],tip.variable= ttr_traits_matrix[,tr], 
                        tree=ttr_tree,main=paste("TTR9 OU",tr),show.tip.label = TRUE)
  plot_nodes_tips_phylo(node.variable=mvanc_mat2[,tr],tip.variable= ttr_traits_matrix[,tr], 
                        tree=ttr_tree,main=paste("TTR9 OU_ML",tr),show.tip.label = TRUE)
}
dev.off()


### look at mvMORPH output reliability across BSMs (for OUM)
results=data.frame(BSM=character(0),aicc=numeric(0),reliability=numeric(0))
mvout_files=list.files(file.path(mv_dir,"sortadate50_merged_TTRtips/mvmorph_output"),pattern="*.Rdata")
for(f in mvout_files){
  mv=readRDS(file.path(mv_dir,"sortadate50_merged_TTRtips/mvmorph_output",f))
  aicc=mv$mvfit$AICc
  reliable=mv$mvfit$hess.values
  newrow=data.frame(BSM=f,aicc=aicc,reliability=reliable)
  results=rbind(results,newrow)
}

