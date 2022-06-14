library(ape)
library(BioGeoBEARS)
library(dplyr)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
BGB_dir=file.path(data_dir,"biogeobears") # parent bgb dir
TTR_dir = file.path(data_dir,"TTR")
dist_dir=file.path(data_dir,"distribution_data")
ms_plot_dir=file.path("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/manuscripts/Uplift_biogeo_paper/figures_tables")
tr_plot_dir=file.path(data_dir,"biogeography_paper/tree_plotting")


source('~/Cambridge/Research/scripts/biogeography_github/biogeobears/my_plot_BioGeoBEARS_results.R')
source('~/Cambridge/Research/scripts/biogeography_github/util/nodiv_plotting_functions_edited.R')
source('~/Cambridge/Research/scripts/biogeography_github/util/nodiv_color_functions.R')

alt_lat=read.csv(file.path(dist_dir,"alt_lat_bands.csv"))
sdm.stats.c=read.csv(file.path(TTR_dir,"Veronica_sdm_stats_corrected.csv"))

## merged areas with TTR tips vs all bgb tips
TTRtips=FALSE

if(TTRtips){
  run_dir=file.path(BGB_dir,"sortadate50_merged_TTRtips")
  trfn = file.path(BGB_dir,"sortadate50_TTRtips.tre")
  resfilename="Veronica_DEC+J_M0_unconstrained_v1.Rdata"
  plotfilename=paste0("DECJ_ML_pie_traits_TTR_",Sys.Date(),".pdf")
} else{
  run_dir=file.path(BGB_dir,"sortadate50_time_expertBands_merged")
  trfn = file.path(BGB_dir,"sortadate50_tree.tre")
  resfilename="Veronica_DEC+J_M0_unconstrained_v1.Rdata"
  plotfilename=paste0("DECJ_ML_pie_traits_",Sys.Date(),".pdf")
}

geogfn = file.path(run_dir,"band_summary_expert_merged.txt")
time_strat=TRUE ## copy time_periods.txt, manual_dispersal_multipliers.txt and areas_allowed.txt to new run_dir/bgb_dir
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size = 6
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))




## adjust visuals
load(file.path(run_dir,resfilename))

tr=read.tree(trfn)
tr$tip.label=gsub("Veronica_","V. ",tr$tip.label)
tr$tip.label=gsub("plano-petiolata", "planopetiolata",tr$tip.label)

## prepare subclade colors
libs_info=read.csv(file.path(tr_plot_dir,"libs_info_detailed.csv"),stringsAsFactors=FALSE)
libs_info[which(libs_info$Species=="V. vernicosa-2926"),"Species"]="V. vernicosa"
libs_info[which(libs_info$Species=="V. hookeriana-BDM12"),"Species"]="V. hookeriana"
#libs_info$Species=gsub("plano-petiolata", "planopetiolata",libs_info$Species)

hebe_groups = read.csv(file.path(tr_plot_dir,"hebe_groups_simple_colors.csv"),stringsAsFactors = FALSE)
hebe_colors = unique(hebe_groups$color)
legend_specs = read.csv(file.path(tr_plot_dir,"combined_legend.csv"),stringsAsFactors = FALSE)

tips = dplyr::left_join(data.frame("tip"=tr$tip.label),libs_info[,c("Species","Sample","natural_group","hebe_group","hebe_group_detail","alpine_status")],by=c("tip"="Species"))
tips=dplyr::left_join(tips,hebe_groups)

## prepare labels
# note: adj is centered relative to the tips at (0.5,0.5) while label.offset starts from the tips, but it's unclear how well they correspond
sdm.stats.c$sp=gsub("_"," ",sdm.stats.c$sp)
alt_lat2=alt_lat
alt_lat2$species=gsub("Veronica_","V. ",alt_lat2$species)
sdm.stats.tips=dplyr::left_join(data.frame(sp=tr$tip.label),sdm.stats.c) %>% 
  left_join(alt_lat2[,c("species","NI","SI")],by=c("sp"="species"))

## prepare colors
temp.mat=as.matrix(sdm.stats.tips[,c("tmin2.back","tmean2.back")])
rownames(temp.mat)=sdm.stats.tips$sp
colscale=choose.colors(temp.mat)
zlim=attr(colscale,"zlim")
tmin.col=create.cols(temp.mat[,1],col=colscale)
tmin.col[is.na(tmin.col)]="grey50"
tmean.col=create.cols(temp.mat[,2],col=colscale)
tmean.col[is.na(tmean.col)]="grey50"

#########
## get posterior probability support values
#########
library(rBt)
phy_nexus_path=file.path(BGB_dir,"rank50_strictEstimate.species.ann.tre")
phy_nexus_rbt=rBt::read.annot.beast(phy_nexus_path)
phy_nexus_rbt$tip.label=paste("V.",phy_nexus_rbt$tip.label)
phy_nexus_rbt$tip.label=gsub("plano-petiolata", "planopetiolata",phy_nexus_rbt$tip.label)

## surrogate ape tree with original tips
full.tr=ape::read.nexus(phy_nexus_path)
full.tr$node.label=paste0("n",1:full.tr$Nnode)
full.tr$tip.label=paste("V.",full.tr$tip.label)
full.tr$tip.label=gsub("plano-petiolata", "planopetiolata",full.tr$tip.label)

# plot(full.tr)
# nodelabels(full.tr$node.label)
# nodelabels()
# tiplabels()
# edgelabels()

## align tips, nodes, edges, and posterior probability values
labels.tipnode.full=data.frame(tip.node.num=1:(Ntip(full.tr)+Nnode(full.tr)),lab=c(full.tr$tip.label,full.tr$node.label),
                               type=c(rep("tip",Ntip(full.tr)),rep("node",Nnode(full.tr))),pp=phy_nexus_rbt$metadata$posterior)

full.edges=as.data.frame(full.tr$edge)
full.edges$edge.num=1:nrow(full.edges)
labels.tipnode.full=left_join(labels.tipnode.full,full.edges[,c("V2","edge.num")],by=c("tip.node.num"="V2"))
## test: identify which edges to thicken and test plot on full tree
# edge.order=dplyr::arrange(labels.tipnode.full,edge.num)
# edge.w=ifelse(edge.order$pp>=0.9,2,1)
# edge.w[is.na(edge.w)]=1
# plot(full.tr,edge.width=edge.w)

## reduce to biogeobears version of tree
tips.to.drop=setdiff(full.tr$tip.label,tr$tip.label)
bgb.tr=drop.tip(full.tr,tips.to.drop)

setdiff(full.tr$node.label,bgb.tr$node.label)

labels.tipnode.bgb=data.frame(tip.node.num=1:(Ntip(bgb.tr)+Nnode(bgb.tr)),lab=c(bgb.tr$tip.label,bgb.tr$node.label),
                               type=c(rep("tip",Ntip(bgb.tr)),rep("node",Nnode(bgb.tr))))
bgb.edges=as.data.frame(bgb.tr$edge)
names(bgb.edges)=c("node1","node2")
bgb.edges$edge.num=1:nrow(bgb.edges)
labels.tipnode.bgb=left_join(labels.tipnode.bgb,bgb.edges[,c("node2","edge.num")],by=c("tip.node.num"="node2"))

## align new tree with original posterior probability values
labels.tipnode.join=left_join(labels.tipnode.full,labels.tipnode.bgb[,c("lab","tip.node.num","edge.num")],by="lab",suffix=c(".full",".bgb"))

## identify which edges to thicken
# edge.order=dplyr::arrange(labels.tipnode.join,edge.num.bgb) %>% filter(!is.na(edge.num.bgb))
# edge.w=ifelse(edge.order$pp>=0.9,2,1)
# edge.w[is.na(edge.w)]=1
# ## test plot
# plot(bgb.tr,edge.width=edge.w)
# nodelabels(bgb.tr$node.label)

## what plot_biogeobears_results() does
bgb.tr.re=reorder(bgb.tr,"pruningwise")
# plot(bgb.tr.re,edge.width=edge.w) ## changes order of edges

## reassign edge numbers
bgb.edges.re = as.data.frame(bgb.tr.re$edge)
names(bgb.edges.re)=c("node1","node2")
bgb.edges.re$edge.num=1:nrow(bgb.edges.re)

bgb.edges.join=left_join(bgb.edges.re,bgb.edges[,c("node2","edge.num")],by=c("node2"),suffix=c(".re",".orig"))
labels.tipnode.join=left_join(labels.tipnode.join,bgb.edges.join[,c("edge.num.orig","edge.num.re")],by=c("edge.num.bgb"="edge.num.orig"))

## identify which edges to thicken
edge.order=dplyr::arrange(labels.tipnode.join,edge.num.re) %>% filter(!is.na(edge.num.re))
edge.w=ifelse(edge.order$pp>=0.9,2,1)
edge.w[is.na(edge.w)]=1
## test plot
plot(bgb.tr.re,edge.width=edge.w)



#############################################################################################
## Version #1
## plot biogeobears tree with adjusted tip labels as base, with colored tips for subclade 
dev.off() ## to clear graphics device since bgb plot gets confused
pdf(file.path(tr_plot_dir,plotfilename),height=10,width=8)
#par(mar=c(5.1,4.1,4.1,2.1),mai=c(1.02,0.82,0.82,0.42))
my_plot_BioGeoBEARS_results(res, analysis_titletxt="DEC+J merged", addl_params=list("j"), 
                            plotwhat="pie", plotlegend=FALSE, label.offset=0.52, tipcex=0.6, pie_tip_statecex=1,
                            statecex=0.7, splitcex=0.6, titlecex=0, label.font=3, tiplabel_adj=0.55,plotsplits=FALSE,
                            tipboxes_text=FALSE,tipcol=hebe_colors[tips$code],edge.width=edge.w,
                            cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

## add labels to tree
tip.pt.adj=0.5+0.05 ## to bump starting adj point based on tiplabel_adj above
adj.width=0.13 ## distance from last column
adj1=adj.width
adj2=2*adj.width
adj3=3*adj.width
adj4=4*adj.width
#adj5=5*adj.width

#tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj1,col=hebe_colors[tips$code])
tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj1,col=ifelse(sdm.stats.tips$NI,"black","white"))
tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj2,col=ifelse(sdm.stats.tips$SI,"black","white"))

tiplabels(pch=15,cex=0.9,adj=tip.pt.adj+adj3,col=tmin.col)
tiplabels(pch=15,cex=0.9,adj=tip.pt.adj+adj4,col=tmean.col)

## add text above the trait columns
par(xpd=TRUE)
xx=max(dispRity::tree.age(tr)$ages) #6.354
text(x=xx+adj1+.03,y=Ntip(tr)*1.015,lab="N",cex=0.75)
text(x=xx+adj2+.04,y=Ntip(tr)*1.015,lab="S",cex=0.75)
text(x=xx+adj3+.08,y=Ntip(tr)*1.03,lab="tmin2",srt=70,cex=0.75) ## angled words have to be adjusted farther up and over
text(x=xx+adj4+.145,y=Ntip(tr)*1.0425,lab="tmeanN2",srt=70,cex=0.75)

legend(x=-.25,y=28.5,legend=legend_specs$text,col=legend_specs$color,pt.bg=legend_specs$pt.bg,
       text.font=legend_specs$font,#text.col=legend_specs$text.color,
       pch=legend_specs$pch,cex=.8)

plot.new()
add_legend(zlim,col=colscale)
text("Temp (\u00B0C)",x=1,y=min(zlim)-2)

dev.off()

####################################################################################################
## Version #2
## plot biogeobears tree with adjusted tip labels as base, with colored dots for subclades
pdf(file.path(tr_plot_dir,plotfilename),height=10,width=8)
#par(mar=c(5.1,4.1,4.1,2.1),mai=c(1.02,0.82,0.82,0.42))
my_plot_BioGeoBEARS_results(res, analysis_titletxt="DEC+J merged", addl_params=list("j"), 
                         plotwhat="pie", plotlegend=FALSE, label.offset=0.65, tipcex=0.6, pie_tip_statecex=1,
                         statecex=0.7, splitcex=0.6, titlecex=0, label.font=3, tiplabel_adj=0.55,plotsplits=FALSE,
                         tipboxes_text=FALSE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

## add labels to tree
tip.pt.adj=0.5+0.05 ## to bump starting adj point based on tiplabel_adj above
adj.width=0.13 ## distance from last column
adj1=adj.width
adj2=2*adj.width
adj3=3*adj.width
adj4=4*adj.width
adj5=5*adj.width

tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj1,col=hebe_colors[tips$code])
tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj2,col=ifelse(sdm.stats.tips$NI,"black","white"))
tiplabels(pch=19,cex=0.6,adj=tip.pt.adj+adj3,col=ifelse(sdm.stats.tips$SI,"black","white"))

tiplabels(pch=15,cex=0.9,adj=tip.pt.adj+adj4,col=tmin.col)
tiplabels(pch=15,cex=0.9,adj=tip.pt.adj+adj5,col=tmean.col)

## add text above the trait columns
par(xpd=TRUE)
xx=max(dispRity::tree.age(tr)$ages) #6.354
text(x=xx+adj2+.03,y=Ntip(tr)*1.015,lab="N",cex=0.75)
text(x=xx+adj3+.04,y=Ntip(tr)*1.015,lab="S",cex=0.75)
text(x=xx+adj4+.08,y=Ntip(tr)*1.03,lab="tmin2",srt=70,cex=0.75) ## angled words have to be adjusted farther up and over
text(x=xx+adj5+.145,y=Ntip(tr)*1.0425,lab="tmeanN2",srt=70,cex=0.75)

legend(x=-.4,y=32,legend=legend_specs$text,col=legend_specs$color,pt.bg=legend_specs$pt.bg,
       text.font=legend_specs$font,#text.col=legend_specs$text.color,
       pch=legend_specs$pch,cex=.8)

plot.new()
add_legend(zlim,col=colscale)
text("Temp (\u00B0C)",x=1,y=min(zlim)-2)

dev.off()

