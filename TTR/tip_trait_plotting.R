library(ape)
library(BioGeoBEARS)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
BGB_dir=file.path(data_dir,"biogeobears") # parent bgb dir
TTR_dir = file.path(data_dir,"TTR")
dist_dir=file.path(data_dir,"distribution_data")
plot_dir=file.path("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/manuscripts/Uplift_biogeo_paper/figures_tables")


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
  plotfilename="DECJ_ML_pie_traits_TTR.pdf"
} else{
  run_dir=file.path(BGB_dir,"sortadate50_time_expertBands_merged")
  trfn = file.path(BGB_dir,"sortadate50_tree.tre")
  resfilename="Veronica_DEC+J_M0_unconstrained_v1.Rdata"
  plotfilename="DECJ_ML_pie_traits.pdf"
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

pdf(file.path(plot_dir,plotfilename),height=10,width=8)
#par(mar=c(5.1,4.1,4.1,2.1),mai=c(1.02,0.82,0.82,0.42))
my_plot_BioGeoBEARS_results(res, analysis_titletxt="DEC+J merged", addl_params=list("j"), 
                         plotwhat="pie", plotlegend=FALSE, label.offset=0.52, tipcex=0.6, pie_tip_statecex=1,
                         statecex=0.7, splitcex=0.6, titlecex=0, label.font=3, tiplabel_adj=0.55,plotsplits=FALSE,
                         tipboxes_text=FALSE,
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

## add labels to tree
# note: adj is centered relative to the tips at (0.5,0.5) while label.offset starts from the tips, but it's unclear how well they correspond
sdm.stats.c$sp=gsub("_"," ",sdm.stats.c$sp)
tr$tip.label = gsub("plano-petiolata","planopetiolata",tr$tip.label)
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

## add to tree
tip.pt.adj=0.5+0.05 ## to bump starting adj point based on tiplabel_adj above
adj.width=0.13 ## distance from last column
adj1=adj.width
adj2=2*adj.width
adj3=3*adj.width
adj4=4*adj.width

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

plot.new()
add_legend(zlim,col=colscale)
text("Temp (°C)",x=1,y=min(zlim)-2)

dev.off()



### old
###############################
# read in and prep trait data #
###############################
trfn_TTR = file.path(BGB_dir,"sortadate50_TTRtips.tre")
phy=read.tree(trfn_TTR)
phy$tip.label = gsub("Veronica","V.",phy$tip.label)
phy$tip.label = gsub("plano-petiolata","planopetiolata",phy$tip.label)

sdm.stats.c=read.csv(file.path(TTR_dir,"Veronica_sdm_stats_corrected.csv"))

# extrasp.phy=setdiff(phy$tip.label,sdm.stats.c$sp)
# phy.red=drop.tip(phy,extrasp.phy)
# extrasp.df=setdiff(sdm.stats.c$sp,phy$tip.label)
# sdm.stats.c=dplyr::filter(sdm.stats.c,!sp %in% extrasp.df)

sdm.stats.tips=dplyr::left_join(data.frame(sp=phy$tip.label),sdm.stats.c)

## picante
library(picante)
color.plot.phylo(phy.red,sdm.stats.tips,"tmin2.back","sp")

## phytools
library(phytools)
temp.mat=as.matrix(sdm.stats.tips[,c("tmin2.back","tmean2.back")])
temp.mat=as.matrix(sdm.stats.tips[,c("tmin1.back","tmin2.back","tmin3.back","tmin4.back","tmean2.back")])

rownames(temp.mat)=sdm.stats.tips$sp
phylo.heatmap(phy,temp.mat,fsize=c(0.5,0.8,0),split=c(0.75,0.25))

##
max(temp.mat)
min(temp.mat)
palette <- colorRampPalette(c("blue", "white", "red"))(n = 300)
show_col(palette)

##########
## nodiv
library(nodiv)
source('~/Cambridge/Research/scripts/biogeography/nodiv_plotting_functions_edited.R')

pdf(file = file.path(TTR_dir,"TTR9_OUvsOUM_traits2.pdf"),height=12,width=12)
for(tr in colnames(mvanc_mat1)){
  plot_nodes_tips_phylo(node.variable=mvanc_mat1[,tr],tip.variable= ttr_traits_matrix[,tr], 
                        tree=ttr_tree,main=paste("TTR9 OU",tr),show.tip.label = TRUE)
  plot_nodes_tips_phylo(node.variable=mvanc_mat2[,tr],tip.variable= ttr_traits_matrix[,tr], 
                        tree=ttr_tree,main=paste("TTR9 OU_ML",tr),show.tip.label = TRUE)
}
dev.off()

