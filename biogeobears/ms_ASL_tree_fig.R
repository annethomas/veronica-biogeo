library(ape)
library(BioGeoBEARS)

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


run_dir=file.path(BGB_dir,"sortadate50_time_expertBands")
trfn = file.path(BGB_dir,"sortadate50_tree.tre")
resfilename="Veronica_BAYAREALIKE+J_M0_unconstrained_v1.Rdata"
plotfilename=paste0("BAYAREAJ_ASL_pie_traits_",Sys.Date(),".pdf")


geogfn = file.path(run_dir,"band_summary_expert_reduced.txt")
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
legend_specs = read.csv(file.path(tr_plot_dir,"combined_legend2.csv"),stringsAsFactors = FALSE)

tips = dplyr::left_join(data.frame("tip"=tr$tip.label),libs_info[,c("Species","Sample","natural_group","hebe_group","hebe_group_detail","alpine_status")],by=c("tip"="Species"))
tips=dplyr::left_join(tips,hebe_groups)



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

plot(full.tr)
nodelabels(full.tr$node.label)
nodelabels()
tiplabels()
edgelabels()

## align tips, nodes, edges, and posterior probability values
labels.tipnode.full=data.frame(tip.node.num=1:(Ntip(full.tr)+Nnode(full.tr)),lab=c(full.tr$tip.label,full.tr$node.label),
                               type=c(rep("tip",Ntip(full.tr)),rep("node",Nnode(full.tr))),pp=phy_nexus_rbt$metadata$posterior)

full.edges=as.data.frame(full.tr$edge)
full.edges$edge.num=1:nrow(full.edges)
labels.tipnode.full=left_join(labels.tipnode.full,full.edges[,c("V2","edge.num")],by=c("tip.node.num"="V2"))

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


## what plot_biogeobears_results() does (changes order of edges)
bgb.tr.re=reorder(bgb.tr,"pruningwise")

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
## plot biogeobears tree with adjusted tip labels as base, with colored tips for subclade 

pdf(file.path(ms_plot_dir,"drafting",plotfilename),height=10,width=6.5)
#par(mar=c(5.1,4.1,4.1,2.1),mai=c(1.02,0.82,0.82,0.42))
my_plot_BioGeoBEARS_results(res, analysis_titletxt="BAYAREA+J ASL", addl_params=list("j"), 
                            plotwhat="pie", plotlegend=FALSE, tipcex=0.6, pie_tip_statecex=1,
                            statecex=0.7, splitcex=0.6, titlecex=0, label.font=3, tiplabel_adj=0.55,plotsplits=FALSE,
                            tipboxes_text=FALSE,tipcol=hebe_colors[tips$code],edge.width = edge.w,
                            cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

legend(x=-.31,y=16,legend=legend_specs$text,col=legend_specs$color,pt.bg=legend_specs$pt.bg,
       text.font=legend_specs$font,#text.col=legend_specs$text.color,
       pch=legend_specs$pch,cex=.65)



dev.off()

####################
## elevation bars ##
####################
dist_dir=file.path(data_dir,"distribution_data")

alt_lat=read.csv(file.path(dist_dir,"alt_lat_bands.csv"))
alt_lat$species=gsub("Veronica_","V. ",alt_lat$species)
pgj_elev=read.csv(file.path(dist_dir,"elevation/annotated_bands_PGJ_elev.csv"))
pgj_elev$species=gsub("Veronica_","V. ",pgj_elev$species)

occ=read.csv(file.path(dist_dir,"Veronica_eflora_strictClean_norepeats_20210226.csv"))
occ$species=gsub("Veronica","V.", occ$species)
occ$species=gsub("cockayniana","cockayneana", occ$species)
occ$species=gsub("plano-petiolata","planopetiolata", occ$species)

occ_stats=occ %>% group_by(species) %>% 
  summarize(min.elev.eflora=min(elevation,na.rm = TRUE),max.elev.eflora=max(elevation,na.rm = TRUE)) %>%
  left_join(pgj_elev[,c("species","min.elev","max.elev")]) %>% left_join(alt_lat[,c("species","category","merged_category")])

occ_stats_tips=left_join(data.frame("species"=tr$tip.label),occ_stats)
occ_stats_tips$species=factor(occ_stats_tips$species,levels=occ_stats_tips$species)

ggplot(occ_stats_tips,aes(species,min.elev)) + geom_segment(aes(xend=species,yend=max.elev),size=1,color="darkgray") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) + coord_flip()

pdf(file=file.path(ms_plot_dir,"drafting/elev_side_plot.pdf"),height=7,width=2)
elev=ggplot(occ_stats_tips,aes(species,min.elev)) + geom_segment(aes(xend=species,yend=max.elev),size=1,color="darkgray") +
  theme_bw() +  coord_flip() + ylab("Elevation range") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid=element_blank(),
        panel.border=element_blank()) 
elev
dev.off()

## with no axes except elevation label
ggplot(occ_stats_tips,aes(species,min.elev)) + geom_segment(aes(xend=species,yend=max.elev),size=1,color="darkgray") +
  theme_void() +  coord_flip() + theme(axis.text.x=element_text(),axis.line.x=element_line())

#### 


#####
## this plots in order of category and elevation
occ_stats$category=factor(occ_stats$category,levels=c("A","AS","ASL","S","SL","L"))
occ_stats = with(occ_stats,occ_stats[order(category,max.elev),])
occ_stats$species=gsub("Veronica_","V. ",occ_stats$species)
occ_stats$species=factor(occ_stats$species,levels=occ_stats$species)

ggplot(occ_stats,aes(species,min.elev)) + geom_segment(aes(xend=species,yend=max.elev,color=category),size=1) +
  theme_bw() +
  ylab("Elevation (m)") +
  theme(axis.text.x=element_text(angle=90,size=8))

