library(dplyr)
library(ggplot2)

script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts"
data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
hyb_dir=file.path(data_dir,"HybPiper")
dist_dir=file.path(data_dir,"distribution_data")
TTR_dir = file.path(data_dir,"TTR")
niche_dir=file.path(data_dir,"biogeography_paper/niche_overlap")
mv_dir=file.path(data_dir,"biogeography_paper/mvMORPH_time_slice")
slice_dir=file.path(mv_dir,"sortadate50_merged_TTRtips/slice_output")

plot_dir=file.path("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/manuscripts/Uplift_biogeo_paper/figures_tables")
tr_plot_dir=file.path(data_dir,"biogeography_paper/tree_plotting")

#dir.create(plot_dir)

########################################
## Fig1: biogeobears + TTR trait tree ##
########################################
## see TTR/ms_tree_fig.R

##############################
## Fig2: TTR trait boxplots ##
##############################
## from investigate_TTR_traits.R
traits_bands_back=read.csv(file.path(TTR_dir,"TTR9_traits_bands_back.csv"))

## get categorywise max for plotting purposes
tmin.summ=traits_bands_back %>% group_by(merged_category) %>%
  summarize(max.tmin=max(tmin2.back))

tmean.summ=traits_bands_back %>% group_by(merged_category) %>%
  summarize(max.tmean=max(tmean2.back))

## for overlaying geom_jitter, add col with outliers excluded
# https://stackoverflow.com/questions/61335789/how-to-exclude-outliers-when-using-geom-boxplot-geom-jitter-in-r#:~:text=You%20can%20hide%20the%20outliers,%25%3E%25%20mutate(cty.
traits_bands_back.filt=traits_bands_back %>%
  select(species,merged_category,tmin2.back,tmean2.back) %>%
  group_by(merged_category) %>%
  mutate(tmin2_filtered = case_when(tmin2.back - quantile(tmin2.back)[4] > 1.5*IQR(tmin2.back) ~ NA_real_,
                                  quantile(tmin2.back)[2] - tmin2.back > 1.5*IQR(tmin2.back) ~ NA_real_,
                                  TRUE ~ tmin2.back),
         tmean2_filtered = case_when(tmean2.back - quantile(tmean2.back)[4] > 1.5*IQR(tmean2.back) ~ NA_real_,
                                    quantile(tmean2.back)[2] - tmean2.back > 1.5*IQR(tmean2.back) ~ NA_real_,
                                    TRUE ~ tmean2.back)) 
  ggplot() + geom_boxplot(aes(drv, cty)) + geom_jitter(aes(drv, cty_filtered))

# #lower 
# max(min(test$tmin2.back,na.rm=T), as.numeric(quantile(test$tmin2.back, 0.25)) - (IQR(test$tmin2.back)*1.5))
# #upper
# min(max(test$tmin2.back,na.rm=T), as.numeric(quantile(test$tmin2.back, 0.75)) + (IQR(test$tmin2.back)*1.5))

## plot tmin and tmean with significant differences from phylANOVA/PGLS 
g=ggplot(data=traits_bands_back.filt,aes(x=factor(merged_category,levels=c("M","ML","L")))) +  
  theme_classic(base_size=15) + theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("M"="Mountain","ML"="Both","L"="Lowland")) +
  scale_color_manual(values=c("green","blue","gray"))#+ theme(plot.margin = margin(t=1,b=.5,unit="cm"))

a=g + geom_boxplot(aes(y=tmin2.back)) + 
  geom_jitter(aes(y=tmin2_filtered,color=merged_category),width=0.2,show.legend = FALSE) + 
  ylab("Minimum temperature (\u00B0C)") + 
  geom_text(data=tmin.summ,aes(label=c("b","a","b"),x=merged_category,
                               y=max.tmin,vjust=-.8),size=5) +
  ylim(NA,max(traits_bands_back$tmin2.back)+2)

b=g + geom_boxplot(aes(y=tmean2.back)) + 
  geom_jitter(aes(y=tmean2_filtered,color=merged_category),width=0.2,show.legend = FALSE) + 
  ylab("Mean temperature (\u00B0C)") + 
  geom_text(data=tmean.summ,aes(label=c("b","a","b"),x=merged_category, 
                                y=max.tmean,vjust=-.8),size=5) +
  ylim(NA,max(traits_bands_back$tmean2.back)+2)

#pdf(file.path(bgb_dir,"TTRtip_traits_bands_plots_ms.pdf"),height=4.5,width=7)
pdf(file.path(plot_dir,"TTRtip_traits_bands_plots_ms.pdf"),height=4.5,width=7)


ggpubr::ggarrange(a,b,labels=c("(a)","(b)"),label.x=0.02, label.y=1,
                  font.label=list(size=14,face="plain"),legend="none") 

dev.off()


####################################
## Fig 3: binned rates (merged M) ##
####################################
## from binning_figures.R
bgb_dir=file.path(data_dir,"biogeobears")

source("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/biogeobears/binning/bin_plot_BSMs.R")

## merged mountains and lowlands
plot_dir=file.path(bgb_dir,"biogeobears_output_100trees/sortadate50_merged/binned_plots/bluegreen/")
load(file.path(bgb_dir,"biogeobears_output_100trees/sortadate50_merged/sortadate50_merged_binned_output_multi.Rdata"))
binned.CI=output.list[[1]]

## change Os to NA for mountains in binned.CI before 3.5 mya
binned.CI$M$insitu.rates.CI=pad.NA(binned.CI$M$insitu.rates.CI,0.5)
binned.CI$M$dispersal.into.rates.CI = pad.NA(binned.CI$M$dispersal.into.rates.CI,0.5)
binned.CI$M$insitu.time.bins.cumsum.CI=pad.NA(binned.CI$M$insitu.time.bins.cumsum.CI,0.5,4)
binned.CI$M$dispersal.into.time.bins.cumsum.CI = pad.NA(binned.CI$M$dispersal.into.time.bins.cumsum.CI,0.5,4)

## divide by two for y-axis to be per million years instead of 0.5 my
binned.CI.My=binned.CI
binned.CI.My$M$insitu.rates.CI=lapply(binned.CI.My$M$insitu.rates.CI,function(x) x/2)
binned.CI.My$M$dispersal.into.rates.CI=lapply(binned.CI.My$M$dispersal.into.rates.CI,function(x) x/2)

# save original par settings
opar=par(no.readonly = TRUE) # restore with par(opar)

##########
## version with two panels
##########
pdf(file=file.path(ms_plot_dir,"drafting/paneled_binned_plots_padded.pdf"),width=8,height=5)

## simple grid
m1=matrix(c(1,2),ncol=2)
layout(m1)
#par(mar=c(4,4,1,.5)) #1=bottom, 2=left, 3=top, 4=right

# combined events plot
num.bins=length(binned.CI$M$dispersal.into.time.bins.cumsum.CI)
bin.width=0.5
trunc.year=8

par(mar=c(4,4,1,.5),oma=c(1,1,1,1))
plot_binned_multi_var(binned.CI,vars=c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),
                      title=NULL,yaxis.title = "Cumulative events",legend.override=c("Lowland cladogenesis","Mountain cladogenesis","Lowland colonization","Mountain colonization"),
                      colors=c("green","blue"), bin.width=0.5,areas=c("L","M"),xlim.low = num.bins-(1/bin.width)*trunc.year)
#mtext("A",side=1,2.5,adj=-.25,cex=2)
mtext("(a)",side=3,line=-1,adj=-.28,cex=1.2)

# mountain dispersal and cladogenesis rates
unit_label="(#/lineage/Ma)"

num.bins=length(binned.CI.My$M$dispersal.into.rates.CI)
bin.width=0.5
trunc.year=3
plot_binned_multi_var(binned.CI.My,vars=c("insitu.rates.CI","dispersal.into.rates.CI"),
                      title=NULL,yaxis.title = paste0("Rate ",unit_label),
                      legend.override=c("Mountain cladogenesis","Mountain colonization"),legend.pos=c(12.5,2),
                      colors=c("blue"), bin.width=0.5,areas=c("M"),xlim.low = num.bins-(1/bin.width)*trunc.year,ylim.pad=.2)


#mtext("B",side=1,2.5,adj=-.25,cex=2)
mtext("(b)",side=3,line=-1,adj=-.28,cex=1.2)

dev.off()
par(opar)

##########
## version with three panels
##########
# ## simple grid
# m1=matrix(c(1,0,2,3),ncol=2)
# 
# pdf(file=file.path(plot_dir,"paneled_binned_plots.pdf"),width=8,height=7)
# 
# layout(m1,heights=c(1,1,.9,1.1))
# #par(mar=c(4,4,1,.5)) #1=bottom, 2=left, 3=top, 4=right
# 
# # combined events plot
# num.bins=length(binned.CI$M$dispersal.into.time.bins.cumsum.CI)
# bin.width=0.5
# trunc.year=8
# 
# par(mar=c(4,4,1,.5),oma=c(1,1,1,1))
# plot_binned_multi_var(binned.CI,vars=c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),
#                       title=NULL,yaxis.title = "Cumulative events",legend.override=c("Lowland cladogenesis","Mountain cladogenesis","Lowland dispersal","Mountain dispersal"),
#                       colors=c("green","blue"), bin.width=0.5,areas=c("L","M"),xlim.low = num.bins-(1/bin.width)*trunc.year)
# mtext("A",side=1,2.5,adj=-.25,cex=2)
# 
# # indivdual dispersal and cladogenesis rates
# par(mar=c(2,4,1,.5))
# #unit_label=paste0("(#/",expression(frac(1,2)),"Ma)")
# unit_label="(#/0.5Ma)"
# plot_binned_single_var(binned.CI,var="insitu.rates.CI",title=NULL,
#                        yaxis.title=paste0("Cladogenesis rate ",unit_label),
#                        legend.override=c("Lowland","Mountain"),
#                        colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),show.x.axis=FALSE,xlim.low = 6)
# mtext("B",side=1,0.5,adj=-.25,cex=2)
# 
# num.bins=length(binned.CI$M$dispersal.into.rates.CI)
# bin.width=0.5
# trunc.year=6
# 
# par(mar=c(4,4,0,.5))
# plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title=NULL,
#                        yaxis.title=paste0("Dispersal rate ",unit_label),
#                        legend.override=c("Lowland","Mountain"),
#                        colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),
#                        ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year)
# mtext("C",side=1,2.5,adj=-.25,cex=2)
# dev.off()
# par(opar)
# 
# 
# trunc.year=3
# plot_binned_multi_var(binned.CI,vars=c("insitu.rates.CI","dispersal.into.rates.CI"),
#                        title=NULL,yaxis.title = "Cumulative events",legend.override=c("Mountain cladogenesis","Mountain dispersal"),
#                        colors=c("blue"), bin.width=0.5,areas=c("M"),xlim.low = num.bins-(1/bin.width)*trunc.year)

##################################
## Fig 4: median niche distance ##
##################################
## from slice_BSMs_OUM
slice_dir=file.path(mv_dir,"sortadate50_merged_TTRtips/slice_output")
slice_results_long=read.csv(file.path(slice_dir,"slice_summary_long_50BSMs.csv"))

pdf(file=file.path(slice_dir,"nichedist_hist_50BSMs.pdf"),width=6,height=5)
ggplot(data=slice_results_long) + geom_histogram(aes(x=med.dist,fill=event_type),binwidth = 0.005) + 
  xlab("Median niche distance") +  ylab("Count") + theme_classic(base_size = 16) + 
  theme(legend.position = c(0.7,0.7)) + scale_fill_discrete(name="Node source",labels=c("Cladogenesis","Colonization"))
dev.off()

################
## Fig 5: DTT ##
################
library(geiger)
source(file.path(script_dir,"biogeography_github/niche_evolution/dtt_abstime.R"))
phy=ape::read.tree(file.path(TTR_dir,"sortadate50_TTRtips.tre"))
phy$tip.label = gsub("Veronica","V.",phy$tip.label)
traits_bands_back=read.csv(file.path(TTR_dir,"TTR9_traits_bands_back.csv"))
row.names(traits_bands_back)=traits_bands_back$species
traitcols=3:11 ## check this

## plot
pdf(file.path(ms_plot_dir,"dtt_plots.pdf"),height=4,width=7)
layout(matrix(1:2,nrow=1))
par(mar=c(4.1,4.2,2,1),mgp=c(2.5,1,0))
## subclade disparity
dtt_TTR9<-dtt.abstime(phy=phy, data=select(traits_bands_back,traitcols), nsim=100, plot=TRUE)
box()
mtext("(a)",side=3,line=-1,adj=-.32,cex=1.2)

## ancestral time slice disparity
BSM_out=readRDS(file.path(slice_dir,"slice_disparity_MLdistdf_output_50BSMs.Rdata"))
mindisp=min(unlist(lapply(BSM_out,function(x) min(x$dist_df$disp))))
maxdisp=max(unlist(lapply(BSM_out,function(x) max(x$dist_df$disp))))
maxtime=max(unlist(lapply(BSM_out,function(x) max(x$dist_df$time_bp))))

plot(BSM_out[[1]]$dist_df[BSM_out[[1]]$dist_df$area=="L","time_bp"],
     BSM_out[[1]]$dist_df[BSM_out[[1]]$dist_df$area=="L","disp"],col="white",
     ylim=c(mindisp,maxdisp),xlim=rev(c(0,maxtime)),xlab="Time (Ma)",ylab="Disparity")
for(i in names(BSM_out)) lines(BSM_out[[i]]$dist_df[BSM_out[[i]]$dist_df$area=="L","time_bp"],
                               BSM_out[[i]]$dist_df[BSM_out[[i]]$dist_df$area=="L","disp"],
                               xlim=rev(c(0,maxtime)),
                               col=adjustcolor("green",alpha.f=0.75))
for(i in names(BSM_out)) lines(BSM_out[[i]]$dist_df[BSM_out[[i]]$dist_df$area=="M","time_bp"],
                               BSM_out[[i]]$dist_df[BSM_out[[i]]$dist_df$area=="M","disp"],
                               xlim=rev(c(0,maxtime)),
                               col=adjustcolor("blue",alpha.f=0.5))
legend('bottomright',c("Lowland","Mountain"),col=c("green","blue"),lty=1,lwd=1.5,cex=.8,bty='n')
mtext("(b)",side=3,line=-1,adj=-.32,cex=1.2)
dev.off()


########################
## Fig 6: niche-range ##
########################
## load ARC results after running on cluster (with Dan Warren's improved ARC)
arc.no.M=readRDS(file.path(niche_dir,"arc.seed123_niche_overlap_phy_M.Rdata"))
arc.ro.M=readRDS(file.path(niche_dir,"arc.seed123_range_overlap_phy_M.Rdata"))

pdf(file.path(ms_plot_dir,"niche_range_ARC_sims_supp.pdf"),height=3.5,width=7)
layout(matrix(1:2,nrow=1))
par(mar=c(4.1,4.1,1,.5)) #1=bottom, 2=left, 3=top, 4=right
## age-range correlation with spatial range overlap and null simulations
plot(arc.ro.M$age.range.correlation,xlab="Age (mya)",ylab="Range overlap")
apply(t(arc.ro.M$MonteCarlo.replicates), 1, abline, lwd = 0.2, col = "grey50")
abline(arc.ro.M$linear.regression$coefficients,col="red")
mtext("(a)",side=3,line=-1,adj=-.3,cex=1.2)

## age-range correlation with TTR niche overlap and null simulations
plot(arc.no.M$age.range.correlation,xlab="Age (mya)",ylab="Niche overlap")
apply(t(arc.no.M$MonteCarlo.replicates), 1, abline, lwd = 0.2, col = "grey50")
abline(arc.no.M$linear.regression$coefficients,col="red")
mtext("(b)",side=3,line=-1,adj=-.3,cex=1.2)

dev.off()


## null median comparisons
pdf(file.path(plot_dir,paste0("niche_range_null_medians_hist_",Sys.Date(),".pdf")),height=3.5,width=7)
layout(matrix(1:2,nrow=1))
par(mar=c(4.5,4.1,2,.5)) #1=bottom, 2=left, 3=top, 4=right

## histogram of median node-averaged range overlap, observed vs null simulations
arc.ro.medians=apply(arc.ro.M$sim.overlaps[,-1],2,median)
obs.med.r=median(arc.ro.M$age.range.correlation[,2]) #[1] 0.07311673
hist(arc.ro.medians,xlim=c(min(arc.ro.medians),round(obs.med.r+.005,2)),main=NULL,
     xlab="Median range overlap",ylab="Count",cex.axis=1,cex.lab=1.2)
abline(v=median(arc.ro.M$age.range.correlation[,2]),lwd=2)
box()
#mtext("A",side=1,line=3.5,adj=-.2,cex=2)
mtext("(a)",side=3,line=.5,adj=-.25,cex=1.2)

## histogram of median node-averaged niche overlap, observed vs null simulations
par(mar=c(4.5,0.6,2,4))
arc.no.medians=apply(arc.no.M$sim.overlaps[,-1],2,median)
obs.med.n=median(arc.no.M$age.range.correlation[,2]) #[1] 0.368052
hist(arc.no.medians,xlim=c(min(arc.no.medians),round(obs.med.n+.01,2)),main=NULL,
     xlab="Median niche overlap",ylab=NULL,yaxt='n',cex.axis=1,cex.lab=1.2)
abline(v=median(arc.no.M$age.range.correlation[,2]),lwd=2)
box()
mtext("(b)",side=3,line=.5,adj=-.05,cex=1.2)
dev.off()

####################
## specimen table ##
#################### 
## from collate_specimen_table.R

## includes some nice extra info, but not updated OLD numbers
sequenced_info_metadata_old=read.csv(file.path(hyb_dir,"sequenced_libs_metadata_20200904.csv"))
sequenced_info_metadata_old[sequenced_info_metadata_old$Species_simple=="V. plano-petiolata",
                            "Species_simple"]="V. planopetiolata"
## has updated OLD numbers but not extra info
sequenced_info_metadata_raw=read.csv(file.path(data_dir,"samples_extraction/sequenced_libs_metadata_full_20210301.csv"))
## has habitat info for species
alt_lat=read.csv(file.path(dist_dir,"alt_lat_bands.csv"))
alt_lat$Species_simple=gsub("Veronica_","V. ",alt_lat$species)
## has occurrence data
sdm_stats=read.csv(file.path(TTR_dir,"Veronica_sdm_stats.csv"))
sdm_stats$Species_simple=gsub("Veronica_","V. ",sdm_stats$sp)
sdm_stats$Occurrences=sdm_stats$d.tpos+sdm_stats$d.fneg

## create combined table
new_acc_table=sequenced_info_metadata_old %>% 
  dplyr::select(Sample,Species,Species_simple,natural_group,chromosome.no,Year.collected) %>%
  left_join(sequenced_info_metadata_raw[,c("Voucher","Sample.lib.Anne","AmJBot_subset")],
            by=c("Sample"="Sample.lib.Anne")) %>%
  left_join(alt_lat[,c("Species_simple","category","merged_category")],by="Species_simple") %>%
  left_join(sdm_stats[,c("Species_simple","Occurrences","d.tpos","d.fpos","d.tneg","d.fneg")],by="Species_simple")

## island species
new_acc_table %>% dplyr::filter(is.na(category)) %>% dplyr::select(Species_simple) 
## species with too few occurrences for TTR
new_acc_table %>% dplyr::filter(is.na(Occurrences),!is.na(category)) %>% dplyr::select(Species_simple) 


new_acc_table_edits=new_acc_table
new_acc_table_edits[is.na(new_acc_table_edits$category),"category"]="Island"
new_acc_table_edits[is.na(new_acc_table_edits$merged_category),"merged_category"]="Island"
new_acc_table_edits[new_acc_table_edits$Species_simple %in% c("V. chamaedrys","V. perfoliata"),"category"]="Outgroup"
new_acc_table_edits[new_acc_table_edits$Species_simple %in% c("V. chamaedrys","V. perfoliata"),"merged_category"]="Outgroup"

write.csv(new_acc_table_edits,file.path(dist_dir,"Species_accessions_occ.csv"))

ms_acc_table = new_acc_table_edits %>% 
  dplyr::select(Species_simple,Voucher,merged_category,category,Occurrences,
                d.tpos,d.fpos,d.tneg,d.fneg,
                natural_group,chromosome.no,Year.collected,AmJBot_subset) %>% 
  dplyr::rename(Species=Species_simple,Habitat_two_band=merged_category,
                Habitat_three_band=category,Subclade=natural_group,
                Chromosome_no=chromosome.no,Year_collected=Year.collected,
                Incl_in_Thomas2021=AmJBot_subset,TTR.tpos=d.tpos,TTR.fpos=d.fpos,
                TTR.tneg=d.tneg,TTR.fneg=d.fneg) %>%
  dplyr::arrange(Species)

write.csv(ms_acc_table,file.path(plot_dir,"Species_metadata.csv"),row.names = FALSE)


###################
## supplementary ##
###################
## Fig S0: PP beast tree
library(rBt)
source('~/Cambridge/Research/scripts/biogeography_github/util/plot.phylo.HPD.fixed.R')
phy_nexus_path=file.path(bgb_dir,"rank50_strictEstimate.species.ann.tre")

phy_nexus_rbt=rBt::read.annot.beast(phy_nexus_path)
tip.idx=which(phy_nexus_rbt$tip.label=="chamaedrys")
node.idx=Ntip(phy_nexus_rbt)+1

phy_red=ape::drop.tip(phy_nexus_rbt,"chamaedrys")
phy_red$posterior=phy_red$posterior[-1]
phy_red$metadata=phy_red$metadata[-c(tip.idx, node.idx),]
phy_red$metadata$node=1:nrow(phy_red$metadata)


libs_info=read.csv(file.path(hyb_dir,"libs_info_detailed.csv"),stringsAsFactors=FALSE)
libs_info[which(libs_info$Species=="V. vernicosa-2926"),"Species"]="V. vernicosa"
libs_info[which(libs_info$Species=="V. hookeriana-BDM12"),"Species"]="V. hookeriana"
libs_info$Species=gsub("plano-petiolata", "planopetiolata",libs_info$Species)

hebe_groups = read.csv(file.path(hyb_dir,"hebe_groups_simple_colors.csv"),stringsAsFactors = FALSE)
hebe_colors = unique(hebe_groups$color)
legend_specs = read.csv(file.path(tr_plot_dir,"combined_legend2.csv"),stringsAsFactors = FALSE)


phy_red$tip.label=paste("V.",phy_red$tip.label)
phy_red$tip.label=gsub("plano-petiolata", "planopetiolata",phy_red$tip.label)

tips = dplyr::left_join(data.frame("tip"=phy_red$tip.label),libs_info[,c("Species","Sample","natural_group","hebe_group","hebe_group_detail","alpine_status")],by=c("tip"="Species"))
tips=dplyr::left_join(tips,hebe_groups)

pdf(file.path(plot_dir,"pp_bar_tree_color2.pdf"),height=10,width=8)
par(mar=c(4,1,1,1),xpd=FALSE)
plot.phylo.HPD.fixed(phy_red,bar.width=.2,bar.col=adjustcolor("purple",alpha.f = 0.5),
                     vline=FALSE,border=NA,cex=.6,label.offset=0.15)
tiplabels(pch=19,col=hebe_colors[tips$code],adj=.6,cex=.5)
node.font=ifelse(phy_red$posterior>=0.9,2,1)
node.color=ifelse(phy_red$posterior>=0.9,"blue","black")
node.color="black"
node.text=ifelse(phy_red$posterior>=0.8,round(phy_red$posterior,2),"")
nodelabels(text=node.text, cex=0.55, frame="none",adj=c(1.1,1.2),col=node.color,font=node.font)

legend(x=0,y=20,legend=legend_specs$text,col=legend_specs$color,
       text.font=legend_specs$font,#text.col=legend_specs$text.color,
       pch=legend_specs$pch,cex=.8)
mtext(text="Time (Ma)", side=1, line=2, cex=1)
dev.off()

#############################
## Fig S2: ASL binned plots #
#############################
## three bands: 100 trees ##
############################
plot2_dir=file.path(BGB_dir,"biogeobears_output_100trees/sortadate50_ASL/binned_plots/plots_2.0")
#dir.create(plot2_dir)
load(file.path(BGB_dir,"biogeobears_output_100trees/sortadate50_ASL/sortadate50_ASL_binned_output_multi.Rdata"))
binned.CI=output.list[[1]]

## change Os to NA for mountains in binned.CI before 3.5 mya
binned.CI$S$insitu.rates.CI=pad.NA(binned.CI$S$insitu.rates.CI,0.5)
binned.CI$S$dispersal.into.rates.CI = pad.NA(binned.CI$S$dispersal.into.rates.CI,0.5)
binned.CI$S$insitu.time.bins.cumsum.CI=pad.NA(binned.CI$S$insitu.time.bins.cumsum.CI,0.5,4)
binned.CI$S$dispersal.into.time.bins.cumsum.CI = pad.NA(binned.CI$S$dispersal.into.time.bins.cumsum.CI,0.5,4)

binned.CI$A$insitu.rates.CI=pad.NA(binned.CI$A$insitu.rates.CI,0.5,1.5)
binned.CI$A$dispersal.into.rates.CI = pad.NA(binned.CI$A$dispersal.into.rates.CI,0.5,1.5)
binned.CI$A$insitu.time.bins.cumsum.CI=pad.NA(binned.CI$A$insitu.time.bins.cumsum.CI,0.5,2)
binned.CI$A$dispersal.into.time.bins.cumsum.CI = pad.NA(binned.CI$A$dispersal.into.time.bins.cumsum.CI,0.5,2)

binned.CI.My=binned.CI
binned.CI.My$M$insitu.rates.CI=lapply(binned.CI.My$M$insitu.rates.CI,function(x) x/2)
binned.CI.My$M$dispersal.into.rates.CI=lapply(binned.CI.My$M$dispersal.into.rates.CI,function(x) x/2)

###############
## panelling ##
###############
# save original par settings
opar=par(no.readonly = TRUE) # restore with par(opar)

## simple grid
m2=matrix(c(1:4),ncol=2)

pdf(file=file.path(ms_plot_dir,"drafting/paneled_binned_plots_abcd.pdf"),width=8,height=7)

layout(m2)
par(mar=c(4,4,1,.5)) #1=bottom, 2=left, 3=top, 4=right

# A top left: in situ events
par(mar=c(2,4,1,.5))
bin.width=.5
num.bins=length(binned.CI$A$insitu.time.bins.cumsum.CI)
trunc.year=7
plot_binned_single_var(binned.CI,var="insitu.time.bins.cumsum.CI",title=NULL,yaxis.title="Cladogenesis events",
                       legend.override=c("Lowland","Subalpine","Alpine"), show.x.axis=FALSE,
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),xlim.low = num.bins-(1/bin.width)*trunc.year)
mtext("(a)",side=3,-1,adj=-.2,cex=1.2)

# c bottom left: dispersal events
par(mar=c(4,4,0,.5))
bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.time.bins.cumsum.CI)
trunc.year=7
plot_binned_single_var(binned.CI,var="dispersal.into.time.bins.cumsum.CI",title=NULL,yaxis.title="Colonization events",
                       legend.override=c("Lowland","Subalpine","Alpine"),colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),
                       xlim.low = num.bins-(1/bin.width)*trunc.year,lty=2)

mtext("(c)",side=3,-1,adj=-.2,cex=1.2)

# b top right: in situ rates
par(mar=c(2,4,1,.5))
bin.width=.5
num.bins=length(binned.CI$A$insitu.rates.CI)
trunc.year=6
plot_binned_single_var(binned.CI.My,var="insitu.rates.CI",title=NULL,yaxis.title="Cladogenesis rate (#/Ma)",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"), legend=FALSE, show.x.axis=FALSE,
                       xlim.low = num.bins-(1/bin.width)*trunc.year)
mtext("(b)",side=3,-1,adj=-.2,cex=1.2)

# d bottom right: dispersal rates
par(mar=c(4,4,0,.5))
bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.rates.CI)
trunc.year=6
plot_binned_single_var(binned.CI.My,var="dispersal.into.rates.CI",title=NULL,yaxis.title="Colonization rate (#/Ma)",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),legend=FALSE,
                       ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year,lty=2)
mtext("(d)",side=3,-1,adj=-.2,cex=1.2)

dev.off()

############
## rpanda ##
############
library(RPANDA)

## paths
source('~/Cambridge/Research/scripts/biogeography_github/rpanda/plot_fit_env_custom_rpanda.R')
data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
bgb_dir=file.path(data_dir,"biogeobears")
paleo_dir=file.path(data_dir,"NZ_geo/NZ_paleoelevation")

res.env.lin = readRDS(file.path(paleo_dir,"rpanda_elev/lin_elev_model.Rdata"))
### elevation data
load(file.path(paleo_dir,"/NZalps_elev.Rdata"))

pdf(file.path(paleo_dir,"rpanda_elev/elevation_spec_rate_suppfig.pdf"),width=7,height=4)
layout(matrix(1:2,ncol=2))

plot(elevation_data[,1],elevation_data[,2],xlim=rev(range(elevation_data[,1])),xlab="Time (Ma)", ylab="Elevation (masl)",cex.lab=1.2,cex.axis=1.2)
lines(smooth_spline_fit,lwd=2,xlim=rev(range(elevation_data[,1])))
mtext("(a)",side=3,line=.75,adj=-.25,cex=1.2)

linear_elev=readRDS(file.path(paleo_dir,"linear_interpolation_elevation.Rdata"))
phy=read.tree(file.path(bgb_dir,"sortadate50_tree.tre"))
tot_time<-max(node.age(phy)$ages)

plot_fit_env_custom(res.env.lin, linear_elev,tot_time,spec.rate.only=TRUE)
mtext("(b)",side=3,line=.75,adj=-.25,cex=1.2)

dev.off()
