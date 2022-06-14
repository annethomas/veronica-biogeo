data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
bgb_dir=file.path(data_dir,"biogeobears")

source("C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/biogeobears/binning/bin_plot_BSMs.R")

#####################################
## two bands (merged and lowlands) ##
#####################################

plot_dir=file.path(bgb_dir,"biogeobears_output_100trees/sortadate50_merged/binned_plots/bluegreen/")
load(file.path(bgb_dir,"biogeobears_output_100trees/sortadate50_merged/sortadate50_merged_binned_output_multi.Rdata"))
binned.CI=output.list[[1]]

## change Os to NA for mountains in binned.CI before 3.5 mya (prevent extra line)
binned.CI$M$insitu.rates.CI=pad.NA(binned.CI$M$insitu.rates.CI,bin.width=0.5)
binned.CI$M$dispersal.into.rates.CI = pad.NA(binned.CI$M$dispersal.into.rates.CI,bin.width=0.5)
binned.CI$M$insitu.time.bins.cumsum.CI=pad.NA(binned.CI$M$insitu.time.bins.cumsum.CI,bin.width=0.5,max.year=4)
binned.CI$M$dispersal.into.time.bins.cumsum.CI = pad.NA(binned.CI$M$dispersal.into.time.bins.cumsum.CI,bin.width=0.5,max.year=4)

## plot cladogenesis and dispersal rates separately
pdf(file=file.path(plot_dir,"sortadate50_merged_insitu_v_dispersal_rates.pdf"),width=6,height=6)
plot_binned_single_var(binned.CI,var="insitu.rates.CI",title=NULL,yaxis.title="In situ cladogenesis rate",
                       colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),xlim.low = 6)
num.bins=length(binned.CI$M$dispersal.into.rates.CI)
bin.width=0.5
trunc.year=4
plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title=NULL,yaxis.title="Dispersal rate",
                       colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),
                       ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year)

dev.off()

## plot clado and dispersal events together
pdf(file=file.path(plot_dir,"sortadate50_merged_insitu_dispersal_events.pdf"),width=6,height=6)
plot_binned_multi_var(binned.CI,vars=c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),
                      title=NULL,yaxis.title = "Cumulative events",var.titles=c("cladogenesis","dispersal"),
                      colors=c("green","blue"), bin.width=0.5,areas=c("L","M"))
dev.off()

###############
## panelling ##
###############
# save original par settings
opar=par(no.readonly = TRUE) # restore with par(opar)

## simple grid
m1=matrix(c(1,0,2,3),ncol=2)

pdf(file=file.path(plot_dir,"paneled_binned_plots.pdf"),width=8,height=7)

layout(m1,heights=c(1,1,.9,1.1))
#par(mar=c(4,4,1,.5)) #1=bottom, 2=left, 3=top, 4=right

# combined events plot
num.bins=length(binned.CI$M$dispersal.into.time.bins.cumsum.CI)
bin.width=0.5
trunc.year=8

par(mar=c(4,4,1,.5),oma=c(1,1,1,1))
plot_binned_multi_var(binned.CI,vars=c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),
                      title=NULL,yaxis.title = "Cumulative events",legend.override=c("Lowland cladogenesis","Mountain cladogenesis","Lowland dispersal","Mountain dispersal"),
                      colors=c("green","blue"), bin.width=0.5,areas=c("L","M"),xlim.low = num.bins-(1/bin.width)*trunc.year)
mtext("A",side=1,2.5,adj=-.25,cex=2)

# indivdual dispersal and cladogenesis rates
par(mar=c(2,4,1,.5))
#unit_label=paste0("(#/",expression(frac(1,2)),"Ma)")
unit_label="(#/0.5Ma)"
plot_binned_single_var(binned.CI,var="insitu.rates.CI",title=NULL,
                       yaxis.title=paste0("Cladogenesis rate ",unit_label),
                       legend.override=c("Lowland","Mountain"),
                       colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),show.x.axis=FALSE,xlim.low = 6)
mtext("B",side=1,0.5,adj=-.25,cex=2)

num.bins=length(binned.CI$M$dispersal.into.rates.CI)
bin.width=0.5
trunc.year=6

par(mar=c(4,4,0,.5))
plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title=NULL,
                       yaxis.title=paste0("Dispersal rate ",unit_label),
                       legend.override=c("Lowland","Mountain"),
                       colors=c("green","blue"),bin.width=0.5,areas=c("L","M"),
                       ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year)
mtext("C",side=1,2.5,adj=-.25,cex=2)
dev.off()

## different sizes (didn't use)
m2=matrix(c(1,1,2,3),ncol=2)
layout(m2, widths=c(2,1,1))
layout(m2)

m3=matrix(c(1,2,1,3),ncol=2)
layout(m3)

############################
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

## plot cladogenesis and dispersal rates separately
pdf(file=file.path(plot2_dir,"sortadate50_ASL_insitu_v_dispersal_rates2.pdf"),width=6,height=6)
plot_binned_single_var(binned.CI,var="insitu.rates.CI",title=NULL,yaxis.title="In situ cladogenesis rate",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),xlim.low = num.bins-(1/bin.width)*trunc.year)
bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.rates.CI)
trunc.year=6
plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title=NULL,yaxis.title="Dispersal rate",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),
                       ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year)

dev.off()

bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.time.bins.cumsum.CI)
trunc.year=7
pdf(file=file.path(plot2_dir,"sortadate50_ASL_insitu_v_dispersal_events2.pdf"),width=6,height=6)
plot_binned_single_var(binned.CI,var="insitu.time.bins.cumsum.CI",title=NULL,yaxis.title="Cumulative cladogenesis events",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),xlim.low = num.bins-(1/bin.width)*trunc.year)

plot_binned_single_var(binned.CI,var="dispersal.into.time.bins.cumsum.CI",title=NULL,yaxis.title="Cumulative dispersal events",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),
                       xlim.low = num.bins-(1/bin.width)*trunc.year)

dev.off()

layout(matrix(1,2,2))
pdf(file=file.path(plot2_dir,"sortadate50_ASL_insitu_v_dispersal_rates_panel.pdf"),width=12,height=12)
dev.off()

###############
## panelling ##
###############
# save original par settings
opar=par(no.readonly = TRUE) # restore with par(opar)

## simple grid
m1=matrix(c(1:4),ncol=2)

pdf(file=file.path(plot2_dir,"paneled_binned_plots.pdf"),width=8,height=7)

layout(m1)
par(mar=c(4,4,1,.5)) #1=bottom, 2=left, 3=top, 4=right

# A top: in situ events
par(mar=c(2,4,1,.5))
bin.width=.5
num.bins=length(binned.CI$A$insitu.time.bins.cumsum.CI)
trunc.year=7
plot_binned_single_var(binned.CI,var="insitu.time.bins.cumsum.CI",title=NULL,yaxis.title="Cladogenesis events",
                       legend.override=c("Lowland","Subalpine","Alpine"), show.x.axis=FALSE,
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),xlim.low = num.bins-(1/bin.width)*trunc.year)

# A bottom: dispersal events
par(mar=c(4,4,0,.5))
bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.time.bins.cumsum.CI)
trunc.year=7
plot_binned_single_var(binned.CI,var="dispersal.into.time.bins.cumsum.CI",title=NULL,yaxis.title="Dispersal events",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),legend=FALSE,
                       xlim.low = num.bins-(1/bin.width)*trunc.year)

mtext("A",side=1,2.5,adj=-.1,cex=2)

# B top: in situ rates
par(mar=c(2,4,1,.5))
bin.width=.5
num.bins=length(binned.CI$A$insitu.rates.CI)
trunc.year=6
plot_binned_single_var(binned.CI,var="insitu.rates.CI",title=NULL,yaxis.title="Cladogenesis rate (#/0.5Ma)",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"), legend=FALSE, show.x.axis=FALSE,
                       xlim.low = num.bins-(1/bin.width)*trunc.year)

# B bottom: dispersal rates
par(mar=c(4,4,0,.5))
bin.width=.5
num.bins=length(binned.CI$A$dispersal.into.rates.CI)
trunc.year=6
plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title=NULL,yaxis.title="Dispersal rate (#/0.5Ma)",
                       colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"),legend=FALSE,
                       ylim=get.ylim(binned.CI,"insitu.rates.CI"),xlim.low = num.bins-(1/bin.width)*trunc.year)
mtext("B",side=1,2.5,adj=-.1,cex=2)

dev.off()
#####




#########################################
### obsolete
## three bands: MCC

load(file.path(bgb_dir,"sortadate50_time_expertBands/test_sortadate50_binned_output_multi.Rdata"))
load(file.path(bgb_dir,"biogeobears_output_100trees/sortadate50_ASL/sortadate50_binned_output_multi.Rdata"))

binned.CI=output.list[[1]]

pdf(file=file.path(plot2_dir,paste0("sortadate50_insitu_dispersal_events.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI,vars=c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),
                      title="In situ and dispersal events",yaxis.title="Cumulative events",
                      colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"))
dev.off()

pdf(file=file.path(plot2_dir,paste0("sortadate50_insitu_dispersal_rates.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI,vars=c("insitu.rates.CI","dispersal.into.rates.CI"),
                      title="In situ and dispersal rates",yaxis.title="Rolling per capita rate",
                      colors=c("red","blue","green"),bin.width=0.5,areas=c("L","S","A"))
dev.off()

vars=c("insitu.rates.CI","dispersal.into.rates.CI")

## split up NI and SI
load(file.path(bgb_dir,"sortadate50_time_expertBands_mergedNS2/sortadate50_mergedNS_binned_output_multi.Rdata"))
binned.CI.NI=output.list[[1]][1:2]
binned.CI.SI=output.list[[1]][3:4]

pdf(file=file.path(bgb_dir,paste0("sortadate50_NI_insitu_dispersal_events.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI.NI,c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),"NI In situ and dispersal events","Cumulative events",c('#0000FF','#00CD00'),1,areas=c("A","B"))
dev.off()

pdf(file=file.path(bgb_dir,paste0("sortadate50_SI_insitu_dispersal_events.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI.SI,c("insitu.time.bins.cumsum.CI","dispersal.into.time.bins.cumsum.CI"),"SI In situ and dispersal events","Cumulative events",c('#FFFF00','#FF0000'),1,areas=c("C","D"))
dev.off()

pdf(file=file.path(bgb_dir,paste0("sortadate50_NI_insitu_dispersal_rates.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI.NI,c("insitu.rates.CI","dispersal.into.rates.CI"),"In situ and dispersal rates","NI Rolling per capita rate",c('#0000FF','#00CD00'),1,areas=c("A","B"))
dev.off()

pdf(file=file.path(bgb_dir,paste0("sortadate50_SI_insitu_dispersal_rates.pdf")),width=6,height=6)
plot_binned_multi_var(binned.CI.SI,c("insitu.rates.CI","dispersal.into.rates.CI"),"In situ and dispersal rates","SI Rolling per capita rate",c('#FFFF00','#FF0000'),1,areas=c("C","D"))
dev.off()


## playing with colors
library(scales)
show_col("blue")
show_col(darken("blue"))
show_col(lighten("blue"))
show_col(lighten("blue",0.7))

## get a blend between blue and green for merged plot
pal=colorRampPalette(c("blue","green"))
pal(3)
#[1] "#0000FF" "#007F7F" "#00FF00"
show_col("#007F7F")


