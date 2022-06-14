##############################################
## prepare and process age-range correlations 

library(ape)
library(phyloclim)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data/"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github/"
dist_dir=file.path(data_dir,"distribution_data")
bgb_dir=file.path(data_dir,"biogeobears")
TTR_dir=file.path(data_dir,"TTR")
niche_dir=file.path(data_dir,"biogeography_paper/niche_overlap")

phy=read.tree(file.path(bgb_dir,"sortadate50_TTRtips.tre"))
alt_lat=read.csv(file.path(dist_dir,"alt_lat_bands.csv"))

## load functions (incl clean_matrix)
source(file.path(script_dir,"niche_overlap/overlap_util_functions.R"))

#######################
## load overlap data ##
#######################
# from nicheOverlap_spaa_phyloclim.R
no.list=readRDS(file=file.path(niche_dir,"TTRoverlap_phyloclim_spaa_01_05.Rdata"))
# from spatialOverlap_phyloclim.R
#ro.list=readRDS(file=file.path(niche_dir,"range_overlap_phyloclim.Rdata"))
ro.list=readRDS(file=file.path(niche_dir,"range_overlap_phyloclim_wgbif.Rdata"))

jacc.scores=as.matrix(read.csv(file.path(niche_dir,"range_overlap_jaccard.csv"),row.names=1))

#############################################
## prune data/phylogeny to specific subset ##
#############################################
## mountain only (M or ML)

tips=data.frame(tips=phy$tip.label)
setdiff(tips$tips,alt_lat$species)
tips=left_join(tips,alt_lat[,c("species","merged_category")],by=c("tips"="species"))
phy.M=drop.tip(phy,tips[!tips$merged_category %in% c("M","ML"),"tips"])
plot(phy)
plot(phy.M)

ro.full=ro.list[["ro.full"]]
ro.M=clean_matrix(ro.full,phy.M$tip.label)
nrow(ro.M)
length(phy.M$tip.label)

no.full=no.list[["TTR.no.full01"]] ## before, this was just "no.full" from the ro.full object, neither of which seemed to be correct
no.M=clean_matrix(no.full,phy.M$tip.label)
nrow(no.M)
length(phy.M$tip.label)

## save for running ARC on cluster (takes a long time, several minutes per MC simulation)
saveRDS(list(phy=phy.M,overlap=no.M),file=file.path(niche_dir,"niche_overlap_phy_M.Rdata"))
# saveRDS(list(phy=phy.M,overlap=ro.M),file=file.path(niche_dir,"range_overlap_phy_M.Rdata"))
saveRDS(list(phy=phy.M,overlap=ro.M),file=file.path(niche_dir,"range_overlap_phy_wgbif_M.Rdata"))


## load ARC results after running on cluster (with Dan Warren's improved ARC)
arc.no.M=readRDS(file.path(niche_dir,"arc.seed123_niche_overlap_phy_M.Rdata"))
#arc.ro.M=readRDS(file.path(niche_dir,"arc.seed123_range_overlap_phy_M.Rdata"))
arc.ro.M=readRDS(file.path(niche_dir,"arc.seed123_range_overlap_phy_wgbif_M.Rdata"))

##############
## plot ARC ##
##############
## age-range correlation with spatial range overlap and null simulations
plot(arc.ro.M$age.range.correlation,ylab="range overlap")
apply(t(arc.ro.M$MonteCarlo.replicates), 1, abline, lwd = 0.2, col = "grey50")
abline(arc.ro.M$linear.regression$coefficients,col="red")

## age-range correlation with TTR niche overlap and null simulations
plot(arc.no.M$age.range.correlation,ylab="niche overlap")
apply(t(arc.no.M$MonteCarlo.replicates), 1, abline, lwd = 0.2, col = "grey50")
abline(arc.no.M$linear.regression$coefficients,col="red")

## histogram of median node-averaged range overlap, observed vs null simulations
arc.ro.medians=apply(arc.ro.M$sim.overlaps[,-1],2,median)
#hist(arc.ro.medians)
obs.med.r=median(arc.ro.M$age.range.correlation[,2]) #[1] 0.07311673
hist(arc.ro.medians,xlim=c(min(arc.ro.medians),round(obs.med.r+.005,2)),main="Median node-avg range overlap",xlab="Null range overlap medians")
abline(v=median(arc.ro.M$age.range.correlation[,2]))
mtext("Observed range overlap median",adj=1)

## histogram of median node-averaged niche overlap, observed vs null simulations
arc.no.medians=apply(arc.no.M$sim.overlaps[,-1],2,median)
#hist(arc.no.medians)
obs.med.n=median(arc.no.M$age.range.correlation[,2]) #[1] 0.368052
hist(arc.no.medians,xlim=c(min(arc.no.medians),round(obs.med.n+.01,2)),main="Median node-avg niche overlap",xlab="Null niche overlap medians")
abline(v=median(arc.no.M$age.range.correlation[,2]))
mtext("Observed niche overlap median",adj=1)

####################
## range vs niche ##
####################
## range vs niche regression/correlation with null simulations
obs.regr= lm(arc.no.M$age.range.correlation[,"overlap"] ~ arc.ro.M$age.range.correlation[,"overlap"])

sims.r=arc.ro.M$sim.overlaps
sims.n=arc.no.M$sim.overlaps

coeff=matrix(nrow=100,ncol=2)
for(n in 2:ncol(sims.r)){
  print(n)
  i=n-1
  fit=lm(sims.n[,n] ~ sims.r[,n])
  coeff[i,]=fit$coefficients
}

## plot niche-range correlation with null simulations
plot(arc.ro.M$age.range.correlation[,"overlap"],arc.no.M$age.range.correlation[,"overlap"],
     xlab="range overlap",ylab="niche overlap")
apply(coeff, 1, abline, lwd = 0.2, col = "grey50")
abline(obs.regr$coefficients,col="red")

#######################
## sister pairs only ##
#######################
library(dispRity)
sch.mat=no.list$sch.mat01

## get sisters with habitat categories
sisters=find_sisters(phy,alt_lat,"merged_category")

## get overlap info and age for sister pairs
sisters$range_overlap=NA
sisters$TTR.D=NA
sisters$age=NA
for(i in 1:nrow(sisters)){
  sisters[i,"range_overlap"]=jacc.scores[sisters[i,"sp1"],sisters[i,"sp2"]]
  sisters[i,"TTR.D"]=sch.mat[sisters[i,"sp1"],sisters[i,"sp2"]]
  sisters[i,"age"]=dispRity::tree.age(phy)[getMRCA(phy,as.character(sisters[i,1:2])),"ages"]
}

plot(sisters$range_overlap,sisters$TTR.D)

## filter by habitat
sisters_M=sisters[sisters$sp1_category=="M"&sisters$sp2_category=="M",]
sisters_ML=sisters[sisters$sp1_category %in% c("M","ML")&sisters$sp2_category%in% c("M","ML"),]
sisters_L=sisters[sisters$sp1_category %in% c("L")|sisters$sp2_category%in% c("L"),]

plot(sisters_ML$range_overlap,sisters_ML$TTR.D)
plot(sisters_L$range_overlap,sisters_L$TTR.D)

plot(sisters$age,sisters$range_overlap)
plot(sisters$age,sisters$TTR.D)

## regressions
sis.fit=lm(sisters$TTR.D ~ sisters$range_overlap)
summary(sis.fit)
abline(sis.fit)
sis.fit.age=lm(sisters$TTR.D ~ sisters$range_overlap + sisters$age)
summary(sis.fit.age)

sisM.fit=lm(sisters_M$TTR.D ~ sisters_M$range_overlap)
summary(sisM.fit)
abline(sisM.fit)
sisM.fit.age=lm(sisters_M$TTR.D ~ sisters_M$range_overlap + sisters_M$age)
summary(sisM.fit.age)

sisML.fit=lm(sisters_ML$TTR.D ~ sisters_ML$range_overlap)
summary(sisML.fit)
abline(sisML.fit)

## compare niche distance of sisters to point distance
veronica.distances=readRDS(file.path(dist_dir,"Veronica_geo_distances.Rdata"))
sisters$pt_dist=NA
for(i in 1:nrow(sisters)){
  sisters[i,"pt_dist"]=veronica.distances[sisters[i,"sp1"],sisters[i,"sp2"]]
  
}
sis.dist.fit=lm(sisters$TTR.D ~ sisters$pt_dist)
summary(sis.dist.fit)


