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
## currently have only done these analyses for merged with OUM
dataset="merged"
out_dir=file.path(mv_dir,"sortadate50_merged_TTRtips")

ML_slice_file=file.path(out_dir,paste0("slice_output/slice_disparity_MLdistdf_output_50BSMs.Rdata"))
M_slice_file=file.path(out_dir,paste0("slice_output/slice_disparity_output_50BSMs.Rdata"))

#####################################
## niche distance of clado vs disp ##
#####################################

## compile BSM slice results (median distance)
BSM_out=readRDS(M_slice_file)
slice_results=data.frame(BSM=character(),med.dist.dispersal=numeric(),med.dist.clado=numeric(),
                         med.diff=numeric(),pct.change=numeric(),pval=numeric())
for(b in names(BSM_out)){
  bsm=BSM_out[[b]]
  med.dist.dispersal=median(bsm$dist_df[bsm$dist_df$event=="dispersal","dist"])
  med.dist.clado=median(bsm$dist_df[bsm$dist_df$event=="clado","dist"])
  med.diff = med.dist.dispersal - med.dist.clado
  pct.change =  (med.dist.dispersal-med.dist.clado)/med.dist.clado
  
  if(is.null(bsm$pvals$p.obsnull.med.diff)){
    # forgot to change the output of the shuffle function from the lm pval default; have to recalculate pval
    obsnull.diff = med.diff - bsm$nulldist$med.dist.nulldist.diff
    p.obsnull.diff=length(which(obsnull.diff <=0))/length(obsnull.diff)
  } else p.obsnull.diff=bsm$pvals$p.obsnull.med.diff

  new.row=data.frame(BSM=b,med.dist.dispersal=med.dist.dispersal,med.dist.clado=med.dist.clado,
                     med.diff=med.diff,pct.change=pct.change,pval=p.obsnull.diff)
  slice_results=rbind(slice_results, new.row)
  
}

## make long version for ggplot
slice_results_long = tidyr::pivot_longer(slice_results,cols=c(med.dist.clado,med.dist.dispersal),
                                         names_to="event_type",names_prefix="med.dist.",values_to="med.dist")


write.csv(slice_results,file.path(out_dir,"slice_output/slice_summary_wide_50.csv"),row.names=FALSE)
write.csv(slice_results_long,file.path(out_dir,"slice_output/slice_summary_long_50BSMs.csv"),row.names=FALSE)

## several options for plotting clado vs dispersal median niche distance across BSMs
slice_dir=file.path(out_dir,"slice_output")
pdf(file=file.path(slice_dir,"nichedist_medians_boxplot.pdf"),width=6,height=4)
ggplot(data=slice_results_long) + geom_boxplot(aes(x=event_type,y=med.dist,fill=event_type)) + theme_classic()
dev.off()

pdf(file=file.path(slice_dir,"nichedist_density.pdf"),width=6,height=4)
ggplot(data=slice_results_long) + geom_density(aes(x=med.dist,fill=event_type))+ 
  xlab("Median niche distance")+ theme_classic(base_size = 16)
dev.off()

pdf(file=file.path(slice_dir,"nichedist_hist_50BSMs.pdf"),width=6,height=5)
ggplot(data=slice_results_long) + geom_histogram(aes(x=med.dist,fill=event_type),binwidth = 0.005) + 
  xlab("Median niche distance") +  ylab("Count") + theme_classic(base_size = 16) + 
  theme(legend.position = c(0.7,0.7)) + scale_fill_discrete(name="Event type",labels=c("cladogenesis","dispersal"))
dev.off()


## mean/CI diff across BSMs (for ms text)
mean(slice_results$pct.change) #0.2114148
quantile(slice_results$pct.change,c(0.05,0.95))
#       5%       95% 
# 0.1136388 0.3048607 


#######################################
## time slice disparity through time ##
#######################################

### plot mountain and lowland disparity for all BSMs
######################################
BSM_out_ML=readRDS(ML_slice_file)

## get min and max for plotting
mindisp=min(unlist(lapply(BSM_out_ML,function(x) min(x$dist_df$disp))))
maxdisp=max(unlist(lapply(BSM_out_ML,function(x) max(x$dist_df$disp))))
maxtime=max(unlist(lapply(BSM_out_ML,function(x) max(x$dist_df$time_bp))))

## current Fig in MS
pdf(file.path(out_dir,"slice_output/timeslice_dtt_50BSMs.pdf"),height=5.5,width=6)
plot(BSM_out_ML[[1]]$dist_df[BSM_out_ML[[1]]$dist_df$area=="L","time_bp"],
     BSM_out_ML[[1]]$dist_df[BSM_out_ML[[1]]$dist_df$area=="L","disp"],col="white",
     ylim=c(mindisp,maxdisp),xlim=rev(c(0,maxtime)),xlab="Time (mya)",ylab="Disparity")
for(i in names(BSM_out_ML)) lines(BSM_out_ML[[i]]$dist_df[BSM_out_ML[[i]]$dist_df$area=="L","time_bp"],
                                  BSM_out_ML[[i]]$dist_df[BSM_out_ML[[i]]$dist_df$area=="L","disp"],
                                  xlim=rev(c(0,maxtime)),
                                  col=adjustcolor("green",alpha.f=0.75))
for(i in names(BSM_out_ML)) lines(BSM_out_ML[[i]]$dist_df[BSM_out_ML[[i]]$dist_df$area=="M","time_bp"],
                                  BSM_out_ML[[i]]$dist_df[BSM_out_ML[[i]]$dist_df$area=="M","disp"],
                                  xlim=rev(c(0,maxtime)),
                                  col=adjustcolor("blue",alpha.f=0.5))
dev.off()
# could try for cleaner version with M v L and mean/CI

### compare each BSM ML disparity to its null distribution
######################################
pdf(file.path(out_dir,"disp_vs_null_ML_BSMs.pdf"))

for(i in 1:length(BSM_out_ML)){
  print(i)
  dist_df=BSM_out_ML[[i]]$dist_df
  null_distdfs=BSM_out_ML[[i]]$nulldist$dist_dfs
  mindisp=min(unlist(lapply(null_distdfs,function(x) min(x$disp))))
  maxdisp=max(unlist(lapply(null_distdfs,function(x) max(x$disp))))
  
  plot(dist_df[dist_df$area=="L","time_bp"],dist_df[dist_df$area=="L","disp"],col="white",ylim=c(mindisp,maxdisp))
  for(i in 1:100) points(null_distdfs[[i]][null_distdfs[[i]]$area=="L","time_bp"],
                         null_distdfs[[i]][null_distdfs[[i]]$area=="L","disp"],col="darkolivegreen1")
  for(i in 1:100) points(null_distdfs[[i]][null_distdfs[[i]]$area=="M","time_bp"],
                         null_distdfs[[i]][null_distdfs[[i]]$area=="M","disp"],col="cadetblue1")
  lines(dist_df[dist_df$area=="M","time_bp"],dist_df[dist_df$area=="M","disp"],col="blue")
  lines(dist_df[dist_df$area=="L","time_bp"],dist_df[dist_df$area=="L","disp"],col="green")
}

dev.off()


## non-ML versions
## compare a single BSM disparity to null distribution (mountains only)
####################################
BSM_out=readRDS(file.path(out_dir,paste0("slice_output/slice_disparity_output_goodBSMs.Rdata")))
null_distdfs=BSM_out$`1`$nulldist$dist_dfs
mindisp=min(unlist(lapply(null_distdfs,function(x) min(x$disp))))
maxdisp=max(unlist(lapply(null_distdfs,function(x) max(x$disp))))
# plot observed BSM values and each of the 100 null sims
plot(BSM_out$`1`$dist_df$time_bp,BSM_out$`1`$dist_df$disp,col="blue",ylim=c(mindisp,maxdisp))
for(i in 1:100) points(null_distdfs[[i]]$time_bp,null_distdfs[[i]]$disp)
points(BSM_out$`1`$dist_df$time_bp,BSM_out$`1`$dist_df$disp,col="yellow")

# plot all BSMs (no null) for mountains
mindisp=min(unlist(lapply(BSM_out,function(x) min(x$dist_df$disp))))
maxdisp=max(unlist(lapply(BSM_out,function(x) max(x$dist_df$disp))))
maxtime=max(unlist(lapply(BSM_out,function(x) max(x$dist_df$time_bp))))
plot(BSM_out$`1`$dist_df$time_bp,BSM_out$`1`$dist_df$disp,col="blue",ylim=c(mindisp,maxdisp),xlim=c(0,maxtime))
for(i in names(BSM_out)) points(BSM_out[[i]]$dist_df$time_bp,BSM_out[[i]]$dist_df$disp)
for(i in names(BSM_out)) lines(BSM_out[[i]]$dist_df$time_bp,BSM_out[[i]]$dist_df$disp)









### look at mvMORPH output
results=data.frame(BSM=character(0),aicc=numeric(0),reliability=numeric(0))
mvout_files=list.files(file.path(out_dir,"mvmorph_output"),pattern="*.Rdata")
for(f in mvout_files){
  mv=readRDS(file.path(out_dir,"mvmorph_output",f))
  aicc=mv$mvfit$AICc
  reliable=mv$mvfit$hess.values
  newrow=data.frame(BSM=f,aicc=aicc,reliability=reliable)
  results=rbind(results,newrow)
}

