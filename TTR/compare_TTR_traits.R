require(ape)
require(nlme)
require(dplyr)
require(ggplot2)
require(phytools)

#script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/"
#install.packages(file.path(script_dir,"TTR.sdm_0.4.tar.gz"),repos=NULL,type="source")
library(TTR.sdm)

data_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/data"
script_dir="C:/Users/aet_a/OneDrive/Documents/Cambridge/Research/scripts/biogeography_github"
TTR_dir = file.path(data_dir,"TTR")
bgb_dir=file.path(data_dir,"biogeobears")
trait_dir=file.path(data_dir,"biogeography_paper/comparative_phylo")
dist_dir=file.path(data_dir,"distribution_data")

sdm.stats.corrected=read.csv(file.path(TTR_dir,"Veronica_sdm_stats_corrected.csv"))

## from phylo_PCA_TTRtraits.R, for selecting traits and matching to phylogeny
traitset="TTR9"
trait_input=readRDS(file=file.path(TTR_dir,paste0("mvmorph_input_",traitset,".Rdata")))
traits_bands_back=read.csv(file.path(TTR_dir,"TTR9_traits_bands_back.csv")) # created in prep_trait_tables.R

########################
## boxplots of traits ##
########################
pdf(file.path(TTR_dir,"stats_figures/TTRtip_traits_bands_plots.pdf"))

## transformed values
g=ggplot(data=traits_bands_back,aes(x=factor(category,levels=c("A","AS","S","SL","L","ASL")))) +  
  theme_classic(base_size=14)

g=ggplot(data=traits_bands_back,aes(x=factor(merged_category,levels=c("M","ML","L")))) +  
  theme_classic(base_size=14) + theme(axis.title.x = element_blank())

g + geom_boxplot(aes(y=tmin3))
g + geom_boxplot(aes(y=q1))
g + geom_boxplot(aes(y=tmax3))
g + geom_boxplot(aes(y=w21))
g + geom_boxplot(aes(y=w11))
g + geom_boxplot(aes(y=tmean22))
g + geom_boxplot(aes(y=tmean2))
g + geom_boxplot(aes(y=tmin2))
g + geom_boxplot(aes(y=tmax1))

dev.off()

## back-transformed values
pdf(file.path(TTR_dir,"stats_figures/TTRtip_traits_bands_plots.pdf"))

g=ggplot(data=traits_bands_back,aes(x=factor(category,levels=c("A","AS","S","SL","L","ASL")))) +  
  theme_classic(base_size=14)

g=ggplot(data=traits_bands_back,aes(x=factor(merged_category,levels=c("M","ML","L")))) +  
  theme_classic(base_size=14) + theme(axis.title.x = element_blank())

g + geom_boxplot(aes(y=tmin3.back))
g + geom_boxplot(aes(y=q1.back))
g + geom_boxplot(aes(y=tmax3.back))
g + geom_boxplot(aes(y=w21.back))
g + geom_boxplot(aes(y=w11.back))
g + geom_boxplot(aes(y=tmean22.back))
g + geom_boxplot(aes(y=tmean2.back))
g + geom_boxplot(aes(y=tmin2.back))
g + geom_boxplot(aes(y=tmax1.back))

dev.off()

## summary values
se=function(x) sd(x)/sqrt(length(x))

trait_summ=traits_bands %>% group_by(merged_category) %>% 
  summarize(mean.tmin2=mean(tmin2),se.tmin2=se(tmin2),
            mean.tmean2=mean(tmean2),se.tmean2=se(tmean2))

##################################
## phylogenetic anova of groups ##
##################################
trait_input$traits_phy$tree$tip.label=gsub("Veronica","V.",trait_input$traits_phy$tree$tip.label)
traits_bands_tips = left_join(data.frame(species=trait_input$traits_phy$tree$tip.label),traits_bands_back)
#test=phylANOVA(trait_input$traits_phy$tree,traits_bands_tips$category,traits_bands_tips$tmean2)

phylanova_list = list()
phylanova_list$ML = list()
phylanova_list$ASL = list()
for(t in colnames(trait_input$traits_matrix)){
  phylanova_list$ML[[t]] = phylANOVA(trait_input$traits_phy$tree,traits_bands_tips$merged_category,traits_bands_tips[,t])
  phylanova_list$ASL[[t]] = phylANOVA(trait_input$traits_phy$tree,traits_bands_tips$category,traits_bands_tips[,t])
  
}

## back transformed version
phylanova_list = list()
phylanova_list$ML = list()
phylanova_list$ASL = list()
for(t in colnames(trait_input$traits_matrix)){
  t=paste0(t,".back")
  phylanova_list$ML[[t]] = phylANOVA(trait_input$traits_phy$tree,traits_bands_tips$merged_category,traits_bands_tips[,t])
  phylanova_list$ASL[[t]] = phylANOVA(trait_input$traits_phy$tree,traits_bands_tips$category,traits_bands_tips[,t])
  
}

phylanova_list$ML[["tmin2.back"]]
phylanova_list$ML[["tmean2.back"]]

for(t in colnames(trait_input$traits_matrix)){
  t=paste0(t,".back")
  print(t)
  print(phylanova_list$ML[[t]])
}

############################
## phylogenetic regression #
############################

library(nlme)
bm<-corBrownian(1, trait_input$traits_phy$tree)
pag=corPagel(1,trait_input$traits_phy$tree)

model01<-gls(tmin2.back~merged_category,data=traits_bands_tips,correlation=bm)
summary(model01)

model02<-gls(tmin2.back~merged_category,data=traits_bands_tips,correlation=pag)
summary(model02)

## AICc??
aic1=AIC(model01)
aic2=AIC(model02)
k1=3      ## penalty for 3 params: 0.2553191; 13 params: 4.333333
k2=k1+1   ## 0.4301075, 5.060241
n=98

aicc1 = aic1 + (2*k1^2 + 2*k1)/(n-k1-1) 
aicc2 = aic2+(2*k2^2 + 2*k2)/(n-k2-1)

## run pgls with both bm and pagel's lambda for 9 traits, checking for non-convergence errors
## this is what's reported in the ms text
pgls_list = list()
pgls_list$ML = list(bm=list(),pag=list())

for(t in colnames(trait_input$traits_matrix)){
  t=paste0(t,".back")
  print(t)
  
  possibleError1 <- tryCatch(
    gls(formula(paste0(t,"~merged_category")),data=traits_bands_tips,correlation=bm),
    error=function(e) e
  )
  possibleError2 <- tryCatch(
    gls(formula(paste0(t,"~merged_category")),data=traits_bands_tips,correlation=pag),
    error=function(e) e
  )
  
  if(!inherits(possibleError1, "error")){
    print("bm")
    pgls_list$ML$bm[[t]] = gls(formula(paste0(t,"~merged_category")),data=traits_bands_tips,correlation=bm)
    
  } else{
    print("bm failed")
    pgls_list$ML$bm[[t]]=NA
  }   
  if(!inherits(possibleError2, "error")) {
    print("pag")
    pgls_list$ML$pag[[t]] = gls(formula(paste0(t,"~merged_category")),data=traits_bands_tips,correlation=pag)
    
  } else{
    print("pag failed")
    pgls_list$ML$pag[[t]]=NA
  }
}

for(t in colnames(trait_input$traits_matrix)){
  t=paste0(t,".back")
  print(t)
  print(summary(pgls_list$ML$bm[[t]]))
  print("")
  print(summary(pgls_list$ML$pag[[t]]))
}

plot(pgls_list$ML$bm[["tmin2.back"]])


################## 
## plots for ms ##
##################
## only plot the significantly different traits
tmin.summ=traits_bands_back %>% group_by(merged_category) %>%
  summarize(max.tmin=max(tmin2.back))

tmean.summ=traits_bands_back %>% group_by(merged_category) %>%
  summarize(max.tmean=max(tmean2.back))


g=ggplot(data=traits_bands_back,aes(x=factor(merged_category,levels=c("M","ML","L")))) +  
  theme_classic(base_size=15) + theme(axis.title.x = element_blank()) +
  scale_x_discrete(labels=c("M"="Mountain","ML"="Both","L"="Lowland")) +
  scale_color_manual(values=c("green","blue","gray"))#+ theme(plot.margin = margin(t=1,b=.5,unit="cm"))

a=g + geom_boxplot(aes(y=tmin2.back)) + 
  geom_jitter(aes(y=tmin2.back,color=merged_category),width=0.2,show.legend = FALSE) + 
  ylab("Minimum temperature (°C)") + 
  geom_text(data=tmin.summ,aes(label=c("b","a","b"),x=merged_category,
                               y=max.tmin,vjust=-.8),size=5) +
  ylim(NA,max(traits_bands_back$tmin2.back)+2)

b=g + geom_boxplot(aes(y=tmean2.back)) + 
  geom_jitter(aes(y=tmean2.back,color=merged_category),width=0.2,show.legend = FALSE) + 
  ylab("Mean temperature (°C)") + 
  geom_text(data=tmean.summ,aes(label=c("b","a","b"),x=merged_category, 
                                y=max.tmean,vjust=-.8),size=5) +
  ylim(NA,max(traits_bands_back$tmean2.back)+2)

pdf(file.path(bgb_dir,"TTRtip_traits_bands_plots_ms.pdf"),height=4.5,width=7)

ggpubr::ggarrange(a,b,labels=c("A","B"),label.x=0.02, label.y=0.05,
                  font.label=list(size=14),legend="none") 

dev.off()



## attempt at trait plotting on tree with phytools; 
## see tip_trait_plotting.R for a different approach
sdm.stats.corrected$sp=gsub("Veronica_","V._",sdm.stats.corrected$sp)
rownames(sdm.stats.corrected)=sdm.stats.corrected$sp
tmin2=sdm.stats.corrected$tmin2.back
names(tmin2)=sdm.stats.corrected$sp
phytools::plotTree.barplot(phy,tmin2)

