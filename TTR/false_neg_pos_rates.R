require(dplyr)

#########################
## false neg/pos rates ##
#########################
sdm.stats=read.csv(file.path(TTR_dir,"Veronica_sdm_stats.csv"))
sdm.stats$fpos.rate=sdm.stats$d.fpos/(sdm.stats$d.fpos+sdm.stats$d.tneg)
sdm.stats$fneg.rate=sdm.stats$d.fneg/(sdm.stats$d.fneg+sdm.stats$d.tpos)

sdm.stats.simple=dplyr::select(sdm.stats,"sp")
sdm.stats.simple=data.frame(species=sdm.stats$sp,n.occ=sdm.stats$d.tpos+sdm.stats$d.fneg,
                            n.abs=sdm.stats$d.fpos+sdm.stats$d.tneg,
                            fpos.rate=sdm.stats$d.fpos/(sdm.stats$d.fpos+sdm.stats$d.tneg),
                            fneg.rate=sdm.stats$d.fneg/(sdm.stats$d.fneg+sdm.stats$d.tpos))

se=function(x) sd(x)/sqrt(length(x))

sdm.stats.summ=sdm.stats.simple %>% 
  summarize(mean.fpos=mean(fpos.rate),se.fpos=se(fpos.rate),
            mean.fneg=mean(fneg.rate),se.fneg=se(fneg.rate))
