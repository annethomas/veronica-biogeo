                    
run_bin_plot_BSM_multi<-function(run_dir,prefix,plotonly=FALSE,areas=c('L','S','A'),
                                 all.codes=c('O','A','S','L','AS','AL','SL','ASL'),bin.width=0.5,colors=c('red','blue','green')){
  if(!plotonly){
    print("running binning for all trees")
    output.list<-get_BSM_rates_events_0.5my_multitrees(run_dir=run_dir,prefix=prefix,areas,all.codes,bin.width)
  } else{
    load(file.path(run_dir,paste0(prefix,"_binned_output_multi.Rdata")))
  }
  outdir=file.path(run_dir,"binned_plots")
  dir.create(outdir)
  print(paste("plotting binned events to",outdir))
  plot_binned_rates_events(output.list[[1]],plots=c("cladogenesis.rates","cladogenesis.events","dispersal.all.rates","dispersal.all.events","lineages.counts","lineages.rates"),outdir=outdir,prefix=prefix,colors=colors,bin.width=bin.width)
  
}

## reworked Javi's code
get_BSM_rates_events_0.5my_multitrees<-function(run_dir,prefix,areas=c('L','S','A'),
                                                all.codes=c('O','A','S','L','AS','AL','SL','ASL'),
                                                bin.width=0.5){
  setwd(run_dir)
  treedirs=list.files(run_dir,pattern="*tree[0-9]+_output")
  if(length(treedirs)==0){
    treedirs=run_dir
  }
  ages=c()
  for(dir in treedirs){
    setwd(dir)
    print(dir)
    if(file.exists("RES_clado_events_tables.Rdata")){
      load("RES_clado_events_tables.Rdata")
      print(length(RES_clado_events_tables))
      ages=c(ages,max(RES_clado_events_tables[[1]]["time_bp"]))
    }
    else{
      stop("RES_...Rdata files missing")
    }
    setwd(run_dir)
  }
  
  max.bins=max(floor(ages/bin.width))+1
  
  #for appending binned ana and clado events of each BSM replicate from each tree
  areas.list.multi<-list()
  for(area.code in areas){
    areas.list.multi[[area.code]]<-list()
    areas.list.multi[[area.code]][["dispersal.into.time.bins.cumsum"]]<-list()
    areas.list.multi[[area.code]][["insitu.time.bins.cumsum"]]<-list()
    areas.list.multi[[area.code]][["dispersal.into.rates"]]<-list()
    areas.list.multi[[area.code]][["insitu.rates"]]<-list()
    areas.list.multi[[area.code]][["lineages.counts"]]<-list()
  }
  
  for(dir in treedirs){
    setwd(dir)
    print(dir)
    if(file.exists("RES_clado_events_tables.Rdata")){
      load("RES_clado_events_tables.Rdata")
      load("RES_ana_events_tables.Rdata")
    }
    
    #for storing binned ana and clado events of each BSM replicate
    areas.list<-list()
    for(area.code in areas){
      areas.list[[area.code]]<-list()
      areas.list[[area.code]][["dispersal.into.time.bins.cumsum"]]<-list()
      areas.list[[area.code]][["insitu.time.bins.cumsum"]]<-list()
      areas.list[[area.code]][["dispersal.into.rates"]]<-list()
      areas.list[[area.code]][["insitu.rates"]]<-list()
      areas.list[[area.code]][["lineages.counts"]]<-list()
      
    }
    
    cat('analysing BSMs','\n')
    for (i in 1:length(RES_clado_events_tables)){
      print(i)
      ## read in BioGeoBEARS 
      a1<-RES_clado_events_tables[[i]]
      a2<-RES_ana_events_tables[[i]]
      
      ## assign events to time bin
      a1$time.bin<-floor(a1$time_bp/bin.width)*bin.width
      a2$time.bin<-floor(a2$time_bp/bin.width)*bin.width
      max.age<-max(a1$time.bin)
      
    
      ## get time bins of each event type in each area
      for(area.code in areas){
        ## bin dispersal events
        print(area.code)
        dispersal.into.time.bins<-sort(c(a1[a1$clado_dispersal_to==area.code,'time.bin'],a2[a2$dispersal_to==area.code,'time.bin']))
        dispersal.into.time.bins.counts<-sapply(seq(from=0,to=max.age,by=bin.width),function(x)length(dispersal.into.time.bins[dispersal.into.time.bins%in%x]))
        dispersal.into.time.bins.counts.rev<-rev(dispersal.into.time.bins.counts)
        dispersal.into.time.bins.cumsum<-cumsum(dispersal.into.time.bins.counts.rev)
        
        ## bin in situ cladogenesis events
        insitu.time.bins<-sort(a1[a1$clado_dispersal_to==''&(a1$sampled_states_AT_nodes%in%grep(area.code,all.codes))&a1$clado_event_txt!='','time.bin'])
        insitu.time.bins.counts<-sapply(seq(from=0,to=max.age,by=bin.width),function(x)length(insitu.time.bins[insitu.time.bins%in%x]))
        insitu.time.bins.counts.rev<-rev(insitu.time.bins.counts)
        insitu.time.bins.cumsum<-cumsum(insitu.time.bins.counts.rev)
        
        ## bin local extinction events
        local.extinctions.time.bins<-sort(a2[a2$extirpation_from==area.code,'time.bin'])
        local.extinctions.time.bins.counts<-sapply(seq(from=0,to=max.age,by=bin.width),function(x)length(local.extinctions.time.bins[local.extinctions.time.bins%in%x]))
        local.extinctions.time.bins.counts.rev<-rev(local.extinctions.time.bins.counts)
        local.extinctions.time.bins.cumsum<-cumsum(local.extinctions.time.bins.counts.rev)

        ## calculate number of lineages in area at each time bin
        lineages.counts.single<-insitu.time.bins.counts.rev+dispersal.into.time.bins.counts.rev-local.extinctions.time.bins.counts.rev
        lineages.counts<-insitu.time.bins.cumsum+dispersal.into.time.bins.cumsum-local.extinctions.time.bins.cumsum
        
        ## calculate rolling per-capita rates
        dispersal.into.rates<-unlist(lapply(2:length(lineages.counts),function(x)dispersal.into.time.bins.counts.rev[x]/lineages.counts[x-1]))
        insitu.rates<-unlist(lapply(2:length(lineages.counts),function(x)insitu.time.bins.counts.rev[x]/lineages.counts[x-1]))
        lineages.rates<-unlist(lapply(2:length(lineages.counts),function(x)lineages.counts.single[x]/lineages.counts[x-1]))
        
        ## add events and rates to BSM list for area
        areas.list[[area.code]][["dispersal.into.time.bins.cumsum"]][[i]]<-dispersal.into.time.bins.cumsum
        areas.list[[area.code]][["insitu.time.bins.cumsum"]][[i]]<-insitu.time.bins.cumsum
        areas.list[[area.code]][["lineages.counts"]][[i]]<-lineages.counts
        
        areas.list[[area.code]][["dispersal.into.rates"]][[i]]<-dispersal.into.rates
        areas.list[[area.code]][["insitu.rates"]][[i]]<-insitu.rates
        areas.list[[area.code]][["lineages.rates"]][[i]]<-lineages.rates
        
      }
    }  
    
    diff.bins=max.bins-length(areas.list[[1]][["dispersal.into.time.bins.cumsum"]][[1]])
    ## append this tree's BSM event lists to overall lists after padding with 0s to match max bins
    for(area.code in areas){
      areas.list[[area.code]][["dispersal.into.time.bins.cumsum"]]<-lapply(areas.list[[area.code]][["dispersal.into.time.bins.cumsum"]],function(x) c(rep(NA,diff.bins),x))
      areas.list.multi[[area.code]][["dispersal.into.time.bins.cumsum"]]<-c(areas.list.multi[[area.code]][["dispersal.into.time.bins.cumsum"]],areas.list[[area.code]][["dispersal.into.time.bins.cumsum"]])
      
      areas.list[[area.code]][["insitu.time.bins.cumsum"]]<-lapply(areas.list[[area.code]][["insitu.time.bins.cumsum"]],function(x) c(rep(NA,diff.bins),x))
      areas.list.multi[[area.code]][["insitu.time.bins.cumsum"]]<-c(areas.list.multi[[area.code]][["insitu.time.bins.cumsum"]],areas.list[[area.code]][["insitu.time.bins.cumsum"]])
      
      areas.list[[area.code]][["lineages.counts"]]<-lapply(areas.list[[area.code]][["lineages.counts"]],function(x) c(rep(0,diff.bins),x))
      areas.list.multi[[area.code]][["lineages.counts"]]<-c(areas.list.multi[[area.code]][["lineages.counts"]],areas.list[[area.code]][["lineages.counts"]])
      
      areas.list[[area.code]][["dispersal.into.rates"]]<-lapply(areas.list[[area.code]][["dispersal.into.rates"]],function(x) c(rep(0,diff.bins),x))
      areas.list.multi[[area.code]][["dispersal.into.rates"]]<-c(areas.list.multi[[area.code]][["dispersal.into.rates"]],areas.list[[area.code]][["dispersal.into.rates"]])
      
      areas.list[[area.code]][["insitu.rates"]]<-lapply(areas.list[[area.code]][["insitu.rates"]],function(x) c(rep(0,diff.bins),x))
      areas.list.multi[[area.code]][["insitu.rates"]]<-c(areas.list.multi[[area.code]][["insitu.rates"]],areas.list[[area.code]][["insitu.rates"]])
      
      areas.list[[area.code]][["lineages.rates"]]<-lapply(areas.list[[area.code]][["lineages.rates"]],function(x) c(rep(0,diff.bins),x))
      areas.list.multi[[area.code]][["lineages.rates"]]<-c(areas.list.multi[[area.code]][["lineages.rates"]],areas.list[[area.code]][["lineages.rates"]])
      
    }

    setwd(run_dir)
    
  }
  
  ### get CI intervals for plots
  
  ## remove Inf
  for(area.code in areas){
    areas.list.multi[[area.code]][["dispersal.into.rates"]]<-lapply(areas.list.multi[[area.code]][["dispersal.into.rates"]],function(x) {x[!is.finite(x)]<-NA;return(x)})
    areas.list.multi[[area.code]][["insitu.rates"]]<-lapply(areas.list.multi[[area.code]][["insitu.rates"]],function(x) {x[!is.finite(x)]<-NA;return(x)})
    areas.list.multi[[area.code]][["lineages.rates"]]<-lapply(areas.list.multi[[area.code]][["lineages.rates"]],function(x) {x[!is.finite(x)]<-NA;return(x)})
  }
  
  
  #get the BSM events and rates CIs
  areas.list.CI<-list()
  for(area.code in areas){
    areas.list.CI[[area.code]]<-list()
    areas.list.CI[[area.code]][["dispersal.into.time.bins.cumsum.CI"]]<-list()
    areas.list.CI[[area.code]][["insitu.time.bins.cumsum.CI"]]<-list()
    areas.list.CI[[area.code]][["lineages.counts.CI"]]<-list()
    areas.list.CI[[area.code]][["dispersal.into.rates.CI"]]<-list()
    areas.list.CI[[area.code]][["insitu.rates.CI"]]<-list()
    areas.list.CI[[area.code]][["lineages.rates.CI"]]<-list()
  }
  
  ## events CI
  num.bins.events=length(areas.list.multi[[1]][["dispersal.into.time.bins.cumsum"]][[1]])
  for (i in 1:num.bins.events){
    for(area.code in areas){
      areas.list.CI[[area.code]][["dispersal.into.time.bins.cumsum.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["dispersal.into.time.bins.cumsum"]],"[[",i)
      areas.list.CI[[area.code]][["insitu.time.bins.cumsum.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["insitu.time.bins.cumsum"]],"[[",i)
      areas.list.CI[[area.code]][["lineages.counts.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["lineages.counts"]],"[[",i)
    }
  }
  for(area.code in areas){
    areas.list.CI[[area.code]][["dispersal.into.time.bins.cumsum.CI"]]<-lapply(areas.list.CI[[area.code]][["dispersal.into.time.bins.cumsum.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
    areas.list.CI[[area.code]][["insitu.time.bins.cumsum.CI"]]<-lapply(areas.list.CI[[area.code]][["insitu.time.bins.cumsum.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
    areas.list.CI[[area.code]][["lineages.counts.CI"]]<-lapply(areas.list.CI[[area.code]][["lineages.counts.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
  }
  
  ## rates CI
  num.bins.rates=length(areas.list.multi[[1]][["dispersal.into.rates"]][[1]])
  for (i in 1:num.bins.rates){
    for(area.code in areas){
      areas.list.CI[[area.code]][["dispersal.into.rates.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["dispersal.into.rates"]],"[[",i)
      areas.list.CI[[area.code]][["insitu.rates.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["insitu.rates"]],"[[",i)
      areas.list.CI[[area.code]][["lineages.rates.CI"]][[i]]<-sapply(areas.list.multi[[area.code]][["lineages.rates"]],"[[",i)
    }
  }
  for(area.code in areas){
    areas.list.CI[[area.code]][["dispersal.into.rates.CI"]]<-lapply(areas.list.CI[[area.code]][["dispersal.into.rates.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
    areas.list.CI[[area.code]][["insitu.rates.CI"]]<-lapply(areas.list.CI[[area.code]][["insitu.rates.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
    areas.list.CI[[area.code]][["lineages.rates.CI"]]<-lapply(areas.list.CI[[area.code]][["lineages.rates.CI"]],function(x)quantile(x,c(0.05,0.95,0.50),na.rm=TRUE))
    
  }
  
   # collate all
 
  output.list<-list("areas.list.CI"=areas.list.CI,"areas.list.multi"=areas.list.multi)
  
  save(output.list,file=file.path(run_dir,paste0(prefix,"_binned_output_multi.Rdata")))
  
  return(output.list)
}

## for regions with a max age, replace padded 0s with NA (bit of a hack)
pad.NA=function(binned.var,bin.width,max.year=3.5){
  num.bins = length(binned.var)
  for(i in 1:num.bins){
    year=(num.bins-i)*bin.width
    print(paste(i,year))
    if(year > max.year){
      print("pad")
      binned.var[[i]] = rep(NA,3)
    }
  }
  return(binned.var)
}


##############
## plotting ##
##############
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

lines.poly.binplot=function(binned.CI,area.code,var,color,lty=1,fade=FALSE,x.start=1){
  if(fade){
    alpha.line=0.3
    alpha.poly=0.05
  }else{
    alpha.line=1
    alpha.poly=0.2
  }
  lines(x=c(1:length(binned.CI[[area.code]][[var]])),
        y=unlist(lapply(binned.CI[[area.code]][[var]],function(x)x[3])),
        col=adjustcolor( color, alpha.f = alpha.line),lty=lty,lwd=2)
  polygon(x=c((1:length(binned.CI[[area.code]][[var]])),(rev(1:length(binned.CI[[area.code]][[var]])))),
          y=c((unlist(lapply(binned.CI[[area.code]][[var]],function(x)x[1]))),
              rev(unlist(lapply(binned.CI[[area.code]][[var]],function(x)x[2])))),
          col=adjustcolor( color, alpha.f = alpha.poly), border=NA)
}



get.ylim<-function(binned.CI,vars,areas=NULL,ylim.pad=0){
  max=0
  if(is.null(areas)){
    for(var in vars){
      y.nums=unlist(lapply(binned.CI,'[[',var)) ## include all areas
      max=max(max,y.nums[is.finite(y.nums)])
      
    }
  } else{
      for(area in areas){
        for(var in vars){
         y.nums=unlist(binned.CI[[area]][[var]])
         max=max(max,y.nums[is.finite(y.nums)])
      
      }
    }
  }
  #return(c(0,max(y.nums[is.finite(y.nums)])))
  ylim.low=-ylim.pad*max
  return(c(ylim.low,max))
}

plot_binned_single_var<-function(binned.CI,var,title,yaxis.title,show.x.axis=TRUE,legend=TRUE,legend.override=NULL,
                                 colors,bin.width,areas,xlim.low=0,ylim=NULL,ylim.pad=0,lty=1){
  num.bin=length(binned.CI[[1]][[var]])
  if(is.null(ylim)) ylim=get.ylim(binned.CI,var,areas,ylim.pad)
  
  if(show.x.axis){
      plot(c(1,1),xlim=c(xlim.low,num.bin), ylim=ylim,
       type='n',xaxt='n',xlab='Time (Ma)',main=title,ylab=yaxis.title,mgp=c(2.5,1,0),cex.axis=1.2,cex.lab=1.5)
  } else{
      plot(c(1,1),xlim=c(xlim.low,num.bin), ylim=ylim,
        type='n',xaxt='n',xlab='',main=title,ylab=yaxis.title,mgp=c(2.5,1,0),cex.axis=1.2,cex.lab=1.5)
  }

  for(area.code in areas){
    lines.poly.binplot(binned.CI=binned.CI,area.code=area.code,var=var,color=colors[which(areas==area.code)],lty=lty)
  }
  if(show.x.axis)  axis(1,at=c(1:num.bin),labels=c(seq(from=num.bin*bin.width-bin.width,to=0,by=-bin.width)),cex.axis=1.2)
  
  if(is.null(legend.override)){
    legend.titles=names(binned.CI)
  } else legend.titles = legend.override
  
  if(legend) legend('topleft',legend.titles,col=colors,lty=lty,lwd=2,cex=1.2,bty='n')
}


## for plotting multiple rates or multiple cumulative events, not mixed
plot_binned_multi_var<-function(binned.CI,vars,title,yaxis.title,var.titles=NULL,legend.override=NULL,legend.pos='topleft',
                                colors,bin.width,areas,fade=c(FALSE,FALSE),darken=TRUE,xlim.low=0,ylim.pad=0){
  colors.vars=list(colors)
  names(colors.vars)=vars[1]
  if(darken){
    for(i in 2:length(vars)){
      colors.vars[[vars[i]]]=sapply(colors.vars[[i-1]], darken,2.4)
    }
  } else{
    for(i in 2:length(vars)){
      colors.vars[[vars[i]]]=colors.vars[[i-1]]
    }
  }

  num.bin=length(binned.CI[[1]][[vars[1]]])
  plot(c(1,1),xlim=c(xlim.low,num.bin),
       ylim=get.ylim(binned.CI,vars,areas,ylim.pad),
       type='n',xaxt='n',xlab='Time (Ma)',ylab=yaxis.title,mgp=c(2.5,1,0),main=title,cex.axis=1.2,cex.lab=1.5)
  for(area.code in areas){
    for(i in 1:length(vars)){
      var=vars[i]
      fadei=fade[i]
      print(var)
      print(fade)
      lines.poly.binplot(binned.CI,area.code,var,colors.vars[[var]][which(areas==area.code)],lty=i,fade=fadei)
    }
      }
  axis(1,at=c(1:num.bin),labels=c(seq(from=num.bin*bin.width-bin.width,to=0,by=-bin.width)),cex.axis=1.2)
  
  #legend('topleft',names(binned.CI),col=colors,lty=1,cex=.7,bty='n')
  if(is.null(legend.override)){
    if(is.null(var.titles)){
      legend.titles=paste(names(binned.CI),sapply(vars,rep,length(names(binned.CI))))
    } else{
      legend.titles=paste(names(binned.CI),sapply(var.titles,rep,length(names(binned.CI))))
    }
  } else { legend.titles=legend.override }
  
  if(is.null(legend.pos)){
    legend('topleft',
           legend.titles,
           col=unlist(colors.vars),
           # lty=unlist(lapply(1:length(vars),rep,length(names(binned.CI)))),
           lty=unlist(lapply(1:length(vars),rep,length(areas))),
           lwd=2,
           cex=1,bty='n')
  }else{
     legend(legend.pos[1],legend.pos[2],
             legend.titles,
             col=unlist(colors.vars),
             # lty=unlist(lapply(1:length(vars),rep,length(names(binned.CI)))),
             lty=unlist(lapply(1:length(vars),rep,length(areas))),
             lwd=2,
             cex=1,bty='n')
  }

}



plot_binned_rates_events<-function(binned.CI,plots,outdir,prefix,colors,bin.width){
  areas=names(binned.CI)
  if('cladogenesis.rates' %in% plots){
    #plot BSM in situ rates
    pdf(file=file.path(outdir,paste0(prefix,"_cladogenesis_rates.pdf")),width=6,height=6)
    plot_binned_single_var(binned.CI,var="insitu.rates.CI",title="In situ cladogenesis",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }
  if('cladogenesis.events' %in% plots){
    pdf(file=file.path(outdir,paste0(prefix,"_cladogenesis_events.pdf")),width=6,height=6)
    plot_binned_single_var(binned.CI,var="insitu.time.bins.cumsum.CI",title="Cumulative in situ cladogenesis events",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }

  if('dispersal.all.rates' %in% plots){
    pdf(file=file.path(outdir,paste0(prefix,"_dispersal_rates.pdf")),width=6,height=6)
    plot_binned_single_var(binned.CI,var="dispersal.into.rates.CI",title="Dispersal into zones",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }
  if('dispersal.all.events' %in% plots){
    pdf(file=file.path(outdir,paste0(prefix,"_dispersal_events.pdf")),width=6,height=6)    
    plot_binned_single_var(binned.CI,var="dispersal.into.time.bins.cumsum.CI",title="Cumulative dispersal events",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }
  if('lineages.rates' %in% plots){
    #plot BSM in situ rates
    pdf(file=file.path(outdir,paste0(prefix,"_lineages_rates.pdf")),width=6,height=6)
    plot_binned_single_var(binned.CI,var="lineages.rates.CI",title="Lineage accumulation rates",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }
  if('lineages.counts' %in% plots){
    pdf(file=file.path(outdir,paste0(prefix,"_cumulative_lineages.pdf")),width=6,height=6)
    plot_binned_single_var(binned.CI,var="lineages.counts.CI",title="Cumulative lineages",colors=colors,bin.width=bin.width,areas=areas)
    dev.off()
  }
}



