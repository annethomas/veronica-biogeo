library(BioGeoBEARS)
library(ape)
library(cladoRcpp)
bgb_dir="~/biogeobears"


## run 6 BioGeoBEARS models (Nick Matzke's script) and generate BSMs for posterior distribution (or other set) of trees
run_BioGeoBEARS_plusBSM_multitrees<-function(tree_dir,geogfn, run_dir, BSM=TRUE,seed=12345){
  names.trees<-list.files(tree_dir)
  prefix="sortadate50"
  #names.trees=paste0(prefix,"_tree",1:length(trees)) ## alternative if reading in trees in one file
  model.list<-list(0)
  model.df<-data.frame(tree=character(0),model=character(0),weight=numeric(0))

  for(i in 1:length(names.trees)){
      name<-gsub("\\.tre","",names.trees[i])
      print(name)
      outdir=file.path(run_dir,paste0(name,"_output"))
      dir.create(outdir)
      if(file.exists(file.path(outdir,paste0(name,"_restable_AICc_rellike_formatted.txt")))){
         print("results files exist")
      }else{
        run_biogeobears(trfn = file.path(tree_dir,names.trees[i]),geogfn=geogfn,outdir=outdir,name=name)
      }
      
      ## store best AICc
      aicc.file<-list.files(outdir,pattern='_AICc_rellike_formatted.txt',full.names=TRUE)
      table<-read.table(aicc.file,header=T,sep='\t')
      table<-table[order(-table$AICc_wt),]
      cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
      model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
      model.df<-rbind(model.df,data.frame(tree=name,model=row.names(table)[1],weight=table[1,'AICc_wt']))

  }
  print(model.df)
  write.csv(model.df,file.path(run_dir,"models.csv"))

   if(BSM){
     results.BSM<-list(0)
     for(i in 1:length(names.trees)){
        name<-gsub("\\.tre","",names.trees[i])
        print(name)
        outdir=file.path(run_dir,paste0(name,"_output"))
        resfn=paste(name,model.list[[i]][1],"M0_unconstrained_v1.Rdata",sep="_")
        print(resfn)
        if(file.exists(file.path(outdir,"RES_clado_events_tables.Rdata"))){
         print("results files exist")
         # add: load exsiting results to results.BSM[[i]]
        }else{

        sink(file.path(outdir,'BSM_output.txt'))
        results.BSM[[i]]<-run_bgb_BSM(resfn=resfn,trfn=file.path(tree_dir,names.trees[i]),geogfn=geogfn,outdir=outdir,name=name,model=model.list[[i]][1],numBSM=100,seed=seed) 
        sink() 
        } 
     
     }
     #setwd(run_dir)
     #assign(paste0(prefix,'alltrees.BSM'),results.BSM)
     #saveRDS(get(paste0(name,'.BSM')),file=paste0(name,'.BSM.RDS'))
  }
}


run_BioGeoBEARS_plusBSM_multitrees_listed<-function(tree_dir,treelist,geogfn, run_dir, BSM=TRUE,seed=12345){
  names.trees<-treelist
  prefix="sortadate50"
  #names.trees=paste0(prefix,"_tree",1:length(trees)) ## alternative if reading in trees in one file
  model.list<-list(0)
  model.df<-data.frame(tree=character(0),model=character(0),weight=numeric(0))

  for(i in 1:length(names.trees)){
      name<-gsub("\\.tre","",names.trees[i])
      print(name)
      outdir=file.path(run_dir,paste0(name,"_output"))
      dir.create(outdir)
      if(file.exists(file.path(outdir,paste0(name,"_restable_AICc_rellike_formatted.txt")))){
         print("results files exist")
      }else{
        run_biogeobears(trfn = file.path(tree_dir,names.trees[i]),geogfn=geogfn,outdir=outdir,name=name)
      }
      
      ## store best AICc
      aicc.file<-list.files(outdir,pattern='_AICc_rellike_formatted.txt',full.names=TRUE)
      table<-read.table(aicc.file,header=T,sep='\t')
      table<-table[order(-table$AICc_wt),]
      cat(i,'   ',table[1,'AICc_wt'],'    ',row.names(table)[1],'\n')
      model.list[[i]]<-c(row.names(table)[1],table[1,'AICc_wt'])
      model.df<-rbind(model.df,data.frame(tree=name,model=row.names(table)[1],weight=table[1,'AICc_wt']))

  }
  print(model.df)
  #write.csv(model.df,file.path(run_dir,"models.csv"))

   if(BSM){
     results.BSM<-list(0)
     for(i in 1:length(names.trees)){
        name<-gsub("\\.tre","",names.trees[i])
        print(name)
        outdir=file.path(run_dir,paste0(name,"_output"))
        resfn=paste(name,model.list[[i]][1],"M0_unconstrained_v1.Rdata",sep="_")
        print(resfn)
        if(file.exists(file.path(outdir,"RES_clado_events_tables.Rdata")) || file.exists(file.path(outdir,"BSM_output.txt"))){
         print("results files exist")
         # add: load exsiting results to results.BSM[[i]]
        }else{

        sink(file.path(outdir,'BSM_output.txt'))
        results.BSM[[i]]<-run_bgb_BSM(resfn=resfn,trfn=file.path(tree_dir,names.trees[i]),geogfn=geogfn,outdir=outdir,name=name,model=model.list[[i]][1],numBSM=100,seed=seed) 
        sink() 
        } 
     
     }
     #setwd(run_dir)
     #assign(paste0(prefix,'alltrees.BSM'),results.BSM)
     #saveRDS(get(paste0(name,'.BSM')),file=paste0(name,'.BSM.RDS'))
  }
}

## run BioGeoBEARS and BSMs in parallel batches
require(parallel)

run_bgb_parallel<-function(batch_size,tree_dir,geogfn, run_dir, BSM,seed=12345){
  
  print(batch_size)
  treelist_full=list.files(tree_dir)
  start <- seq(1,length(treelist_full),by=batch_size)
  print(start)
  run_bgb_parallel <- function(start){
    treelist=treelist_full[start:(start+batch_size-1)] 
    print(treelist)
    run_BioGeoBEARS_plusBSM_multitrees_listed(tree_dir=tree_dir,treelist=treelist,geogfn=geogfn, run_dir=run_dir,BSM=BSM,seed=seed)
    return(0)
  }
  
  out<-mclapply(mc.cores=ceiling(length(treelist_full)/batch_size),X=start,FUN=run_bgb_parallel)
  
}


setup_bgb_runobj = function(trfn,geogfn,run_dir,max_range_size,time_strat=TRUE){
  # Intitialize a default model (DEC model)
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  
  # Give BioGeoBEARS the location of the phylogeny Newick file
  BioGeoBEARS_run_object$trfn = trfn
  
  # Give BioGeoBEARS the location of the geography text file
  BioGeoBEARS_run_object$geogfn = geogfn
  
  # Input the maximum range size
  BioGeoBEARS_run_object$max_range_size = max_range_size
  
  BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
  BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
  
  # Set up a time-stratified analysis:
  # 1. Here, un-comment ONLY the files you want to use.
  # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
  # 3. For example files see (a) extdata_dir, 
  #  or (b) http://phylo.wikidot.com/biogeobears#files
  #  and BioGeoBEARS Google Group posts for further hints)
  #
  # Uncomment files you wish to use in time-stratified analyses:
  if(time_strat){
    BioGeoBEARS_run_object$timesfn = file.path(run_dir,"time_periods.txt")
    BioGeoBEARS_run_object$dispersal_multipliers_fn = file.path(run_dir,"manual_dispersal_multipliers.txt")
    BioGeoBEARS_run_object$areas_allowed_fn = file.path(run_dir,"areas_allowed.txt")
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.
  }
  
  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
  BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
  BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx()
  BioGeoBEARS_run_object$num_cores_to_use = 1
  
  BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
  
  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
  if(time_strat){
    BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table
  }
  
  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  return(BioGeoBEARS_run_object)
}


run_biogeobears<-function(trfn,geogfn,run_dir,name,max_range_size = 6,time_strat){
  tr=read.tree(trfn)
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  
  setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  setwd(run_dir)

  # Run DEC
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  # Set up DEC model
  # (nothing to do; defaults)
  
  # Run this to check inputs. Read the error messages if you get them!
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # For a slow analysis, run once, then set runslow=FALSE to just 
  # load the saved result.
  runslow = TRUE
  resfn = paste0(genus,"_DEC_M0_unconstrained_v1.Rdata")
  if (runslow){
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resDEC = res
  } else {
    # Loads to "res"
    load(resfn)
    resDEC = res
  }
  
  #######################################################
  # Run DEC+J
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDEC$outputs@params_table["d","est"]
  estart = resDEC$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste0(genus,"_DEC+J_M0_unconstrained_v1.Rdata")
  
  
  if (runslow)
  {
    res = bears_optim_run(BioGeoBEARS_run_object)
    res    
    save(res, file=resfn)
    resDECj = res
  } else {
    # Loads to "res"
    load(resfn)
    resDECj = res
  }
  #######################################################
  # PDF plots
  #######################################################
  pdffn = paste0(genus,"_DEC_vs_DEC+J_M0_unconstrained_v1.pdf")
  pdf(pdffn, width=12, height=14)
  
  #######################################################
  # Plot ancestral states - DEC
  #######################################################
  analysis_titletxt = paste0("BioGeoBEARS DEC on ",genus, " M0_unconstrained")
  
  # Setup
  results_object = resDEC
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), 
                                  plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                                  splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                  include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), 
                           plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                           splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  #######################################################
  # Plot ancestral states - DECJ
  #######################################################
  analysis_titletxt = paste0("BioGeoBEARS DEC+J on ",genus, " M0_unconstrained")
  
  # Setup
  results_object = resDECj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()  # Turn off PDF
  
  
  #######################################################
  #######################################################
  # DIVALIKE AND DIVALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
  # Ronquist (1997)'s parsimony DIVA. It is a likelihood
  # interpretation of DIVA, constructed by modelling DIVA's
  # processes the way DEC does, but only allowing the 
  # processes DIVA allows (widespread vicariance: yes; subset
  # sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
  #
  # DIVALIKE is a likelihood interpretation of parsimony
  # DIVA, and it is "like DIVA" -- similar to, but not
  # identical to, parsimony DIVA.
  #
  # I thus now call the model "DIVALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run DIVALIKE
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  # Set up DIVALIKE model
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  
  resfn = paste0(genus, "_DIVALIKE_M0_unconstrained_v1.Rdata")
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
  
  #######################################################
  # Run DIVALIKE+J
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  # Set up DIVALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resDIVALIKE$outputs@params_table["d","est"]
  estart = resDIVALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Remove subset-sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  # Allow classic, widespread vicariance; all events equiprobable
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  # Add jump dispersal/founder-event speciation
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste0(genus,"_DIVALIKE+J_M0_unconstrained_v1.Rdata")
  runslow = TRUE
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
  
  pdffn = paste(genus, "_DIVALIKE_vs_DIVALIKE+J_M0_unconstrained_v1.pdf")
  pdf(pdffn, width=12, height=12)
  
  #######################################################
  # Plot ancestral states - DIVALIKE
  #######################################################
  analysis_titletxt = paste0("BioGeoBEARS DIVALIKE on ", genus, " M0_unconstrained")
  
  # Setup
  results_object = resDIVALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  #######################################################
  # Plot ancestral states - DIVALIKE+J
  #######################################################
  analysis_titletxt =paste0("BioGeoBEARS DIVALIKE+J on ", genus, "M0_unconstrained")
  
  # Setup
  results_object = resDIVALIKEj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()
  #######################################################
  #######################################################
  # BAYAREALIKE AND BAYAREALIKE+J ANALYSIS
  #######################################################
  #######################################################
  # NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is 
  # not identical with the full Bayesian model implemented 
  # in the "BayArea" program of Landis et al. (2013). 
  #
  # Instead, this is a simplified likelihood interpretation
  # of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
  # "d" and "e" work like they do in the DEC model of Lagrange 
  # (and BioGeoBEARS), and then BayArea's cladogenesis assumption
  # (which is that nothing in particular happens at cladogenesis) is 
  # replicated by BioGeoBEARS.
  #
  # This leaves out 3 important things that are in BayArea:
  # 1. Distance dependence (you can add this with a distances 
  #    matrix + the "x" parameter in BioGeoBEARS, however)
  # 2. A correction for disallowing "e" events that drive
  #    a species extinct (a null geographic range)
  # 3. The neat Bayesian sampling of histories, which allows
  #    analyses on large numbers of areas.
  #
  # The main purpose of having a "BAYAREALIKE" model is 
  # to test the importance of the cladogenesis model on 
  # particular datasets. Does it help or hurt the data 
  # likelihood if there is no special cladogenesis process?
  # 
  # BAYAREALIKE is a likelihood interpretation of BayArea,
  # and it is "like BayArea" -- similar to, but not
  # identical to, Bayesian BayArea.
  # I thus now call the model "BAYAREALIKE", and you should also. ;-)
  #######################################################
  #######################################################
  
  #######################################################
  # Run BAYAREALIKE
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  # Set up BAYAREALIKE model
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # No jump dispersal/founder-event speciation
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # Check the inputs
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  runslow = TRUE
  resfn = paste0(genus,"_BAYAREALIKE_M0_unconstrained_v1.Rdata")
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
  
  #######################################################
  # Run BAYAREALIKE+J
  #######################################################
  BioGeoBEARS_run_object=setup_bgb_runobj(trfn,geogfn,run_dir,max_range_size,time_strat)
  
  # Set up BAYAREALIKE+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  dstart = resBAYAREALIKE$outputs@params_table["d","est"]
  estart = resBAYAREALIKE$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # No subset sympatry
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  # No vicariance
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  # Adjust linkage between parameters
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  # Only sympatric/range-copying (y) events allowed, and with 
  # exact copying (both descendants always the same size as the ancestor)
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows 
  # machines. I can't replicate this on my Mac machines, but it is almost certainly
  # just some precision under-run issue, when optim/optimx tries some parameter value 
  # just below zero.  The "min" and "max" options on each parameter are supposed to
  # prevent this, but apparently optim/optimx sometimes go slightly beyond 
  # these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" 
  # slightly for each parameter:
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  resfn = paste0(genus,"_BAYAREALIKE+J_M0_unconstrained_v1.Rdata")
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resBAYAREALIKEj = res
  
  pdffn = paste0(genus,"_BAYAREALIKE_vs_BAYAREALIKE+J_M0_unconstrained_v1.pdf")
  pdf(pdffn, width=12, height=14)
  
  #######################################################
  # Plot ancestral states - BAYAREALIKE
  #######################################################
  analysis_titletxt = paste0("BioGeoBEARS BAYAREALIKE on ",genus, " M0_unconstrained")
  
  # Setup
  results_object = resBAYAREALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  #######################################################
  # Plot ancestral states - BAYAREALIKE+J
  #######################################################
  #png(file="test.png",2000,1200)
  analysis_titletxt = paste0("BioGeoBEARS BAYAREALIKE+J on ",genus, " M0_unconstrained")
  
  # Setup
  results_object = resBAYAREALIKEj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=1, statecex=0.9, splitcex=0.3, titlecex=2, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  # plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
  #                          label.offset=0.45, tipcex=2.3, statecex=1.8, titlecex=2.5, 
  #                          axiscex = 1.8,axispadj = .7,plotsplits=FALSE, cornercoords_loc=scriptdir, 
  #                          include_null_range=TRUE, tr=tr, tipranges=tipranges)
  # #dev.off()
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()
  
  ## adjust visuals
  # pdf("BAYAREALIKE+J_ASL_pie_legend.pdf",height=14,width=12)
  # tr$tip.label=gsub("Veronica_","V. ",tr$tip.label)
  # plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", plotlegend=TRUE, label.offset=0.1, tipcex=1, statecex=0.7, splitcex=0.6, titlecex=1, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  # dev.off()
  
  #########################################################################
  #########################################################################
  #########################################################################
  #########################################################################
  # 
  # CALCULATE SUMMARY STATISTICS TO COMPARE
  # DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  # 
  #########################################################################
  #########################################################################
  #########################################################################
  #########################################################################
  
  #########################################################################
  #########################################################################
  # REQUIRED READING:
  #
  # Practical advice / notes / basic principles on statistical model 
  #    comparison in general, and in BioGeoBEARS:
  # http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
  #########################################################################
  #########################################################################
  
  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL
  
  #######################################################
  # Statistics -- DEC vs. DEC+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DEC, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DEC+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
  # confer the same likelihood on the data. See: Brian O'Meara's webpage:
  # http://www.brianomeara.info/tutorials/aic
  # ...for an intro to LRT, AIC, and AICc
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- DIVALIKE vs. DIVALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # DIVALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #######################################################
  # Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
  #######################################################
  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
  
  numparams1 = 3
  numparams2 = 2
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  stats
  
  # BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  conditional_format_table(stats)
  
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
  
  #########################################################################
  # ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
  #########################################################################
  teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
  teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
  row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
  restable = put_jcol_after_ecol(restable)
  restable
  
  # Look at the results!!
  restable
  teststable
  
  #######################################################
  # Save the results tables for later -- check for e.g.
  # convergence issues
  #######################################################
  
  # Loads to "restable"
  save(restable, file= paste0(name,"_restable_v1.Rdata"))
  load(file=paste0(name,"_restable_v1.Rdata"))
  
  # Loads to "teststable"
  save(teststable, file= paste0(name,"_teststable_v1.Rdata"))
  load(file= paste0(name,"_teststable_v1.Rdata"))
  
  # Also save to text files
  write.table(restable, file= paste0(name,"_restable.txt"), quote=FALSE, sep="\t")
  write.table(unlist_df(teststable), file= paste0(name,"_teststable.txt"), quote=FALSE, sep="\t")
  
  #######################################################
  # Model weights of all six models
  #######################################################
  restable2 = restable
  
  # With AICs:
  AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable = cbind(restable, AICtable)
  restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
  restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
  restable_AIC_rellike
  
  # With AICcs -- factors in sample size
  samplesize = length(tr$tip.label)
  AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
  restable2 = cbind(restable2, AICtable)
  restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
  restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
  restable_AICc_rellike
  
  # Also save to text files
  write.table(restable_AIC_rellike, file= paste0(name,"_restable_AIC_rellike.txt"), quote=FALSE, sep="\t")
  write.table(restable_AICc_rellike, file= paste0(name,"_restable_AICc_rellike.txt"), quote=FALSE, sep="\t")
  
  # Save with nice conditional formatting
  write.table(conditional_format_table(restable_AIC_rellike), file= paste0(name,"_restable_AIC_rellike_formatted.txt"), quote=FALSE, sep="\t")
  write.table(conditional_format_table(restable_AICc_rellike), file= paste0(name,"_restable_AICc_rellike_formatted.txt"), quote=FALSE, sep="\t")
}


##################################################################################

run_bgb_BSM<-function(resfn,trfn,geogfn,run_dir,name,model,numBSM,max_range_size=6,time_strat=TRUE,seed=12345){
  tr=read.tree(trfn)
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = max_range_size
  
  timesfn = file.path(run_dir,"time_periods.txt")
  dispersal_multipliers_fn = file.path(run_dir,"manual_dispersal_multipliers.txt")
  areas_allowed_fn = file.path(run_dir,"areas_allowed.txt")

  setwd(run_dir)
  getwd()
  model_name = paste0(name,model)
  if(file.exists(file.path(run_dir,resfn))){
     load(file.path(run_dir,resfn)) ## loads res
  } else {print(paste(resfn,"does not exist")) }
  
  
  #######################################################
  # Stochastic mapping on specified model stratified with elevational zones
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr) #Anne's note: doesn't seem to match actual structure of stochastic_mapping_inputs_list
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  #set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=numBSM, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=seed, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  #head(clado_events_tables[[1]])
  #head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  

     
  #######################################################
  # Plot all stochastic maps to PDF
  #######################################################
  # Setup  

  include_null_range = TRUE
  areanames = names(tipranges@df)
  areas = areanames
  max_range_size = max_range_size
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified=time_strat
  
  # Loop through the maps and plot to PDF
  pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
  pdf(file=pdffn, width=12, height=14)
  
  nummaps_goal = numBSM
  for (i in 1:nummaps_goal)
  {
    clado_events_table = clado_events_tables[[i]]
    analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
    plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
  } # END for (i in 1:nummaps_goal)
  
  dev.off()
  
  return(BSM_output)
}


## version that includes summarizing (could make separate function)
run_bgb_BSM_full<-function(resfn,trfn,geogfn,run_dir,name,model,numBSM){
  tr=read.tree(trfn)
  tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
  max_range_size = 6
  
  timesfn = file.path(run_dir,"time_periods.txt")
  dispersal_multipliers_fn = file.path(run_dir,"manual_dispersal_multipliers.txt")
  areas_allowed_fn = file.path(run_dir,"areas_allowed.txt")

  setwd(run_dir)
  model_name = paste0(name,model)
  load(file.path(run_dir,resfn)) ## loads res
  
  #######################################################
  # Plot ancestral states - old
  #######################################################
#  pdffn = paste0(genus,"_", model_name, "_v1.pdf")
#  pdf(pdffn, width=6, height=6)
#  
#  analysis_titletxt = paste0(model_name, " on ",genus)
#  
#  # Setup
  results_object = res
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
#  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
#  
#  # Pie chart
#  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
#  
#  dev.off()  # Turn off PDF
#  cmdstr = paste("open ", pdffn, sep="")
#  system(cmdstr) # Plot it
  
  #######################################################
  # Stochastic mapping on DEC M3b stratified with islands coming up
  #######################################################
  clado_events_tables = NULL
  ana_events_tables = NULL
  lnum = 0
  
  #######################################################
  # Get the inputs for Biogeographical Stochastic Mapping
  # Note: this can be slow for large state spaces and trees, since 
  # the independent likelihoods for each branch are being pre-calculated
  # E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
  # for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
  # for storage of "BSM_inputs_file.Rdata".
  # Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
  # the same settings will be used for get_inputs_for_stochastic_mapping().
  #######################################################
  BSM_inputs_fn = "BSM_inputs_file.Rdata"
  runInputsSlow = TRUE
  if (runInputsSlow)
  {
    stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
    save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
  } else {
    # Loads to "stochastic_mapping_inputs_list"
    load(BSM_inputs_fn)
  } # END if (runInputsSlow)
  
  # Check inputs (doesn't work the same on unconstr) #Anne's note: doesn't seem to match actual structure of stochastic_mapping_inputs_list
  names(stochastic_mapping_inputs_list)
  stochastic_mapping_inputs_list$phy2
  stochastic_mapping_inputs_list$COO_weights_columnar
  stochastic_mapping_inputs_list$unconstr
  set.seed(seed=as.numeric(Sys.time()))
  
  runBSMslow = TRUE
  if (runBSMslow == TRUE)
  {
    # Saves to: RES_clado_events_tables.Rdata
    # Saves to: RES_ana_events_tables.Rdata
    BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=numBSM, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
    
    RES_clado_events_tables = BSM_output$RES_clado_events_tables
    RES_ana_events_tables = BSM_output$RES_ana_events_tables
  } else {
    # Load previously saved...
    
    # Loads to: RES_clado_events_tables
    load(file="RES_clado_events_tables.Rdata")
    # Loads to: RES_ana_events_tables
    load(file="RES_ana_events_tables.Rdata")
    BSM_output = NULL
    BSM_output$RES_clado_events_tables = RES_clado_events_tables
    BSM_output$RES_ana_events_tables = RES_ana_events_tables
  } # END if (runBSMslow == TRUE)
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  head(clado_events_tables[[1]])
  head(ana_events_tables[[1]])
  length(clado_events_tables)
  length(ana_events_tables)
  
  #######################################################
  # Plot one stochastic map, manual method
  #######################################################
  # (we have to convert the stochastic maps into event
  #  maps for plotting)
  
  ######################
  # Get the color scheme
  ######################
  include_null_range = TRUE
  areanames = names(tipranges@df)
  areas = areanames
  max_range_size = 4
  
  # Note: If you did something to change the states_list from the default given the number of areas, you would
  # have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  
  ############################################
  # Setup for painting a single stochastic map
  ############################################
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified=TRUE
  clado_events_table = clado_events_tables[[1]]
  ana_events_table = ana_events_tables[[1]]
  
  # cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
  # colnums = match(cols_to_get, names(ana_events_table))
  # ana_events_table_cols_to_add = ana_events_table[,colnums]
  # anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
  # ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
  # rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
  # master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)
  
  ############################################
  # Open a PDF
  ############################################
  pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  # Convert the BSM into a modified res object
  master_table_cladogenetic_events = clado_events_tables[[1]]
  resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)
  
  # Paint on the branch states
  paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)
  
  plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)
  
  ############################################
  # Close PDF
  ############################################
  dev.off()
  cmdstr = paste("open ", pdffn, sep="")
  system(cmdstr)
  
  #######################################################
  # Plot all stochastic maps to PDF
  #######################################################
  # Setup
  include_null_range = include_null_range
  areanames = areanames
  areas = areanames
  max_range_size = max_range_size
  states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
  colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  stratified = stratified
  
  # Loop through the maps and plot to PDF
  pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
  pdf(file=pdffn, width=6, height=6)
  
  nummaps_goal = numBSM
  for (i in 1:nummaps_goal)
  {
    clado_events_table = clado_events_tables[[i]]
    analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
    plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
  } # END for (i in 1:nummaps_goal)
  
  dev.off()
  cmdstr = paste("open ", pdffn, sep="")
  system(cmdstr)
  
  #######################################################
  # Summarize stochastic map tables
  #######################################################
  length(clado_events_tables)
  length(ana_events_tables)
  
  head(clado_events_tables[[1]][,-20])
  tail(clado_events_tables[[1]][,-20])
  
  head(ana_events_tables[[1]])
  tail(ana_events_tables[[1]])
  
  areanames = names(tipranges@df)
  actual_names = areanames
  actual_names
  
  # Get the dmat and times (if any)
  dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
  dmat_times
  
  # Extract BSM output
  clado_events_tables = BSM_output$RES_clado_events_tables
  ana_events_tables = BSM_output$RES_ana_events_tables
  
  # Simulate the source areas
  BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
  clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
  ana_events_tables = BSMs_w_sourceAreas$ana_events_tables
  
  # Count all anagenetic and cladogenetic events
  counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)
  
  summary_counts_BSMs = counts_list$summary_counts_BSMs
  print(conditional_format_table(summary_counts_BSMs))
  
  # Histogram of event counts
  hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))
  
  #######################################################
  # Print counts to files
  #######################################################
  tmpnames = names(counts_list)
  cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
  for (i in 1:length(tmpnames))
  {
    cmdtxt = paste0("item = counts_list$", tmpnames[i])
    eval(parse(text=cmdtxt))
    
    # Skip cubes
    if (length(dim(item)) != 2)
    {
      next()
    }
    
    outfn = paste0(tmpnames[i], ".txt")
    if (length(item) == 0)
    {
      cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
      cat("\n")
    } else {
      cat(outfn)
      cat("\n")
      write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    } # END if (length(item) == 0)
  } # END for (i in 1:length(tmpnames))
  cat("...done.\n")
  
  #######################################################
  # Check that ML ancestral state/range probabilities and
  # the mean of the BSMs approximately line up
  #######################################################
  library(MultinomialCI)    # For 95% CIs on BSM counts
  check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)
}

