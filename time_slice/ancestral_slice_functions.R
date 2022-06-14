

prep_traits_tree <- function(trait_file,tree_file){
  require(ape)
  require(dispRity)
  require(Claddis)
  full_traits <- c("tmax1","tmax2","tmax3","tmax4","q1","q2","w11","w12","ns1","ns2","tmean1","tmean2",
                "w21","w22","w23","w24","nsoil1","nsoil2","tmin1","tmin2","tmin3","tmin4","tmean21","tmean22")
  
  phy<-read.tree(tree_file)
  ## read traits from nexus file into claddis object as required by downstream functions
  ttr_traits_claddis<-Claddis::read_nexus_matrix(file=trait_file) #equalize_weights = TRUE causes error
  ## prune tree and data to match tips
  ttr_traits_phy<-dispRity::clean.data(ttr_traits_claddis$matrix_1$matrix,phy)
  ## replace data matrix with pruned data
  ttr_traits_claddis$matrix_1$matrix <- ttr_traits_phy$data
  ## convert continuous traits to numeric
  ttr_traits_claddis$matrix_1$matrix<-apply(ttr_traits_claddis$matrix_1$matrix, 2, as.numeric)
  ## add row names back in, column names from traits
  row.names(ttr_traits_claddis$matrix_1$matrix)<-row.names(ttr_traits_phy$data)
  colnames(ttr_traits_claddis$matrix_1$matrix)<-full_traits
  ## make character weights equal
  ttr_traits_claddis$matrix_1$character_weights <- rep(1,length(ttr_traits_claddis$matrix_1$character_weights))
  
  ## remove rogue values that were read in as E+ instead of E-, change to 0
  rogue <- which(ttr_traits_claddis$matrix_1$matrix>100)
  if(length(rogue)>0) print(paste("replacing rogue large values with 0 at indices",rogue))
  ttr_traits_claddis$matrix_1$matrix[rogue] <- 0
  
  ## add info for later functions
  ttr_traits_phy$tree$root.time <- max(tree.age(ttr_traits_phy$tree)[, 1])
  ttr_traits_phy$tree <-makeNodeLabel(ttr_traits_phy$tree, method = "number", prefix = "n")
  
  return(list(traits_matrix=ttr_traits_claddis,traits_phy=ttr_traits_phy))
}


subset_traits <- function(trait_input,eigen_file,num_to_keep="all"){
  full_traits <- c("tmax1","tmax2","tmax3","tmax4","q1","q2","w11","w12","ns1","ns2","tmean1","tmean2",
                   "w21","w22","w23","w24","nsoil1","nsoil2","tmin1","tmin2","tmin3","tmin4","tmean21","tmean22")
  ## read in annotated eigenvector loadings csv
  eigen_ranks <- read.csv(eigen_file)
  ## filter and rank by annotations and loading score
  eigen_ranks.f<-eigen_ranks[eigen_ranks$include,]
  eigen_ranks.f<-eigen_ranks.f[-which(is.na(eigen_ranks.f$variable)),]
  eigen_ranks.f<-eigen_ranks.f[order(eigen_ranks.f$variable.ranking,decreasing = TRUE),]
  if(num_to_keep=="all"){
    traits_to_keep<-eigen_ranks.f$variable
  } else if(is.numeric(num_to_keep)){
    traits_to_keep<-eigen_ranks.f$variable[1:num_to_keep]
  } else { stop("num_to_keep should be \"all\" or a number")}
  traits_to_keep_idx<-which(full_traits %in% traits_to_keep)
  
  ## remove rejected trait columns
  traits_sub<-trait_input$traits_matrix
  traits_sub$matrix_1$matrix <- traits_sub$matrix_1$matrix[,traits_to_keep]
  traits_sub$matrix_1$ordering <- traits_sub$matrix_1$ordering[traits_to_keep_idx]
  traits_sub$matrix_1$character_weights <- traits_sub$matrix_1$character_weights[traits_to_keep_idx]
  traits_sub$matrix_1$minimum_values <- traits_sub$matrix_1$minimum_values[traits_to_keep_idx]
  traits_sub$matrix_1$maximum_values <- traits_sub$matrix_1$maximum_values[traits_to_keep_idx]
  trait_input$traits_matrix <- traits_sub
  return(trait_input)
}


get_morphospace <- function(traits_matrix,best_model,plot=FALSE){
  mat<-best_model$mvanc$estimates
  if(ncol(mat)!=ncol(traits_matrix$matrix_1$matrix )) stop("traits matrix and ancestral states' col numbers differ")
  
  ## add nodes to tip matrix
  traits_matrix$matrix_1$matrix <- rbind(traits_matrix$matrix_1$matrix, mat) 

  ## distance matrix
  matrix_dist <- calculate_morphological_distances(traits_matrix)

  ## Ordinating the matrix
  cmdout<-cmdscale(matrix_dist$distance_matrix, k = nrow(matrix_dist$distance_matrix) - 2, add = TRUE)
  morphospace <- cmdout$points
  col<-c("red","blue")
  if(plot) plot(morphospace,col=col[grepl("V",row.names(morphospace))*1+1],main="morphospace: red=nodes,blue=tips") #red=nodes, blue=tips
  
  return(morphospace)
}


get_morphospace_tips <- function(traits_matrix,cols){
    ## distance matrix
  matrix_dist <- calculate_morphological_distances(traits_matrix)
  
  ## Ordinating the matrix
  cmdout<-cmdscale(matrix_dist$distance_matrix, k = nrow(matrix_dist$distance_matrix) - 2, add = TRUE)
  morphospace <- cmdout$points
  col<-c("red","blue")
  plot(morphospace,col=cols,pch=19) #red=nodes, blue=tips
  
  return(morphospace)
}

get_morphospace_nmds <- function(traits_matrix,best_model){
  mat<-best_model$mvanc$estimates
  if(ncol(mat)!=ncol(traits_matrix$matrix_1$matrix )) stop("traits matrix and ancestral states' col numbers differ")
  
  ## add nodes to tip matrix
  traits_matrix$matrix_1$matrix <- rbind(traits_matrix$matrix_1$matrix, mat) 
  
  ## distance matrix
  matrix_dist <- dist(traits_matrix$matrix_1$matrix)
  
  ## Ordinating the matrix
  nmdsout<-monoMDS(matrix_dist, k = 3,maxit=1000)
 
  return(nmdsout)
}

check_dispersal=function(anc,desc,target_areas=c('A','S')){
  
  if(length(anc)==length(desc)){
    out<-c()
    for(i in 1:length(anc)){
      anci<-unlist(strsplit(anc[i],split=""))
      desci<-unlist(strsplit(desc[i],split=""))
      target_desc<-desci[desci %in% target_areas]
      out<-c(out,(any(target_areas %in% desci) & any(!target_desc %in% anci)))
    }
  } else stop("length of ancestors does not equal descendants")
  
  return(out)
}

check_dispersal_value<-function(anc,desc){
  
  if(length(anc)==length(desc)){
    out<-c()
    for(i in 1:length(anc)){
      anci<-unlist(strsplit(anc[i],split=""))
      desci<-unlist(strsplit(desc[i],split=""))
      out<-c(out,paste0(desci[which(!desci %in% anci)],collapse=""))
    }
  } else stop("length of ancestors does not equal descendants")
  
  return(out)
}


prep_node_info <- function(clado_table,areas,area_codes,focal_areas,tree,morphospace){
  require(dplyr)
  ## simplify clado table and code changes from ancestral node's region
  node_info<-clado_table[,c("node","node.type","ancestor","time_bp","label","clado_event_txt","sampled_states_AT_nodes","sampled_states_AT_brbots")]
  node_info<-dplyr::rename(node_info,current_area=sampled_states_AT_nodes,ancestor_area=sampled_states_AT_brbots)
  ## remove duplicated nodes from different time strata
  dup<-node_info[duplicated(node_info$node),"node"]
  node_info<-node_info[-which(node_info$clado_event_txt=='' & node_info$node %in% dup),]
  
  ## identify cases of dispersal to a focal area
  node_info$ancestor_area_txt<-area_codes[node_info$ancestor_area]
  node_info$current_area_txt<-area_codes[node_info$current_area]
  node_info$dispersed_to_focal<-check_dispersal(node_info$ancestor_area_txt,node_info$current_area_txt,focal_areas)
  node_info$dispersed_to_txt<-check_dispersal_value(node_info$ancestor_area_txt,node_info$current_area_txt)
  
  ## match matrix and tree labels/matrix order
  num.inNodes<-length(tree$node.label)
  num.tips<-length(tree$tip.label)
  node_info2<-data.frame(mat_label=c(tree$tip.label,tree$node.label),tree_label=1:(num.tips+num.inNodes))
  node_info2<-dplyr::left_join(node_info2,data.frame(mat_label=row.names(morphospace),mat_num=1:nrow(morphospace)))

  ## combine info
  node_info<-dplyr::left_join(node_info,node_info2,by=c("node"="tree_label"))
  node_info<-node_info[order(node_info$time_bp,decreasing=TRUE),]
  return(node_info)
}


slice_distance <- function(morphospace,tree,node_info,areas,area_codes,focal_areas,plot_title_suffix,tips=FALSE,distdf_only=FALSE){
  require(ggplot2)
  ## slice tree at each node and select ancestral nodes of intersected branches
  ## deltran is supposed to select descendants but it selects ancestors
  slice.anc<-chrono.subsets(morphospace,tree,method="continuous",time=c(node_info[node_info$node.type=="internal","time_bp"],0),inc.nodes=TRUE,model="deltran")
  
  subsets<-slice.anc$subsets
  names(subsets)<-c(node_info[node_info$node.type=="internal","mat_label"],"tips")
  
  ## filter subsets by region and calculate focal node area stats
  for(subset in names(subsets)){
    #print(subset)
    ## check there are enough nodes in slice to calcluate disparity
    if(length(subsets[[subset]]$elements)<2) next
    if(subset=="tips") next
    # name of subset is the focal.node (sliced at its branching time) for which to calcluate distance from the slice
    focal.node<-subset
    # extract nodes in this slice
    nodes<-row.names(morphospace)[subsets[[subset]]$elements]
    # identify area of focal node (number code) and save as text (could change to go directly to text from current_area_txt)
    focal.node.area<-node_info[node_info$mat_label==focal.node,"current_area"]
    subsets[[subset]]$focal.area<-area_codes[focal.node.area]
    
    for(area in areas){
      # subset nodes from this slice found in area
      nodes.in.area<-dplyr::filter(node_info,mat_label %in% nodes & current_area %in% grep(area,area_codes)) %>% dplyr::select(mat_label) %>% unlist()
      subsets[[subset]][[area]]$nodes <- nodes.in.area
      # calculate disparity and distance of the area subset only if the focal node is in this area and there are enough other nodes in the area
      if(focal.node.area %in% grep(area,area_codes) & length(nodes.in.area) > 1){
        # extract morphospace row of the focal node
        focal.node.mat<-matrix(morphospace[focal.node,],nrow=1,ncol=length(morphospace[focal.node,]))
        # calculate distance of each node in area from centroid
        nodes.in.area.disp<-centroids(morphospace[nodes.in.area,])
        # calculate mean distance of focal node from each of the other points
        subsets[[subset]][[area]]$focal.node.point.dist<-mean(point.dist(morphospace[nodes.in.area,],focal.node.mat))
        # calculate distance of focal node from the centroid of the other points
        subsets[[subset]][[area]]$focal.node.cen.dist<-point.dist(focal.node.mat,morphospace[nodes.in.area,])
        # calculate overall disparity as mean of distance of points from centroid
        subsets[[subset]][[area]]$focal.area.disp<-mean(nodes.in.area.disp)
        # calcluate ratio between distance of focal node from centroid and the overall disparity (just for reference)
        subsets[[subset]][[area]]$focal.node.ratio<-subsets[[subset]][[area]]$focal.node.cen.dist/mean(nodes.in.area.disp)
      }
    }
  }
  
  ## plot disparity by area and event type in single dataframe
  dist_df<-data.frame(node=character(0),area=character(0),event=character(0),dist=numeric(0),disp=numeric(0),time_bp=numeric(0),num_species=numeric(0))
  
  for(node in names(subsets)){
    if(length(subsets[[node]]$elements)<2) next
    if(node=="tips") next
    for (area in focal_areas){
      ## include focal node if area of the focal node contains overall focal area (e.g. ML contains M)
      ## note: there will be multiple entries for a node with multiple focal areas included (e.g. AS)
      if(area %in% unlist(strsplit(subsets[[node]]$focal.area,split=""))){
        if(length(subsets[[node]][[area]]$nodes)<2) next
        dist<-subsets[[node]][[area]]$focal.node.cen.dist
        disp<-subsets[[node]][[area]]$focal.area.disp
        # if focal node was from dispersal to focal area, count as dispersal, otherwise clado (note: if S->AS, there will be a clado entry for S and disp entry for A)
        event<-ifelse(node_info[node_info$mat_label==node,"dispersed_to_focal"] & 
                       grepl(area,node_info[node_info$mat_label==node,"dispersed_to_txt"]),"dispersal","clado")
        num_species<-length(subsets[[node]][[area]]$nodes)
        new_row<-data.frame(node=node,area=area,event=event,dist=dist,disp=disp,time_bp=node_info[node_info$mat_label==node,"time_bp"],num_species=num_species)
        dist_df<-rbind(dist_df,new_row)
      }
    }
  }
  
  ## if desired, skip remaining calculations and only return dist_df with individual node measurements
  if(distdf_only) return(list(dist_df=dist_df))
  
  ## tips
  if(tips){
    for(tip in node_info[node_info$node.type=="tip","mat_label"]){
      print(tip)
      focal.node<-tip
      nodes<-row.names(morphospace)[subsets[["tips"]]$elements]
      focal.node.area<-node_info[node_info$mat_label==focal.node,"current_area"]
      subsets[["tips"]][[tip]]$focal.area<-area_codes[focal.node.area]
      for(area in areas){
        nodes.in.area<-dplyr::filter(node_info,mat_label %in% nodes & current_area %in% grep(area,area_codes)) %>% dplyr::select(mat_label) %>% unlist()
        subsets[["tips"]][[area]]$nodes <- nodes.in.area
        if(focal.node.area %in% grep(area,area_codes) & length(nodes.in.area) > 1){
          focal.node.mat<-matrix(morphospace[focal.node,],nrow=1,ncol=length(morphospace[focal.node,]))
          nodes.in.area.disp<-centroids(morphospace[nodes.in.area,])
          subsets[["tips"]][[tip]][[area]]$focal.node.point.dist<-mean(point.dist(morphospace[nodes.in.area,],focal.node.mat))
          subsets[["tips"]][[tip]][[area]]$focal.node.cen.dist<-point.dist(focal.node.mat,morphospace[nodes.in.area,])
          subsets[["tips"]][[tip]][[area]]$focal.area.disp<-mean(nodes.in.area.disp)
        }
      }
    }
    ## probably will create separate df
    for(tip in names(subsets$tips)){
      if(!tip %in% node_info[node_info$node.type=="tip","mat_label"]) next
      
      for (area in focal_areas){
        ## check if area of the focal node contains overall focal area (e.g. ML contains M)
        if(area %in% unlist(strsplit(subsets$tips[[tip]]$focal.area,split=""))){
          dist<-subsets$tips[[tip]][[area]]$focal.node.cen.dist
          disp<-subsets$tips[[tip]][[area]]$focal.area.disp
          event<-ifelse(node_info[node_info$mat_label==tip,"dispersed_to_focal"],"dispersal","clado")
          num_species<-length(subsets$tips[[area]]$nodes)
          new_row<-data.frame(node=tip,area=area,event=event,dist=dist,disp=disp,time_bp=0,num_species=num_species)
          dist_df<-rbind(dist_df,new_row)
        }
      }
    }
  }
 
  ## examine relationships of node distance/disparity with time for cladogenesis vs dispersal with plots, linear model, and correlation coefficients
  
  time.scatter<-ggplot(data=dist_df) + geom_point(aes(x=time_bp,y=dist,col=event,pch=area)) + 
    scale_x_reverse() + theme_classic() + ggtitle(paste("Node distance from mountain niche,",plot_title_suffix))
 # print(time.scatter)
  disp.points<-ggplot(data=dist_df) + geom_point(aes(x=time_bp,y=disp,col=event,pch=area)) +
    scale_x_reverse() + theme_classic()
  event.box<-ggplot(data=dist_df) + geom_boxplot(aes(x=area,y=dist,col=event)) + theme_classic()
  
  lm.clado=lm(dist ~ time_bp,data=dist_df[which(dist_df$event=="clado"),])
  #plot(lm.clado)
  summary(lm.clado)
  
  lm.disp=lm(dist ~ time_bp,data=dist_df[which(dist_df$event=="dispersal"),])
  #plot(lm.disp)
  summary(lm.disp)
  
  cor.clado=cor(dist_df[which(dist_df$event=="clado"),"time_bp"],dist_df[which(dist_df$event=="clado"),"dist"])
  cor.dispersal=cor(dist_df[which(dist_df$event=="dispersal"),"time_bp"],dist_df[which(dist_df$event=="dispersal"),"dist"])  
  
  cor.spear.clado=cor(dist_df[which(dist_df$event=="clado"),"time_bp"],dist_df[which(dist_df$event=="clado"),"dist"],method="spearman")
  cor.spear.dispersal=cor(dist_df[which(dist_df$event=="dispersal"),"time_bp"],dist_df[which(dist_df$event=="dispersal"),"dist"],method="spearman")  
  
  ## calculate median node distance for cladogenesis and dispersal (no time element)
  med.dist.clado=median(dist_df[dist_df$event=="clado" & dist_df$area %in% focal_areas,"dist"])
  med.dist.dispersal=median(dist_df[dist_df$event=="dispersal" & dist_df$area %in% focal_areas,"dist"])
  
  ## output
  slice.plot.list<-list()
  slice.plot.list$node_info<-node_info
  slice.plot.list$subsets<-subsets
  slice.plot.list$dist_df<-dist_df
  # to add: tips
  
  slice.plot.list$ggplots<-list()
  slice.plot.list$ggplots$time.scatter<-time.scatter
  slice.plot.list$ggplots$disp.points<-disp.points
  slice.plot.list$ggplots$event.box<-event.box
  
  slice.plot.list$lm.clado=lm.clado
  slice.plot.list$lm.dispersal=lm.disp
  
  slice.plot.list$med.dist.clado=med.dist.clado
  slice.plot.list$med.dist.dispersal=med.dist.dispersal
  
  slice.plot.list$cor_vals=list(cor.clado=cor.clado,cor.dispersal=cor.dispersal,
                                cor.spear.clado=cor.spear.clado,cor.spear.dispersal=cor.spear.dispersal)


  return(slice.plot.list)

}

shuffle_slice <- function(morphospace,tree,node_info,areas,area_codes,focal_areas,plot_title_suffix,num.perm=100,seed=123,tips=FALSE,distdf_only=FALSE){
  shuffled.list<-list(out=list(),nulldist=list())
  set.seed(seed)
  for(i in 1:num.perm){
    print(i)
    shuffled.list$out[[i]]<-list()
    
    ## shuffle
    num.tips<-length(tree$tip.label)
    morphospace.nodes<-morphospace[(num.tips+1):nrow(morphospace),]
    morphospace.shuffle<-rbind(morphospace[1:num.tips,],morphospace.nodes[sample(nrow(morphospace.nodes)),])
    row.names(morphospace.shuffle)<-row.names(morphospace)
    shuffled.list$out[[i]]$morphospace<-morphospace.shuffle
    
    ## slice
    slice_output <- slice_distance(morphospace=morphospace.shuffle,tree=ttr_tree,node_info=node_info,
                                   areas=areas, area_codes=area_codes,focal_areas=focal_areas,
                                   plot_title_suffix=paste(plot_title_suffix,"null",i),tips=tips,distdf_only=distdf_only)
    shuffled.list$out[[i]]$slice_output <- slice_output
  }
  
  ## get median distance for clado and dispersal for each null 
  med.dist.nulldist.clado <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$med.dist.clado))
  med.dist.nulldist.dispersal <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$med.dist.dispersal))
  
  ## get lm slope against time for clado and dispersal for each null 
  slope.nulldist.clado <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$lm.clado$coefficients[2]))
  slope.nulldist.dispersal <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$lm.dispersal$coefficients[2]))
  #slope.dist.clado.p <- unlist(lapply(shuffled.list,function(x) summary(x$slice_output$lm.clado)$coefficients[2,4]))
  
  ## get correlation coefficients against time for clado and dispersal for each null 
  cor.nulldist.clado <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$cor_vals$cor.clado))
  cor.nulldist.dispersal <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$cor_vals$cor.dispersal))
  cor.spear.nulldist.clado <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$cor_vals$cor.spear.clado))
  cor.spear.nulldist.dispersal <- unlist(lapply(shuffled.list$out,function(x) x$slice_output$cor_vals$cor.spear.dispersal))
  
  ## get original nodewise dist_df for each null
  null_distdfs <- lapply(shuffled.list$out,function(x) x$slice_output$dist_df)
  
  ## add key output/summary stats to nulldist list
  shuffled.list$nulldist$dist_dfs <- null_distdfs
  shuffled.list$nulldist$med.dist.nulldist.clado <- med.dist.nulldist.clado
  shuffled.list$nulldist$med.dist.nulldist.dispersal <- med.dist.nulldist.dispersal
  shuffled.list$nulldist$med.dist.nulldist.diff <- med.dist.nulldist.dispersal - med.dist.nulldist.clado
  
  shuffled.list$nulldist$slope.nulldist.clado <- slope.nulldist.clado
  shuffled.list$nulldist$slope.nulldist.dispersal <- slope.nulldist.dispersal
  
  shuffled.list$nulldist$cor.nulldist.clado <- cor.nulldist.clado
  shuffled.list$nulldist$cor.nulldist.dispersal <- cor.nulldist.dispersal
  shuffled.list$nulldist$cor.spear.nulldist.clado <- cor.spear.nulldist.clado
  shuffled.list$nulldist$cor.spear.nulldist.dispersal <- cor.spear.nulldist.dispersal
  
  return(shuffled.list)
}

# shuffle only within dispersal or cladogenesis
shuffle_slice_categories <- function(morphospace,tree,node_info,areas,area_codes,focal_areas,plot_title_suffix,num.perm=100,seed=123,tips=FALSE,distdf_only=FALSE){
  shuffled.list<-list()
  set.seed(seed)
  for(i in 1:num.perm){
    print(i)
    shuffled.list[[i]]<-list()
    
    ## shuffle
    num.tips<-length(tree$tip.label)
    morphospace.nodes<-morphospace[(num.tips+1):nrow(morphospace),]
    
    morphospace.nodes.dispersal <- morphospace.nodes[which(row.names(morphospace.nodes) %in% node_info[node_info$dispersed_to_focal,"mat_label"]),]
    morphospace.nodes.dispersal.sh <- morphospace.nodes.dispersal[sample(nrow(morphospace.nodes.dispersal)),]
    row.names(morphospace.nodes.dispersal.sh)<-row.names(morphospace.nodes.dispersal)
    
    morphospace.nodes.clado <- morphospace.nodes[which(!row.names(morphospace.nodes) %in% node_info[node_info$dispersed_to_focal,"mat_label"]),]
    morphospace.nodes.clado.sh <- morphospace.nodes.clado[sample(nrow(morphospace.nodes.clado)),]
    row.names(morphospace.nodes.clado.sh)<-row.names(morphospace.nodes.clado)
    
    morphospace.shuffle<-rbind(morphospace[1:num.tips,],morphospace.nodes.dispersal.sh,morphospace.nodes.clado.sh)
    shuffled.list[[i]]$morphospace<-morphospace.shuffle
    
    ## slice
    slice_output <- slice_distance(morphospace=morphospace.shuffle,tree=ttr_tree,node_info=node_info,
                                   areas=areas, area_codes=area_codes,focal_areas=focal_areas,
                                   plot_title_suffix=paste(plot_title_suffix,"null",i),tips=tips,distdf_only=distdf_only)
    shuffled.list[[i]]$slice_output <- slice_output
    
  }
  slope.dist.clado <- unlist(lapply(shuffled.list,function(x) x$slice_output$lm.clado$coefficients[2]))
  slope.dist.dispersal <- unlist(lapply(shuffled.list,function(x) x$slice_output$lm.dispersal$coefficients[2]))
  #slope.dist.clado.p <- unlist(lapply(shuffled.list,function(x) summary(x$slice_output$lm.clado)$coefficients[2,4]))
  

  
  shuffled.list$slope.dist.clado <- slope.dist.clado
  shuffled.list$slope.dist.dispersal <- slope.dist.dispersal
  
  return(shuffled.list)
}


run_slice <- function(out_dir,clado_table, areas, area_codes,focal_areas,tree,morphospace,
                      title_base,num.perm=100,tips=FALSE,distdf_only=FALSE,plot=FALSE){
  #prep node info
  print("preparing node info")
  node_info <- prep_node_info(clado_table=clado_table, areas=areas, area_codes=area_codes,
                              focal_areas=focal_areas,tree=ttr_tree,morphospace=morphospace)
  
  ## slice and calculate distance for each node and its slice (also plots)
  print("slicing and calculating disparity")
  slice_output <- slice_distance(morphospace=morphospace,tree=ttr_tree,node_info=node_info,
                                 areas=areas, area_codes=area_codes,focal_areas=focal_areas,
                                 plot_title_suffix=title_base,tips=tips,distdf_only=distdf_only)
  
  ## shuffled null model comparison
  print("creating null distribution and saving to pdf")
  shuffled_null <- shuffle_slice(morphospace=morphospace,tree=ttr_tree,node_info=node_info,
                                 areas=areas, area_codes=area_codes,focal_areas=focal_areas,
                                 plot_title_suffix=title_base,num.perm=num.perm,tips=tips,distdf_only = distdf_only)
  
  out.list=list()  
  
  if(distdf_only){
    print("saving only dist_df for observed and null distributions")
    out.list$dist_df <- slice_output$dist_df 
    out.list$nulldist <- shuffled_null$nulldist
  } else{
    ## plots for observed output
    if(plot) plot_single_dist_lm(slice_output)
    ## plots for observed and null distribution to pdf (if pdf=TRUE)
    ## note there are several toggles in the function arguments for which pvals to calculate
    pvals <- plot_dist_obs_null(slice_output,shuffled_null,title_base=title_base,out_dir,plot=plot,pdf=FALSE)
  
    out.list$dist_df <- slice_output$dist_df 
    out.list$cor_vals <- slice_output$cor_vals
    out.list$lm.clado <-slice_output$lm.clado
    out.list$lm.dispersal <-slice_output$lm.dispersal
    out.list$nulldist <- shuffled_null$nulldist
    out.list$pvals <- pvals
    if(plot) out.list$ggplots <- slice_output$ggplots
  }
  return(out.list)
}



plot_single_dist_lm <- function(slice_output,lm.diag=FALSE){
  require(ggpubr)
  #ggarrange(slice_output$ggplots$time.scatter,slice_output$ggplots$disp.points,ncol=2,nrow=1)
  print(ggarrange(plotlist=slice_output$ggplots,ncol=2,nrow=2))
  
  ## linear model diagnostics
  if(lm.diag){
    layout(matrix(c(1:4),2,2))
    plot(slice_output$lm.clado,ask=FALSE)
    mtext("clado",outer=TRUE,line=-2)
    plot(slice_output$lm.dispersal,ask=FALSE)
    mtext("dispersal",outer=TRUE,line=-2)
  }

}


plot_single_scatter <- function(dist_df,pvals,plot_title_suffix){
    time.scatter<-ggplot(data=dist_df) + geom_point(aes(x=time_bp,y=dist,col=event),size=3) + 
      scale_x_reverse() + theme_classic(base_size = 20) + ggtitle(paste("Node distance from mountain niche,",plot_title_suffix)) + ylab("distance from niche centroid")+
      annotate("text",x=2.6,y=max(dist_df$dist),label=paste("p.clado=",pvals$p.clado.diff,"  p.dispersal=",pvals$p.dispersal.diff)) 
  
 
  return(time.scatter)

}



calc_nulldist_pvals <- function(BSMout_i){

  ## plot histograms of difference between slopes of observed and null for clado and dispersal niche distances
  layout(matrix(c(1:2),1,2))
  # multiply by -1 to reflect reversed time axis toward present
  clado.diff = -1*(BSMout_i$lm.clado$coefficients[2] - BSMout_i$nulldist$slope.nulldist.clado)
  p.clado.diff=length(which(clado.diff <=0))/length(clado.diff)

  
  dispersal.diff=-1*(BSMout_i$lm.dispersal$coefficients[2] - BSMout_i$nulldist$slope.nulldist.dispersal)
  p.dispersal.diff=length(which(dispersal.diff <=0))/length(dispersal.diff)

  return(list(p.clado.diff=p.clado.diff,p.dispersal.diff=p.dispersal.diff))
}


plot_dist_obs_null <- function(slice_output,shuffled_output,title_base,out_dir,plot=TRUE,lm.diag=FALSE,null.plot=FALSE,
                               lm.p=FALSE,cor.p=FALSE,med.diff.p=TRUE,cor.type="spearman",pdf=FALSE){
  p.list=list()
  
  if(pdf){
    plot=TRUE
    pdf(file=file.path(out_dir,paste0("Node_centroid_distance_plots_",title_base,".pdf")),height=8,width=12)
  } 
  ## plot the observed dataset plots
  plot_single_dist_lm(slice_output,lm.diag=lm.diag)
  
  if(med.diff.p){ 
    ## plot histograms of difference between slopes of observed and null for clado and dispersal niche distances
    med.diff = slice_output$med.dist.dispersal - slice_output$med.dist.clado
    obsnull.diff = med.diff - shuffled_output$nulldist$med.dist.nulldist.diff
    p.obsnull.diff=length(which(obsnull.diff <=0))/length(obsnull.diff)
    print(paste("Median distance diff p = ", p.obsnull.diff))
    p.list$p.obsnull.med.diff=p.obsnull.diff
   
    if(plot){
    hist(obsnull.diff,main="Obs-null median niche distance difference (dispersal-clado)")
    mtext(paste0("p=",p.obsnull.diff),adj=0.1)
    }
  }
  
  if(lm.p){ 
    ## plot histograms of difference between slopes of observed and null for clado and dispersal niche distances
    # multiply by -1 to reflect reversed time axis toward present
    clado.diff = -1*(slice_output$lm.clado$coefficients[2] - shuffled_output$nulldist$slope.nulldist.clado)
    p.clado.diff=length(which(clado.diff <=0))/length(clado.diff)
    print(paste("clado slope p = ", p.clado.diff))
    
    if(plot){
      layout(matrix(c(1:2),1,2))
      hist(clado.diff,main="Obs-null clado niche distance slope")
      mtext(paste0("p=",p.clado.diff),adj=0.1)
    }
    
    dispersal.diff=-1*(slice_output$lm.dispersal$coefficients[2] - shuffled_output$nulldist$slope.nulldist.dispersal)
    p.dispersal.diff=length(which(dispersal.diff <=0))/length(dispersal.diff)
    print(paste("dispersal slope p = ", p.dispersal.diff))
   
     if(plot){
      hist(dispersal.diff,main="Obs-null dispersal niche distance slope")
      mtext(paste0("p=",p.dispersal.diff),adj=0.1)
    }
  }
  
  if(cor.p){
    ## plot histograms of difference between slopes of observed and null for clado and dispersal niche distances
    
    # multiply by -1 to reflect reversed time axis toward present
    if(cor.type=="spearman"){
      clado.diff = -1*(slice_output$cor_vals$cor.spear.clado - shuffled_output$nulldist$cor.spear.nulldist.clado)
    } else{
      clado.diff = -1*(slice_output$cor_vals$cor.clado - shuffled_output$nulldist$cor.nulldist.clado)
    }
    print(clado.diff)
    p.clado.diff=length(which(clado.diff <=0))/length(clado.diff)
    print(paste("clado cor p = ", p.clado.diff))
    p.list$p.clado.diff=p.clado.diff
    
    if(plot){
      layout(matrix(c(1:2),1,2))
      hist(clado.diff,main="Obs-null clado niche distance r")
      mtext(paste0("p=",p.clado.diff),adj=0.1)
    }

    if(cor.type=="spearman"){
      dispersal.diff=-1*(slice_output$cor_vals$cor.spear.dispersal - shuffled_output$nulldist$cor.spear.nulldist.dispersal)
    } else{
      dispersal.diff=-1*(slice_output$cor_vals$cor.dispersal - shuffled_output$nulldist$cor.nulldist.dispersal)
    }
    p.dispersal.diff=length(which(dispersal.diff <=0))/length(dispersal.diff)
    print(paste("dispersal cor p = ", p.dispersal.diff))
    p.list$p.dispersal.diff=p.dispersal.diff
    
    if(plot){
      hist(dispersal.diff,main="Obs-null dispersal niche distance r")
      mtext(paste0("p=",p.dispersal.diff),adj=0.1)
    }
  }
  
  if(null.plot){
    ## plot null distribution plots
    for(i in 1:100){
      plot_single_dist_lm(shuffled_output[[i]]$slice_output,lm.diag=lm.diag)
    }
  }
  
  if(pdf) dev.off()
  
  return(p.list)
}

get_slopes_bsm <- function(BSM_out,lm=TRUE){
  df <- data.frame(clado.slope=numeric(),dispersal.slope=numeric(),
                   p.clado=numeric(),p.dispersal=numeric())
  for(i in 1:length(BSM_out)){
    if(lm.diag){
      new.row=data.frame(clado.slope=BSM_out[[i]]$lm.clado$coefficients[2],
                         dispersal.slope=BSM_out[[i]]$lm.dispersal$coefficients[2],
                         p.clado=BSM_out[[i]]$pvals.slope$p.clado.diff,
                         p.dispersal=BSM_out[[i]]$pvals.slope$p.dispersal.diff)
    } else{
          new.row=data.frame(clado.slope=BSM_out[[i]]$lm.clado$coefficients[2],
                       dispersal.slope=BSM_out[[i]]$lm.dispersal$coefficients[2],
                       p.clado=BSM_out[[i]]$pvals$p.clado.diff,
                       p.dispersal=BSM_out[[i]]$pvals$p.dispersal.diff)
    }

    df <- rbind(df,new.row)
  }
  return(df)
}


get_cor_bsm <- function(BSM_out,type="spearman"){
  df <- data.frame(clado.cor=numeric(),dispersal.cor=numeric(),
                   p.clado=numeric(),p.dispersal=numeric())
  for(i in 1:length(BSM_out)){
    if(type=="spearman"){
      new.row=data.frame(clado.cor=BSM_out[[i]]$cor_vals$cor.spear.clado,
                         dispersal.cor=BSM_out[[i]]$cor_vals$cor.spear.dispersal,
                         p.clado=BSM_out[[i]]$pvals$p.clado.diff,
                         p.dispersal=BSM_out[[i]]$pvals$p.dispersal.diff)
    } else{
      new.row=data.frame(clado.cor=BSM_out[[i]]$cor_vals$cor.clado,
                         dispersal.cor=BSM_out[[i]]$cor_vals$cor.dispersal,
                         p.clado=BSM_out[[i]]$pvals$p.clado.diff,
                         p.dispersal=BSM_out[[i]]$pvals$p.dispersal.diff)
    }

    df <- rbind(df,new.row)
  }
  return(df)
}

## experimental
calc_area_disparity <- function(morphospace,tree,node_info,areas,area_codes,focal_areas){
    area_disp = list()
    for(area in areas){
      area_disp[[area]]=list()
      # subset nodes from this slice found in area
      tips.in.area<-dplyr::filter(node_info,node.type=="tip" & current_area %in% grep(area,area_codes)) %>% dplyr::select(mat_label) %>% unlist()
      area_disp[[area]]$tips.in.area=tips.in.area
      # calculate distance of each node in area from centroid
      area_disp[[area]]$tips.in.area.disp<-centroids(morphospace[tips.in.area,])
      # calculate overall disparity as mean of distance of points from centroid
      area_disp[[area]]$area.disp<-mean(area_disp[[area]]$tips.in.area.disp)
      
    }
    node_info = dplyr::filter(node_info,node.type=="tip")
    dist_df=data.frame(node=character(),area=character(),event=character(),dist.from.centroid=numeric())
    for(tip in node_info_tips$mat_label){
      for (area in focal_areas){
        tips.disp=area_disp[[area]]$tips.in.area.disp
        ## check if area of the focal node contains overall focal area (e.g. ML contains M)
        if(area %in% unlist(strsplit(node_info[which(node_info$mat_label==tip),"current_area_txt"],split=""))){
          dist<-tips.disp[which(names(tips.disp)==tip)]
          event<-ifelse(node_info[node_info$mat_label==tip,"dispersed_to_focal"],"dispersal","clado")
          new_row<-data.frame(node=tip,area=area,event=event,dist.from.centroid=dist)
          dist_df<-rbind(dist_df,new_row)
        }
      }
    }
    area_disp$dist_df=dist_df
    return(area_disp)
}
  


## deprecated
plot_dist_lm_plusnull_old <- function(slice_output,shuffled_output,title_base,out_dir,lm.diag=FALSE,null.plot=TRUE){
  pdf(file=file.path(out_dir,paste0("Node_centroid_distance_plots_",title_base,"_plus100null.pdf")),height=8,width=12)
  ## plot the observed dataset plots
  plot_single_dist_lm(slice_output,lm.diag=lm.diag)
  
  ## plot histograms of difference between slopes of observed and null for clado and dispersal niche distances
  layout(matrix(c(1:2),1,2))
  # multiply by -1 to reflect reversed time axis toward present
  clado.diff = -1*(slice_output$lm.clado$coefficients[2] - shuffled_output$nulldist$slope.nulldist.clado)
  p.clado.diff=length(which(clado.diff <=0))/length(clado.diff)
  print(paste("clado p = ", p.clado.diff))
  hist(clado.diff,main="Obs-null clado niche distance slope")
  mtext(paste0("p=",p.clado.diff),adj=0.1)
  
  dispersal.diff=-1*(slice_output$lm.dispersal$coefficients[2] - shuffled_output$nulldist$slope.nulldist.dispersal)
  p.dispersal.diff=length(which(dispersal.diff <=0))/length(dispersal.diff)
  print(paste("dispersal p = ", p.dispersal.diff))
  hist(dispersal.diff,main="Obs-null dispersal niche distance slope")
  mtext(paste0("p=",p.dispersal.diff),adj=0.1)
  
  if(null.plot){
    ## plot null distribution plots
    for(i in 1:100){
      plot_single_dist_lm(shuffled_output[[i]]$slice_output,lm.diag=lm.diag)
    }
  }
  
  dev.off()
  return(list(p.clado.diff=p.clado.diff,p.dispersal.diff=p.dispersal.diff))
}
