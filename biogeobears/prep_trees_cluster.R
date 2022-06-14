#http://coleoguy.blogspot.com/2012/09/randomly-sampling-trees.html
library(ape)

sample_trees<-function(treefile, burnin, final.number, format, prefix, outdir){ 
  trees=read.nexus(treefile) ## might need to add "END;"
  
  #NUMBER OF TREES IN ORIGINAL FILE
  original.number<-as.numeric(length(trees))                           
  
  #THIS CREATES THE POST BURNIN PORTION
  post.burnin.trees<-trees[(burnin*original.number):original.number]   
  
  #THIS DOWNSAMPLES THE COLLECTION OF TREES
  final.trees<-sample(post.burnin.trees, final.number)                 
  
  #THIS SAVES THEM AS NEWICK FORMAT
  if(format=="new"){
	outfile=file.path(outdir,paste0(prefix,"_",final.number,"trees.nwk"))
 	write.tree(final.trees, file=outfile)
	return(outfile)
  }         
  
  #THIS SAVES THEM AS NEXUS FORMAT
  if(format=="nex"){
	outfile=file.path(outdir,paste0(prefix,"_",final.number,"trees.nex"))
	write.nexus(final.trees, file=outfile)	
	return(outfile)
  }
}


prep_trees_BioGeoBEARS<-function(treefile,geogfn,prefix,treeoutdir){
  trees=read.nexus(treefile) ## 100 trees
  geog=read.table(geogfn,skip=1,header = FALSE)
  for(i in 1:length(trees)){
      print(i)
      tree=trees[[i]]
      tree$tip.label=paste0("Veronica_",tree$tip.label)
      print(setdiff(geog$V1,tree$tip.label))
      tips.to.drop=setdiff(tree$tip.label,geog$V1)
      tree_edit=drop.tip(tree,tips.to.drop)
      #trees[[i]]=tree_edit
      write.tree(tree_edit, file.path(treeoutdir,paste0(prefix,"_tree",i,".tre")))
  }

  #write.tree(trees,outfile)
}
