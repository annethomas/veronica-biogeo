## util functions for niche/range overlap prep and processing

clean_matrix=function(x,labels_to_match=NULL){
  rownames(x)[which(rownames(x)=="Veronica_plano-petiolata")] = "Veronica_planopetiolata"
  colnames(x)[which(colnames(x)=="Veronica_plano-petiolata")] = "Veronica_planopetiolata"
  rownames(x)[which(rownames(x)=="Veronica_cockayniana")] = "Veronica_cockayneana"
  colnames(x)[which(colnames(x)=="Veronica_cockayniana")] = "Veronica_cockayneana"
  
  if(!is.null(labels_to_match)){
    rows.to.remove = which(rownames(x) %in% setdiff(rownames(x), labels_to_match))
    print("removing")
    print(rows.to.remove)
    x=x[-rows.to.remove,-rows.to.remove]
  }
  
  return(x)
}

## for preparing TTR niche projection raster
prep_species_raster <- function(r,env.grid,pred.grd,type="pred"){
  grd.df <- data.frame(long=env_grid$long,lat=env_grid$lat)
  
  if(type == "pred"){
    grd.df$p <- pred.grd
  }
  else if(type == "prob")
  {
    grd.df$p <- prob.grd
  }
  
  coordinates(grd.df)<- ~ long + lat
  proj4string(grd.df) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  rast <- rasterize(grd.df,r,field="p",fun=mean)
  return(rast)
  
}


plot_species_raster <- function(rast,title,type="pred",col="#00A600FF",presence=NULL,absence=NULL){
  print(type)
  if(type=="pred"){
    plot(rast,col=c("gray87",col),legend=FALSE, xaxt='n', yaxt='n', ann=FALSE)
    if(!is.null(absence)) {points(absence,cex=.2,col="red")}
  }
  else if(type=="prob"){
    col = heat.colors(10)
    plot(rast,col=col,legend=FALSE, xaxt='n', yaxt='n', ann=FALSE)
    if(!is.null(absence)) {points(absence,cex=.2,col="green")}
  }
  
  if(!is.null(presence)) {points(presence,cex=.2,col="blue")}
  
  mtext(title,padj = .5,side=1,ps=28,font=3)
}

## find sisters
#http://blog.phytools.org/2013/02/function-to-get-sisters-of-node-or-tip.html
find_sisters=function(phy,alt_lat=NULL,category=NULL){
  require(ape)
  require(phytools)
  require(phangorn)
  require(dplyr)

  dd<-lapply(1:phy$Nnode+Ntip(phy),function(n,t)
  Descendants(t,n)[[1]],t=phy)
  nodes<-c(1:phy$Nnode+Ntip(phy))[which(sapply(dd,length)==2)]
  sisters<-as.data.frame(t(sapply(nodes,function(n,t)  
    t$tip.label[Descendants(t,n)[[1]]],t=phy)))
  sisters=dplyr::rename(sisters,sp1=V1,sp2=V2)
  
  if(!is.null(alt_lat)){
    if(length(setdiff(sisters$V1,alt_lat$species))>0){
      stop("sister species labels and alt_lat labels don't match, please fix")
    } else if(is.null(category)){
      stop("please provide category column name from alt_lat")
    } else {
      sisters=left_join(sisters,alt_lat[,c("species",category)],by=c("sp1"="species"))
      sisters=dplyr::rename(sisters,sp1_category=eval(parse(text='category')))
      sisters=left_join(sisters,alt_lat[,c("species",category)],by=c("sp2"="species"))
      sisters=dplyr::rename(sisters,sp2_category=eval(parse(text='category')))
    }
  }
  return(sisters)
}


## phylogeny heatmap (matrix only; doesn't work with niolap or dist objects; latter can be converted to matrix)
heatmap_phylo=function(mat,phy,pdf.name=NULL){
  mat=clean_matrix(mat,phy$tip.label)
  phy.dend=phylogram::as.dendrogram(phy)

  # reorder phylogeny
  mat.ord=mat[phy$tip.label,phy$tip.label]
  if(!is.null(pdf.name)){
    pdf(file=file.path(niche_dir,pdf.name),height = 14,width = 14)
    hm=heatmap(mat.ord,Rowv = phy.dend,Colv="Rowv")
    dev.off()
  } else  hm=heatmap(mat.ord,Rowv = phy.dend,Colv="Rowv")

}

