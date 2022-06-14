plot.phylo.HPD.fixed=function (x, nodes = NULL, pb = FALSE, bar.width = 0.3, bar.col = NA, 
          border = NULL, at = NULL, minor = NULL, vline = TRUE, ...) 
{
  op <- par(no.readonly = TRUE)
  plot(x, plot = F, ...)
  ppenv <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  N <- length(x$tip.label)
  yycrds <- ppenv$yy[(N + 1):length(ppenv$yy)]
  yycrdsu <- yycrds + bar.width
  yycrdsl <- yycrds - bar.width
  maxxlim <- max(ppenv$x.lim)
  xxmax <- max(ppenv$xx)
  if (pb) {
    hpd <- sub("_", ",", x$node.label)
    hpd <- matrix(as.numeric(unlist(strsplit(hpd, ","))), 
                  ncol = 2, byrow = T)
    hpd <- rbind(matrix(0, ncol = 2, nrow = length(x$tip.label)), 
                 hpd)
    hpd <- hpd[, c(2, 1)]
  }
  else {
    if (is.null(x$metadata)) 
      stop("phylo object was not read by read.annot.beast; try with pb=T, perhaps?")
    hpd <- x$metadata[, "height_95%_HPD"]
    hpd <- matrix(as.numeric(unlist(strsplit(hpd, ","))), 
                  ncol = 2, byrow = T)
  }
  hpdl <- hpd[(N + 1):dim(hpd)[1], 1]
  hpdu <- hpd[(N + 1):dim(hpd)[1], 2]
  xxu <- -(hpdu - xxmax)
  xxl <- -(hpdl - xxmax)
  par(op, new=TRUE) ## added new=TRUE to prevent blank plot in pdf http://blog.phytools.org/2017/09/possible-solution-for-functions-that.html
  plot(x, x.lim = c(-min(xxl, na.rm = TRUE), maxxlim), ...) ## added negative to min(xxl) to prevent plotting past margin
  if (!is.null(nodes)) {
    rect(xxu[nodes - N], yycrdsl[nodes - N], xxl[nodes - 
                                                   N], yycrdsu[nodes - N], border = border, col = bar.col)
  }
  else {
    rect(xxu, yycrdsl, xxl, yycrdsu, border = border, col = bar.col)
  }
  ax <- simpleAxisPhylo(at = at, minor = minor)
  if (vline) 
    abline(v = ax, lty = 2, lwd = 0.5)
}
