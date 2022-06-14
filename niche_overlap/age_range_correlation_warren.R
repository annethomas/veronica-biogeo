# modified from https://github.com/danlwarren/arc-extensions/blob/master/age.range.correlation.SES.R

age.range.correlation.SES <- function (phy, overlap, tri = "upper", n = 1000)
{
  age <- branching.times(phy)
  ovlap <- overlap
  if (tri == "upper")
    ovlap[lower.tri(ovlap)] <- t(ovlap)[lower.tri(ovlap)]
  if (tri == "lower")
    ovlap[upper.tri(ovlap)] <- t(ovlap)[upper.tri(ovlap)]
  id <- match(phy$tip.label, rownames(ovlap))
  #print(phy$tip.label)
  #print(rownames(ovlap))
  #print(id)
  ovlap <- ovlap[id, id]
  overlap <- sapply(names(age), nested.mean.overlap, phy = phy,
                    olap = ovlap)
  x <- cbind(age, overlap)
  x.lm <- lm(overlap ~ age)
  randomization <- function(phy, o, n, age) {
    id <- sample(seq(along = o[, 1]))
    print(id)
    rownames(o) <- colnames(o) <- colnames(o)[id]
    rand.o <- sapply(names(age), nested.mean.overlap, phy = phy,
                     olap = o)
    o <- cbind(age, rand.o)
    o <- lm(o[, 2] ~ age)
    list(coef = o$coefficients, fitted.values = o$fitted.values, random.overlap = rand.o)
  }
  ### MODIFIED BIT ###
  random.x <-array(dim=c(n,2))
  random.overlap <-array(dim=c(n,length(age)))
  
  for (i in 1:n){
    print(paste(i,"of",n,"iterations",sep=" "))
    print(Sys.time())
    thisrand  <- randomization(o = ovlap, phy = phy,age = age)
    random.x[i,] <- thisrand$coef
    random.overlap[i,] <- thisrand$random.overlap
  }
  ### END MODIFIED BIT ###
  
  f.intercept <- length(which(random.x[,1] > x.lm$coefficients[1]))/n
  f.slope <- length(which(random.x[,2] > x.lm$coefficients[2]))/n
  f <- c(f.intercept, f.slope)
  p <- sapply(f, function(x) 2 * min(x, 1 - x))
  sig <- cbind(f, p)
  rownames(sig) <- c("intercept", "slope")
  ## more modification:
  sds <- apply(random.overlap, 2, sd)
  means <- apply(random.overlap, 2, mean)
  meds <- apply(random.overlap, 2, median) ## added by AT
  print(means)
  print(sds)
  print(x[,"overlap"])
  ses <- (x[,"overlap"] - means)/sds
  random.overlap <- cbind(age, t(random.overlap))
  colnames(random.overlap) <- c("age", paste("rep", seq(1,n)))
  list(age.range.correlation = x, linear.regression = x.lm,
       sig = sig, sim.overlaps = random.overlap,
       MonteCarlo.replicates = t(random.x), standard.effects = ses, means = means, medians = meds, sds = sds)
}

## from internal phyloclim:
nested.mean.overlap <- function(phy, node, olap){
  
  # match ordering of phy and olap
  # ------------------------------
  id <- match(phy$tip.label, rownames(olap))
  olap <- olap[id, id]
  
  # get daughter nodes
  # ------------------
  d2 <- phy$edge[phy$edge[, 1] == node, 2]
  
  # get descendents of both daughter nodes
  # --------------------------------------
  C1 <- descendants(phy, d2[1])
  C2 <- descendants(phy, d2[2])
  
  # calculate mean overlap
  # ----------------------
  o <- 0
  for (j in C1){
    for (k in C2){
      n <- nbConnectingNodes(phy, c(j, k))
      o <- o + 0.5 ^ (n - 1) * olap[j, k]
    }
  }
  o
}

descendants <-
  function(tree, node, internal = FALSE, string = FALSE){
    
    tips <- seq(along = tree$tip.label)
    x <- tree$edge[,2][tree$edge[,1] == node]
    repeat{
      xx <- x
      x <- sort(unique(c(x, tree$edge[,2][tree$edge[,1] %in% x])))
      if (identical(x, xx)) break
    }
    # return tip number if input is tip number:
    # -----------------------------------------
    if (length(x) == 0) x <- node
    if (!internal)
      x <- x[x %in% tips]
    if (string)
      x <- tree$tip.label[x]
    x
  }

nbConnectingNodes <- function(phy, npair){
  ntips <- length(phy$tip.label)
  nds <- getMRCA(phy, npair)
  nds <- descendants(phy, nds, internal = TRUE)
  if (identical(sort(nds), sort(npair)))
    nb <- 1										else {
      nds <- nds[nds > ntips]
      check <- function(x, npair)
        any(npair %in% descendants(phy, x))
      id <- sapply(nds, check, npair = npair)
      nds <- nds[id] 
      nb <- length(nds) + 1 
    }
  nb
}