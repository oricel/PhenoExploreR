#' Design permutation
#' @param init.design a vector of labels with two phenotypes
#' @param permutation.n a number of permutations
#' @param permuation.ratio interval specifiying the percentage of overlap with original design
#' @export
#' @return permutations of labels
permuteDesign <- function(init.design, permutation.n=1000, permutation.ratio=c(0.4, 0.6)){
  ind <- which(!is.na(init.design))
  n <- length(ind)
  new.design <- matrix(NA, nrow=permutation.n * 50, ncol=length(init.design))
  new.design <- unique(t(apply(new.design, 1, function(x) {
    x[ind] <- init.design[ind[sample(1:n, n, replace=FALSE)]]
    x
  })))
  keep <- apply(new.design, 1, function(x)
    length(which(x[ind] == init.design[ind])) >= permutation.ratio[1]*n &
      length(which(x[ind] == init.design[ind])) <= permutation.ratio[2]*n)

  return(rbind(init.design, new.design[which(keep)[1:min(length(which(keep)), permutation.n)],]))
}

#' Pathway Distribution Distance
#' @param m a matrix of pathway activities, rows are pathways, columns samples
#' @param design a vector of labels with phenotypes
#' @param min minimum value for histogram scale; default -1, based on activity values
#' @param max maximum value for histogram scale; default 1, based on activity values
#' @param breaks a number of breaks to use to create pathway distribution histograms
#' @param asGDD whether to return the summary as Global Distribution Distance; summed over PDDs
#' @import HistogramTools
#' @import reshape2
#' @export
#' @return Pathway Distribution Distance
#'
computePDD <- function(m, design, min=-1,max=1, breaks=50, asGDD = FALSE){
  locations <- unique(design)
  comparisons <- combn(locations, 2)
  design.m <- t(apply(comparisons, 2, function(x) {
    d <- design
    d[!(d%in%x)] <- NA
    d
  }))
  dist.list <- lapply(seq(1:nrow(design.m)), function(i){
    dist <- apply(m, 1, function(mm){
      dd <- design.m[i,]
      comparison <- setdiff(unique(dd), NA)
      mm1 <- mm[dd==comparison[1]]
      mm2 <- mm[dd==comparison[2]]
      h1 <- hist(mm1, breaks=seq(min,max,l=breaks), plot=FALSE)
      h2 <- hist(mm2, breaks=seq(min,max,l=breaks), plot=FALSE)
      minkowski.dist(h1,h2,p=2)
    })
  })
  dist.df <- Reduce(cbind, dist.list)
  colnames(dist.df) <- apply(comparisons, 2, function(x) paste0(x[1],"_",x[2]))
  if (asGDD){
    return(apply(dist.df, 2, sum))
  }
  return(dist.df)
}


