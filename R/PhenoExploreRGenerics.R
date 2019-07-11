prune <- function(m, samplefreq){
  seq(1:(ncol(m)/samplefreq))*samplefreq
}

#' pw is a phenoWanderer object
#' @export getWalks
getWalks <- function(pw){
  pw@w
}


#' pw is a phenoWanderer object
#' @export getGeneWalks
getGeneWalks <- function(pw){
  lapply(getWalks(pw), function(w){
    w@x
  })
}
#' pw is a phenoWanderer object
#' @export getPathwayWalks
getPathwayWalks <- function(pw){
  lapply(getWalks(pw), function(w){
    w@y
  })
}

#' pw is a phenoExplorer object
#' @export getDistanceWalks
getDistanceWalks <- function(pw){
  lapply(getWalks(pw), function(w){
    w@d
  })
}


