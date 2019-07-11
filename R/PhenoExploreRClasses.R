#' Phenotype
#' @export phenotype
phenotype <- setClass("phenotype", slots=c("expr", "design", "motif", "pathways", "activity", "locations", "colors"))

#' @export walk
walk <- setClass("walk", slots=c("x", "y", "d"))

#' @export phenoWanderer
phenoWanderer <- setClass("phenoWanderer", slots=c(w="list"))

#' @export phenoExplorer
phenoExplorer <- setClass("phenoExplorer", slots=c("goal"), contains="phenoWanderer")
