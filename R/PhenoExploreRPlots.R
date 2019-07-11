#' p is a phenotype object
#' @param expr is a flag for whether to plot gene level (TRUE) or pathway activity level
#' @param method is a choice of "PCA" or "BGA"; default is PCA
#' @import made4
#' @import factoextra
#' @import FactoMineR
#' @export plotPhenotypePCA
plotPhenotypePCA <- function(pheno, expr=TRUE, method=c("PCA","BGA")){
  if (expr){
    m <- pheno@expr
  } else {
    m <- pheno@activity
  }
  colors <- pheno@colors
  if (is.null(colors)){
    colors <- rainbow(length(pheno@locations))
  }
  names(colors) <- pheno@locations
  if (method == "PCA"){
    d.pca <- PCA(t(m), graph=FALSE)
    p <- fviz_pca_ind(d.pca, label="none", habillage = as.factor(pheno@design), palette=colors, addEllipses = TRUE)
  }
  if (method== "BGA"){
    d.bga <- bga(m, classvec=pheno@design)
    p <- plot(d.bga, arraycol=colors[levels(d.bga$fac)])
  }
  p
}
