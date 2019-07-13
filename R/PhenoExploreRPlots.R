#' prune matrix by keeping only steps at given frequency
prune <- function(m, samplefreq){
  m[,seq(1:(ncol(m)/samplefreq))*samplefreq]
}

#' plot phenotype trajectories as PCA/BGA
#' @param pheno is a phenotype object
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

#' plots walks of select individual genes
#' @param x is a matrix of gene walks, with genes on rows and time steps on columns
#' @param indrange specifies which genes to plot
#' @param ncol specifies how many columns to wrap on
#' @param scales specifies whether to have "fixed" or "free" scale
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @export genesWalk
genesWalk <- function(x, indrange=1:10, ncol=5, scales="free"){
  xm <- melt(x[indrange,])
  colnames(xm) <- c("gene","time","expression")
  xm %>% ggplot(aes(x = time, y = expression)) +
    geom_line(size = 0.8) +
    facet_wrap(~ gene, ncol=ncol, scales=scales) +
    theme_minimal()
}

#' plots walks of select genes in the same plot
#' @param x is a matrix of gene walks with genes on rows and time steps on columns
#' @param indrange specifies which genes to plot
#' @param ncol specifies how many columns to wrap on
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @export allGenesWalk
allGenesWalk <- function(x, indrange=1:10, ymin=-40, ymax=40){
  xm <- melt(x[indrange,])
  colnames(xm) <- c("gene","time","expression")
  xm %>% ggplot(aes(x = time, y = expression, colour=gene)) +
    geom_line(size = 0.5, show.legend = FALSE) +
    ylim(ymin, ymax) +
    theme_minimal()
}

#' plots walks of select indivdual pathways
#' @param y is a matrix of pathway walks with pathways on rows and time steps on columns
#' @param indrange is the indices of pathways to plot
#' @param bounds is the quantile of target destination
#' @param ncol specifies how many columns to wrap on
#' @param scales whether to have "fixed" or "free" scales
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @export pathwaysWalk
pathwaysWalk <- function(y, indrange=1:10, bounds, ncol=5, scales="free"){
  ym <- data.frame(melt(y[indrange,]), stringsAsFactors = FALSE)
  colnames(ym) <- c("pathway","time","activity")
  colnames(bounds) <- c("lb","ub")
  g <- data.frame(pathway=rownames(bounds),bounds)
  ym <- merge(ym, g)

  ym %>% ggplot(aes(x = time, y = activity)) +
    geom_line(size = 0.8) +
    geom_line(aes(x=time, y=lb), color="blue", linetype="longdash") +
    geom_line(aes(x=time, y=ub), color="red", linetype="longdash") +
    facet_wrap(~ pathway, scales=scales,
               labeller = labeller(pathway = label_wrap_gen(10)), ncol=ncol) +
    theme_minimal()
}


#' plots walks of select pathways in the same plot
#' @param y is a matrix of pathway walks with pathways on rows and time steps on columns
#' @param indrange is the indices of pathways to plot
#' @param ncol specifies how many columns to wrap on
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @export allPathwaysWalk
allPathwaysWalk <- function(y, indrange=1:10, ymin=-1, ymax=1){
  ym <- melt(y[indrange,])
  colnames(ym) <- c("pathway","time","activity")
  ym %>% ggplot(aes(x = time, y = activity, colour=pathway)) +
    geom_line(size = 0.5, show.legend = FALSE) +
    ylim(ymin, ymax) +
    theme_minimal()
}


#' plots distribution walks of select pathways
#' @param w is a list of pathway walks matrices
#' @param startLoc specifies which is the start location
#' @param destLoc specifies target location
#' @param pheno specifies the phenotype object
#' @param indrange specifies the indices of pathways to plot
#' @param samplefreq specifies at what frequency to sample the walks
#' @param ncol specifies how many columns to wrap on
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' @import ggridges
#' @export distributionsWalk
distributionsWalk <- function(w, startLoc, destLoc, pheno, indrange, samplefreq=50, ncol=5){
  yl <- lapply(w, function(ww) ww@y)
  yl <- lapply(yl, function(y) prune(y[indrange,],samplefreq))
  ym <- data.frame(melt(yl)[,1:3], stringsAsFactors = FALSE)
  colnames(ym) <- c("pathway", "time", "activity")
  ym$location <- startLoc

  dest.gsva <- pheno@activity[indrange, which(pheno@design==destLoc)]
  dm <- data.frame(melt(dest.gsva)[,c(1,3)], stringsAsFactors = FALSE)
  colnames(dm) <- c("pathway", "activity")
  dm$location <- destLoc
  dm$time <- rep(max(ym$time), nrow(dm))
  #dm <- Reduce(rbind, lapply(unique(ym$time), function(x) cbind(dm, time=rep(x, length(indrange)))))

  df <- rbind(ym, dm[,colnames(ym)])
  locations <- c(startLoc, destLoc)
  df$location <- factor(df$location, levels=locations)
  df$time <- df$time*samplefreq
  df$time <- factor(df$time)


  p <- ggplot(data=df, aes(x=activity, y=time))+
    geom_density_ridges(aes(x=activity, fill=location), alpha=0.8, color="white") +
    facet_wrap(~pathway, labeller = labeller(pathway = label_wrap_gen(10)), ncol=ncol) +
    theme_minimal() +
    scale_x_continuous() +
    scale_fill_manual(labels=locations,
                      values=unname(pheno@colors[locations]))
  p
}


#' plot distances
#' @param ds is a list of distances
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @export allDistancesWalk
allDistancesWalk <- function(ds, ymin=0, ymax=35){
  dmat <- as.matrix(Reduce(cbind, ds))
  colnames(dmat) <- seq(1:ncol(dmat))

  dmat <- cbind(dmat, mean=apply(dmat, 1, mean))
  dmat_long <- melt(dmat)
  dmat_long$fill <- "individual"
  dmat_long$fill[which(dmat_long$Var2 =="mean")] <- "mean"
  colnames(dmat_long) <- c("step", "pathway", "distance","stat")

  p <- ggplot(data=dmat_long, aes(x=step, y=distance, colour=stat)) +
    geom_line(size=1) +
    theme_minimal() +
    scale_colour_manual(labels=c("individual","mean"), values=c("slategray1","steelblue")) +
    ylim(ymin,ymax)
  p
}
