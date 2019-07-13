#' Phenotype as pathway activity
#'
#' @import GSVA
pheno <- function(m, pathways, method="ssgsea", kcdf="Gaussian"){
  gsva(m, pathways, method=method, kcdf=kcdf)
}

#' Distance of individual from destination distribution
#'
distance <- function(from, to){
  mlist <- sapply(seq(1:length(from)), function(i)
    ifelse (from[i] >= to[i,1] & from[i] <= to[i,2], 0, min(abs(to[i,1]-from[i]),abs(from[i]-to[i,2]))))
  sum(mlist)
}

#' Mismatch as sum of distances
#'
mismatch <- function(from, to, epsilon=3){
  sum(distance(from, to)) <= epsilon
}

#' Converged pathways
converged <- function(from, to, pathways){
  unique(unname(unlist(pathways[sapply(seq(1:length(from)), function(i) from[i] >= to[i,1] & from[i] <= to[i,2])])))
}

#' Phenotype Wanderer
#'
#' @param expr Gene Expression dataset, can be matrix or data.frame of expression values
#' @param design Binary vector indicating init location information
#' @param motif Regulatory data, can be a binary square matrix or data frame where motif[i,j]=1 indicates TF i regulates gene j
#' @param pathways list of genesets to use of phenotype
#' @param logTransform indicates whether log transformation of expression is done initally
#' @param maxIter numer of steps to simulate
#' @param deltaT time scaling parameter
#' @param alpha step scaling parameter
#' @param outputDir character vector specifying a directory where to save walker results, default is NA
#' @param numMaxCores requires doParallel, foreach.  Runs walker in parallel computing environment.  Set to 1 to avoid parallelization.
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import GSVA
#' @importFrom methods new
#' @return An object of class "phenoWanderer" containing results
#' @examples
#' data(IVY)
#' design <- c(rep(0, 10), rep(NA,20),rep(1,10))
#' phenoWandererRes <- wander(IVY$expr[,0:40],design,IVY$motif,IVY$pathways,maxIter=10, outputDir="simulations",numMaxCores=4)
wander <- function(expr,
                   design,
                   motif,
                   pathways,
                   logTransform=TRUE,
                   maxIter=100,
                   deltaT = 4,
                   alpha=0.5,
                   outputDir=NA,
                   numMaxCores=1){
  # Data type checking
  expr <- checkDataType(expr)

  if(!is.na(numMaxCores)){
    # Calculate the number of cores
    numCores <- detectCores()
    numCores <- min(numCores, numMaxCores)

    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    print(paste(numCores,"cores used"))
  }

  #start time
  strt  <- Sys.time()

  if(!is.na(outputDir)){
    outputDir <- file.path(outputDir, strt)
    dir.create(file.path(outputDir), recursive = TRUE)
  }

  if (logTransform)
    expr <- log2(expr+0.01)
  expr <- t(scale(t(as.matrix(expr))))
  expr[is.nan(expr)] = 0

  # compute the cohort pathway activity scores
  all.gsva <- pheno(expr, pathways)

  # Remove unassigned data
  expr.init <- expr[,design%in%c(0)]
  a.init <- all.gsva[,design%in%c(0)]

  #loop

  # compute initial J
  J0 <- cor(t(expr.init))
  J0[is.na(J0)]=0
  n <- nrow(J0)

  genes <- rownames(motif)

  walks <- foreach(it=1:ncol(expr.init),#ncol(expr.init),
                   .packages=c("GSVA","reshape2")) %dopar% {
                     print(paste0("Running sample ", it))
                     xlist <- c()
                     deltaxlist <- c()
                     ylist <- c()
                     J <- J0
                     X = as.matrix(expr.init[,it])
                     for (j in 1:maxIter){
                       #print(paste0("Walk iteration ", j))
                       # update J
                       y <- pheno(X, pathways)

                       # update X
                       W <- motif * J
                       # compute saturation function Fi
                       Fi <- tanh(X)

                       deltaX <- 1./deltaT * (W %*% Fi)
                       resX <- (1.-alpha/deltaT) * X
                       X <- deltaX + resX

                       # iterate
                       ylist <- cbind(ylist, y)
                       xlist <- cbind(xlist, X)
                     }
                     w <- walk(x=xlist,y=ylist)
                     saveRDS(w,file.path(outputDir,paste0('walk_',colnames(expr.init)[it],'_maxiter_',maxIter,'.rds')))
                     w
                   }

  if(!is.na(numMaxCores)){
    stopCluster(cl)
  }
  print(Sys.time()-strt)
  return(phenoWanderer(w=walks))
}

#' Phenotype Explorer
#'
#' @param expr Gene Expression dataset, can be matrix or data.frame of expression values
#' @param design Binary vector indicating init and dest location information
#' @param motif Regulatory data, can be a binary square matrix or data frame where motif[i,j]=1 indicates TF i regulates gene j
#' @param pathways list of genesets to use of phenotype
#' @param logTransform indicates whether log transformation of expression is done initally
#' @param maxIter numer of steps to simulate
#' @param percent the proportion of the destination distributon (centered around the mean) that should be reached for convergence
#' @param Dsquare noise scaling parameter
#' @param deltaT time scaling parameter
#' @param alpha step scaling parameter
#' @param Jscaling indicates whether a gradient descent is taken: Jscaling 0 < 1 achieves for gradient descent; Jscaling=1 for no scaling
#' @param outputDir character vector specifying a directory where to save walker results, default is NA
#' @param numMaxCores requires doParallel, foreach.  Runs walker in parallel computing environment.  Set to 1 to avoid parallelization.
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import GSVA
#' @importFrom methods new
#' @return An object of class "phenoExplorer" containing results
#' @examples
#' data(IVY)
#' design <- c(rep(0, 10), rep(NA,20),rep(1,10))
#' phenoExplorerRes <- explore(IVY$expr[,0:40],design,IVY$motif,IVY$pathways,maxIter=10, outputDir="simulations",numMaxCores=4)
explore <- function(expr,
                   design,
                   motif,
                   pathways,
                   logTransform=TRUE,
                   maxIter=100,
                   percent=0.8,
                   Dsquare = 0.1, # am not square rooting it below
                   deltaT = 4,
                   alpha=0.5,
                   Jscaling=1,
                   outputDir=NA,
                   numMaxCores=1){

    # Data type checking
    expr <- checkDataType(expr)

    if(!is.na(numMaxCores)){
      # Calculate the number of cores
      numCores <- detectCores()
      numCores <- min(numCores, numMaxCores)

      cl <- makeCluster(numCores)
      registerDoParallel(cl)
      print(paste(numCores,"cores used"))
    }

    #start time
    strt  <- Sys.time()
    if(!is.na(outputDir)){
        outputDir <- file.path(outputDir, strt)
        dir.create(file.path(outputDir), recursive = TRUE)
    }

    #
    if (logTransform)
      expr <- log2(expr+0.01)
    expr <- t(scale(t(as.matrix(expr))))
    expr[is.nan(expr)] = 0

    all.gsva <- pheno(expr, pathways)

    # Remove unassigned data
    expr.init <- expr[,design%in%c(0)]
    a.init <- all.gsva[,design%in%c(0)]
    a.final <- all.gsva[,design%in%c(1)]
    a.quantiles <- t(apply(a.final,1, function(x) unname(quantile(x, probs=c(0.5-percent/2,0.5+percent/2)))))

    #loop

    # compute initial J
    J0 <- cor(t(expr.init))
    J0[is.na(J0)]=0
    n <- nrow(J0)

    genes <- rownames(motif)

    walks <- foreach(it=1:ncol(expr.init),#ncol(expr.init),
                     .packages=c("GSVA","reshape2")) %dopar% {
      print(paste0("Running sample ", it))
      xlist <- c()
      jlist <- list()
      dlist <- c()
      deltaxlist <- c()
      mlist <- c()
      ylist <- c()
      J <- J0
      DeltaJM <- motif
      X = as.matrix(expr.init[,it])
      for (j in 1:maxIter){
        #print(paste0("Walk iteration ", j))
        # update J
        y <- pheno(X, pathways)
        d <- distance(y, a.quantiles)
        dlist <- c(dlist, d)
        m <- mismatch(y, a.quantiles)


        deltaJ <- Dsquare * (1-m) * matrix(rnorm(n*n),n)
        if (Jscaling < 1) {
          clist <- intersect(converged(y, a.quantiles, pathways), genes)
          DeltaJM[clist, ] = Jscaling*DeltaJM[clist,]
          DeltaJM[setdiff(genes,clist),clist] = Jscaling*DeltaJM[setdiff(genes,clist),clist]
          deltaJ <- deltaJ * DeltaJM
        }

        J <- J + deltaJ
        # update X
        W <- motif * J
        # compute saturation function Fi
        Fi <- tanh(X)

        deltaX <- 1./deltaT * (W %*% Fi)
        resX <- (1.-alpha/deltaT) * X
        X <- deltaX + resX

        # iterate
        ylist <- cbind(ylist, y)
        xlist <- cbind(xlist, X)
        mlist <- c(mlist, m)
      }
      w <- walk(x=xlist, y=ylist, d=dlist)
      saveRDS(w,file.path(outputDir,paste0('walk_',colnames(expr.init)[it],'_maxiter_',maxIter,'.rds')))
      w
    }

    if(!is.na(numMaxCores)){
      stopCluster(cl)
    }
    print(Sys.time()-strt)
    return(phenoExplorer(w=walks,goal=a.quantiles))
}

#' Checks that data is in the right format
#'
#' @param expr Gene Expression dataset
#' @return expr Gene Expression dataset in the proper form (may be the same as input)
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' expr.matrix <- matrix(rnorm(100),ncol=15)
#' checkDataType(expr.matrix)
#' #TRUE
#'
checkDataType <- function(expr){
    assert_that(is.data.frame(expr)||is.matrix(expr))
    if(is.data.frame(expr)){
        expr <- as.matrix(expr)
    }
    expr
}
