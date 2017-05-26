makeMADIIfromFieldMap <- function(fieldMap){
  madII <- data.frame(Entry=fieldMap$cname.final, Range=fieldMap$Range, Column=fieldMap$Column)
  return(madII)
}

# Generate something like AR1 error
# The wikipedia entry on AR1 (Autoregressive_model) gives 
# var = \sigma^2_e / (1 - \phi^2) and autocovariance as cov = \phi^dist * \sigma^2_e / (1 - \phi^2)
# This function parameters: field size in terms of numbers of rows and columns; plot size; \phi; \sigma^2_e
# You can supply a MADII design supplied by MADIIdgn instead of nRows and nCols
# plotSize is a two-vector with plotSize[1] the width (along the rows) and plotSize[2] the length (along the columns)
# ar1error returns a vector with residuals that follow the autoregressive model with a variance of 1.
# That is, sigma^2_e = 1 - \phi^2 and so \sigma^2_e does not need to be supplied.
# if a MADII design was supplied, the residuals correspond to the plot order given in the design.
# if nRows and nCols were supplied, the residuals can be formed into a matrix by matrix(ar1error.out, nRows, nCols)

GausCovMat <- function(madII=NULL, nRows=NULL, nCols=NULL, plotSize=c(1, 2), phi){
  if (!is.null(madII)){
    rowVec <- madII$Range
    colVec <- madII$Column
    nRows <- max(rowVec)
    nCols <- max(colVec)
  } else{
    rowVec <- rep(1:nRows, times=nCols)
    colVec <- rep(1:nCols, each=nRows)    
  }
  # Calculate the distance matrix
  rowDistMat <- sapply(rowVec, function(rowNum) abs(rowNum - rowVec))
  colDistMat <- sapply(colVec, function(colNum) abs(colNum - colVec))
  distMat <- sqrt(plotSize[1]^2 * rowDistMat^2 + plotSize[2]^2 * colDistMat^2)
  Gaus <- exp(-(distMat/phi)^2)
  diag(Gaus) <- diag(Gaus)+1e-6
  return(Gaus)
}

GausError <- function(madII=NULL, nRows=NULL, nCols=NULL, plotSize=c(1, 2), phi){
  sqrtCov <- chol(GausCovMat(madII, nRows, nCols, plotSize, phi))
  theErr <- c(t(rnorm(nrow(sqrtCov))) %*% sqrtCov)
  return(theErr / sd(theErr))
}

# NOTE: plot sizes may not be 1 x 1.  In fact, that is not the default.
# fracSpatial is the fraction of the error variation that has spatial correlation.
# The rest is the nuggetVar
simGausFieldTrial <- function(madII, entryEffects=NULL, heritability, fracSpatial, plotSize=c(1, 2), phi){
  # If no entry effects, simulate them
  if (is.null(entryEffects)){
    entries <- sort(unique(madII$Entry)) # Note: if there is a Fill, it will also be given an effect
    entryEffects <- rnorm(length(entries))
    names(entryEffects) <- entries
    errVar <- (1 - heritability) / heritability
  } else{
    errVar <- var(entryEffects) * (1 - heritability) / heritability
  }
  nuggetVar <- errVar * (1 - fracSpatial)
  spatialVar <- errVar * fracSpatial
  spatialEff <- sqrt(spatialVar) * GausError(madII, plotSize=plotSize, phi=phi)
  nuggetEff <- rnorm(nrow(madII), sd=sqrt(nuggetVar))
  fieldObs <- cbind(entryEffects[madII$Entry], spatialEff, nuggetEff)
  rownames(fieldObs) <- NULL
  fieldObs <- data.frame(madII$Entry, madII$Range, madII$Column, fieldObs, apply(fieldObs, 1, sum))
  colnames(fieldObs) <- c("Entry", "Range", "Column", "genoVal", "spatialEff", "nuggetEff", "phenoVal")
  return(list(fieldObs=fieldObs, entryEffects=entryEffects, heritability=heritability, fracSpatial=fracSpatial, plotSize=plotSize))
}

normalDensityPerturbations <- function(madII=NULL, nRows=NULL, nCols=NULL, plotSize=c(1, 2), nPerturb=10, perturbStdDev=sqrt(sum(plotSize^2)), ratioGenPert=0.4){
  if (!is.null(madII)){
    rowVec <- madII$Range
    colVec <- madII$Column
    nRows <- max(rowVec)
    nCols <- max(colVec)
  } else{
    rowVec <- rep(1:nRows, times=nCols)
    colVec <- rep(1:nCols, each=nRows)    
  }
  # Calculate the distance matrix
  rowDistMat <- sapply(rowVec, function(rowNum) abs(rowNum - rowVec))
  colDistMat <- sapply(colVec, function(colNum) abs(colNum - colVec))
  distMat <- sqrt(plotSize[1]^2 * rowDistMat^2 + plotSize[2]^2 * colDistMat^2)
  # Choose plots where the perturbations are centered
  coord <- cbind(sample(nRows, nPerturb, replace=TRUE), sample(nCols, nPerturb, replace=TRUE))
  if (!is.null(madII)){
    plotRows <- apply(coord, 1, function(plotCoord) which(plotCoord[1] == madII$Range & plotCoord[2] == madII$Column))
  } else{
    plotRows <- apply(coord, 1, function(plotCoord) nRows * plotCoord[2] + plotCoord[1])
  }
  # Distribute the perturbation across all the plots
  nPlots <- length(rowVec)
  pertErr <- numeric(nPlots)
  dummy <- sapply(plotRows, function(rowNum){pertErr <<- pertErr + rnorm(1) * dnorm(distMat[rowNum,], sd=runif(1, perturbStdDev/2, 2*perturbStdDev)); return(NULL)})
  pertErr <- pertErr / sd(pertErr)
  totalErr <- rnorm(nPlots, sd=sqrt(ratioGenPert)) + pertErr
  return(totalErr / sd(totalErr))
}

# fracSpatial is the fraction of the error variation that has spatial correlation.
# The rest is the nuggetVar
simNormPertFieldTrial <- function(madII, entryEffects=NULL, heritability=0.5, fracSpatial=0.6, plotSize=c(1, 2), nPerturb=10, perturbStdDev=sqrt(sum(plotSize^2)), ratioGenPert=0.4){
  # If no entry effects, simulate them
  if (is.null(entryEffects)){
    entries <- sort(unique(madII$Entry)) # Note: if there is a Fill, it will also be given an effect
    entryEffects <- rnorm(length(entries))
    names(entryEffects) <- entries
    errVar <- (1 - heritability) / heritability
  } else{
    errVar <- var(entryEffects) * (1 - heritability) / heritability
  }
  nuggetVar <- errVar * (1 - fracSpatial)
  spatialVar <- errVar * fracSpatial
  spatialEff <- sqrt(errVar) * normalDensityPerturbations(madII, plotSize=plotSize, nPerturb=nPerturb, perturbStdDev=perturbStdDev, ratioGenPert=ratioGenPert)
  nuggetEff <- rnorm(nrow(madII), sd=sqrt(nuggetVar))
  fieldObs <- cbind(entryEffects[madII$Entry], spatialEff, nuggetEff)
  rownames(fieldObs) <- NULL
  fieldObs <- data.frame(madII$Entry, madII$Range, madII$Column, fieldObs, apply(fieldObs, 1, sum))
  colnames(fieldObs) <- c("Entry", "Range", "Column", "genoVal", "spatialEff", "nuggetEff", "phenoVal")
  return(list(fieldObs=fieldObs, entryEffects=entryEffects, heritability=heritability, fracSpatial=fracSpatial, plotSize=plotSize))
}

# START HERE
doSim <- FALSE
if (doSim){
  fieldMap <- read.csv("Ibadan_fieldcoord_methodI.csv", stringsAsFactors = F)
  madII <- makeMADIIfromFieldMap(fieldMap)
  ar1.h205.phi05 <- simAR1FieldTrial(madII)
  normPert.h205 <- simNormPertFieldTrial(madII)
  save(ar1.h205.phi05, normPert.h205, file="simForAni.RData")
  
  plotSize <- c(2,1)
  ar1.h205.phi05 <- simAR1FieldTrial(madII, phi=0.6, plotSize=plotSize, fracSpatial=0.5)
  fieldObs <- ar1.h205.phi05$fieldObs
  
  library(ggplot2)
  ggplot(fieldObs, aes(x = Col, y = Row)) + 
    geom_point(mapping = aes(col=fieldObs[,"spatialEff"]), size = 5) +
    scale_color_gradient(low = "red", high = "green") + 
    ggtitle("Spatial Effect")
  ggplot(fieldObs, aes(x = Col, y = Row)) + 
    geom_point(mapping = aes(col=fieldObs[,"nuggetEff"]), size = 5) +
    scale_color_gradient(low = "red", high = "green") +
    ggtitle("Nugget Effect")
  ggplot(fieldObs, aes(x = Col, y = Row)) + 
    geom_point(mapping = aes(col=fieldObs[,"phenoVal"]), size = 5) +
    scale_color_gradient(low = "red", high = "green") +
    ggtitle("Phenotype")
}