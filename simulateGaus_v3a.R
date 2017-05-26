makeMADIIfromFieldMap <- function(fieldMap){
  madII <- data.frame(Entry=fieldMap$cname.final, Range=fieldMap$Range, Column=fieldMap$Column)
  return(madII)
}

# Generate something like Gaus error
# This function parameters: field size in terms of numbers of rows and columns; plot size; \phi; \sigma^2_e
# You can supply a MADII design supplied by MADIIdgn instead of nRows and nCols
# plotSize is a two-vector with plotSize[1] the width (along the rows) and plotSize[2] the length (along the columns)
# returns a vector with residuals that follow the Gaussian model with a variance of 1.
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

