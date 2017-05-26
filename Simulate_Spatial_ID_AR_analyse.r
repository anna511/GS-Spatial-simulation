###Simulation and analysis of simulated data##
##Simulation code from Jean-Luc Jannink ##December, 2015###
##Ani A. Elias ##December, 2015## modified on May, 2016 ##


##Simulation code to run in loop for phi vector,heritability, and fractSpatial...
##get values for base, model1, and model2 
##plot the correlaiton and mse values of all simulations in a single plot


library(nlme) # to fit the first base model
library(LMERConvenienceFunctions) # remove outliers 
library(regress) #to fit models 
library(reshape2) #to melt matrix
library(lattice) # to visualie cor and err matrix
library(latticeExtra) #add grid lines to plot
library(ggplot2) #for heatmap

#use simulated data madII
#wrapper function to retrieve effects from simulated data
#remove outliers (2.5 times the sd) before proceeding to analysis

#this functon assumes that field coordinates are used for spatial analysis. It is required to change the 
  #method of distance calculation etc. if geocoordinates are used. 
#assuming that range is on y axis and column on x axix as field coordinates; so for AR1, Range is the row and column is the column
#plotSize is a vector of plotdimensions, plotSize[1] is the length of plot in the range direction, 
         #plotSize[2] is width of plot in the column direction
          # doesn't really matter if only isotropic correlation is used
#madII is a matrix from real data - contain genotype name and field coordinates
#phi is the true autocorrelation value to be used for simulating data - a vector of phi values
  # phi  <- seq(0.2,0.8,0.2)
#fracSpatial is the fraction of spatial effect value
#heritability is the value for heritability
#repl indicates status of replication - 1 for only checks replicated, 50 and 100 for 50% and 100% genotypes replicated respectively


geno_spatial_simulate_AR1 <- function(madII, phi, fracSpatial, heritability,plotSize, repl){


entry.cor <- matrix(NA, 20, 10)
colnames(entry.cor) <- c("phi", "Base", "Exp", "Gaus", "Sph","AR","phi.E","phi.G","phi.S","phi.AR")

entry.mse <- matrix(NA, 20, 6)
colnames(entry.mse) <- c("phi", "Base", "Exp", "Gaus", "Sph", "AR")

sp.cor <- matrix(NA, 20, 9)
colnames(sp.cor) <- c("phi", "sp_Exp", "sp_Gaus","sp_Sph","sp_AR", "phi.E","phi.G","phi.S","phi.AR")

sp.mse <- matrix(NA, 20, 5)
colnames(sp.mse) <- c("phi", "sp_Exp", "sp_Gaus","sp_Sph","sp_AR")

resid.cor <- matrix(NA, 20, 5)
colnames(resid.cor) <- c("phi", "res_Exp", "res_Gaus","res_Sph","res_AR")

resid.mse <- matrix(NA, 20, 5)
colnames(resid.mse) <- c("phi", "res_Exp","res_Gaus","res_Sph","res_AR")

fra.Sp <- matrix(NA, 20, 6)
colnames(fra.Sp) <- c("phi", "fra.Sp", "FrSp_Exp","FrSp_Gaus","FrSp_Sph","FrSp_AR")

herit <- matrix(NA,20,6)
colnames(herit) <- c("phi", "h2", "h2_Exp","h2_Gaus","h2_Sph","h2_AR")

model.fail <- matrix(NA,20,4)
colnames(model.fail) <- c("fail_Exp","fail_Gaus","fail_Sph","fail_AR")



entry.cor.2 <- matrix(NA,length(phi),10)
colnames(entry.cor.2) <- c("phi", "Base", "Exp", "Gaus", "Sph","AR","phi.E","phi.G","phi.S","phi.AR")

entry.mse.2 <- matrix(NA,length(phi),6)
colnames(entry.mse.2) <- c("phi", "Base", "Exp", "Gaus", "Sph", "AR")

sp.cor.2 <- matrix(NA, length(phi), 9)
colnames(sp.cor.2) <- c("phi", "sp_Exp", "sp_Gaus","sp_Sph","sp_AR", "phi.E","phi.G","phi.S","phi.AR")

sp.mse.2 <- matrix(NA, length(phi), 5)
colnames(sp.mse.2) <- c("phi", "sp_Exp", "sp_Gaus","sp_Sph","sp_AR")

resid.cor.2 <- matrix(NA, length(phi), 5)
colnames(resid.cor.2) <- c("phi", "res_Exp", "res_Gaus","res_Sph","res_AR")

resid.mse.2 <- matrix(NA, length(phi), 5)
colnames(resid.mse.2) <- c("phi", "res_Exp","res_Gaus","res_Sph","res_AR")

fra.Sp.2 <- matrix(NA, length(phi), 6)
colnames(fra.Sp.2) <- c("phi", "fra.Sp", "FrSp_Exp","FrSp_Gaus","FrSp_Sph","FrSp_AR")

herit.2 <- matrix(NA, length(phi), 6)
colnames(herit.2) <- c("phi", "h2", "h2_Exp","h2_Gaus","h2_Sph","h2_AR")

model.fail.2 <- matrix(NA, length(phi), 4)
colnames(model.fail) <- c("fail_Exp","fail_Gaus","fail_Sph","fail_AR")


for(i in c(1:length(phi))){ 
  for(s in c(1:20)){
#create data1
  ar1.h2.phi <- simAR1FieldTrial(madII, phi=phi[i], plotSize=plotSize, fracSpatial=fracSpatial, heritability = heritability)
  data1 <- ar1.h2.phi$fieldObs
   
  fixed <- formula("phenoVal ~ 1")
  base <- lme(fixed = fixed, random = ~ 1|Entry, data = na.omit(data1), method = "ML")
  data2 <- (romr.fnc(base, na.omit(data1), trim = 2.5))$data #outlier removal

#re-running the model using dataset with outliers removed
  base <- try(regress(fixed, ~ as.factor(data2[,"Entry"]), pos=rep(TRUE,2), data = data2),silent=TRUE)
  fixed.effect <- base$formula
  
#calculate the distance matrices
  #modified field coordinates and distance matrices
  rowVec<- data2$Range
  colVec <- data2$Column

  rowVec2 <- rowVec * plotSize[1] 
  rowDistMat2 <- sapply(rowVec2, function(rowNum) abs(rowNum - rowVec2))

  colVec2 <- colVec * plotSize[2]
  colDistMat2 <- sapply(colVec2, function(colNum) abs(colNum - colVec2))

  distMat <- sqrt(rowDistMat2^2 + colDistMat2^2)

#distance design matrix
  data2$Slno <- c(1:nrow(data2))
  rownames(distMat) <- colnames(distMat) <- data2$Slno

  idDist <- factor(as.character(data2$Slno),levels =rownames(distMat))
    Z.spDist <- model.matrix(~idDist-1)  

#model 1 - add correlation structures to selected base model
    #Isotropic correlation structures  
 phi.1 <- seq(0.5,max(distMat), 6)
  model1.E.E <- model1.E.fail <- vector()
  model1.G.E <- model1.G.fail <- vector()
  model1.S.E <- model1.S.fail <- vector()
  
  for (j in 1:length(phi.1)){
    #correlation structure
    spE <-  exp(-(distMat/phi.1[j]))
    
    spG <- exp(-(distMat/phi.1[j])^2)
      
    distMat.s <- distMat
    distMat.s[distMat.s > phi.1[j]] <- NA     
    spS <- 1-(1.5*(distMat.s/phi.1[j])) + 0.5*((distMat.s/phi.1[j])^3)
    spS[is.na(spS)] <- 0 

    spE.R <- Z.spDist%*%spE%*%t(Z.spDist)
    spG.R <- Z.spDist%*%spG%*%t(Z.spDist)
    spS.R <- Z.spDist%*%spS%*%t(Z.spDist)

    #check models and select
    model1.E <- try(regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spE.R, pos= rep(TRUE,3), tol = 1e-4, data = data2),silent = TRUE)
    if (class(model1.E) != "try-error") {
      model1.E.E[j] <- mean(abs((data2[,"genoVal"] + data2[,"spatialEff"]) - model1.E$predicted))
    }else{
      model1.E.fail[j] <- 1
    }
    model1.G <- try(regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spG.R, pos= rep(TRUE,3), tol = 1e-4, data = data2),silent = TRUE)
    if (class(model1.G) != "try-error") {
      model1.G.E[j] <- mean(abs((data2[,"genoVal"] + data2[,"spatialEff"]) - model1.G$predicted))
    }else{
      model1.G.fail[j] <- 1
    }
    model1.S <- try(regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spS.R, pos= rep(TRUE,3), tol = 1e-4, data = data2),silent = TRUE)
    if (class(model1.S) != "try-error") {
      model1.S.E[j] <- mean(abs((data2[,"genoVal"] + data2[,"spatialEff"]) - model1.S$predicted))
    }else{
      model1.S.fail[j] <- 1
    }
  }
  
  model1.E.min <- which(model1.E.E == min(model1.E.E, na.rm = TRUE))
  model1.E.failed <- sum(model1.E.fail, na.rm = TRUE)
  model1.G.min <- which(model1.G.E == min(model1.G.E, na.rm = TRUE))
  model1.G.failed <- sum(model1.G.fail, na.rm = TRUE)
  model1.S.min <- which(model1.S.E == min(model1.S.E, na.rm = TRUE))
  model1.S.failed <- sum(model1.S.fail, na.rm = TRUE)

  phi.1.E <- phi.1[model1.E.min]
  phi.1.G <- phi.1[model1.G.min]
  phi.1.S <- phi.1[model1.S.min]

  #correlation structures
  spE <-  exp(-(distMat/phi.1.E))
    
  spG <- exp(-(distMat/phi.1.G)^2)
      
  distMat.s <- distMat
  distMat.s[distMat.s > phi.1.S] <- NA     
  spS <- 1-(1.5*(distMat.s/phi.1.S)) + 0.5*((distMat.s/phi.1.S)^3)
  spS[is.na(spS)] <- 0 

  spE.R <- Z.spDist%*%spE%*%t(Z.spDist)
  spG.R <- Z.spDist%*%spG%*%t(Z.spDist)
  spS.R <- Z.spDist%*%spS%*%t(Z.spDist)

  model1.E <- regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spE.R, pos= rep(TRUE,3), tol = 1e-4, data = data2)
  model1.G <- regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spG.R, pos= rep(TRUE,3), tol = 1e-4, data = data2)
  model1.S <- regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + spS.R, pos= rep(TRUE,3), tol = 1e-4, data = data2)
   


#model2 - AR1 model
  phi.2 <- seq(0.1,0.9,0.1)  
    model2.E <- model2.fail <-  vector() 

  for(j in 1:length(phi.2)){
    rcAR <- phi.2[j]^distMat

    rcAR.R <- Z.spDist%*%rcAR%*%t(Z.spDist)

    model2 <- try(regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + rcAR.R, pos=rep(TRUE,3),tol = 1e-4, data = data2), silent = TRUE)
    if (class(model2) != "try-error") {
    model2.E[j] <- mean(abs((data2[,"genoVal"] + data2[,"spatialEff"]) - model2$predicted))
    }else {
      model2.fail[j] <- 1
    }
  }

  model2.min <- which(model2.E == min(model2.E, na.rm = TRUE)) 
  model2.failed <- sum(model2.fail, na.rm = TRUE)
 
  phi.2a <- phi.2[model2.min]
  rcAR <- phi.2a^distMat
  rcAR.R <- Z.spDist%*%rcAR%*%t(Z.spDist)

  model2 <- regress(fixed.effect, ~ as.factor(data2[,"Entry"]) + rcAR.R, pos = rep(TRUE,3), tol = 1e-4, data = data2)

  

 ##blup and correlation##
  #blup from base and models
  baseBLUP <- BLUP(base)$Mean
  model1.E.BLUP <- BLUP(model1.E)
  model1.G.BLUP <- BLUP(model1.G)
  model1.S.BLUP <- BLUP(model1.S) 
  model2.BLUP <- BLUP(model2)
  

  #separating spatial and entry blups for models
   model1.E.spBLUP <- model1.E.BLUP$Mean[grep("spE", names(model1.E.BLUP$Mean), fixed=TRUE)]
   model1.G.spBLUP <- model1.G.BLUP$Mean[grep("spG", names(model1.G.BLUP$Mean), fixed=TRUE)]
   model1.S.spBLUP <- model1.S.BLUP$Mean[grep("spS", names(model1.S.BLUP$Mean), fixed=TRUE)]
   model2.spBLUP <- model2.BLUP$Mean[grep("rcAR", names(model2.BLUP$Mean), fixed=TRUE)]
  
   data2.sp <- cbind(data2, model1.E.spBLUP, model1.G.spBLUP,model1.S.spBLUP, model2.spBLUP)
    

  base.BLUP <- BLUP(base)$Mean
  model1.E.EB <- model1.E.BLUP$Mean[grep("Entry", names(model1.E.BLUP$Mean), fixed=TRUE)]
    namesEB <- names(model1.E.BLUP$Mean)[grep("Entry", names(model1.E.BLUP$Mean), fixed=TRUE)]
    namesEB <- matrix(unlist(strsplit(namesEB, ".", fixed=TRUE)), 3)[3,]
    names(model1.E.EB) <- namesEB  
  model1.G.EB <- model1.G.BLUP$Mean[grep("Entry", names(model1.G.BLUP$Mean), fixed=TRUE)]
    namesEB <- names(model1.G.BLUP$Mean)[grep("Entry", names(model1.G.BLUP$Mean), fixed=TRUE)]
    namesEB <- matrix(unlist(strsplit(namesEB, ".", fixed=TRUE)), 3)[3,]
    names(model1.G.EB) <- namesEB  
  model1.S.EB <- model1.S.BLUP$Mean[grep("Entry", names(model1.S.BLUP$Mean), fixed=TRUE)]
    namesEB <- names(model1.S.BLUP$Mean)[grep("Entry", names(model1.S.BLUP$Mean), fixed=TRUE)]
    namesEB <- matrix(unlist(strsplit(namesEB, ".", fixed=TRUE)), 3)[3,]
    names(model1.S.EB) <- namesEB  
  model2.EB <- model2.BLUP$Mean[grep("Entry", names(model2.BLUP$Mean), fixed=TRUE)]
    namesEB <- names(model2.BLUP$Mean)[grep("Entry", names(model2.BLUP$Mean), fixed=TRUE)]
    namesEB <- matrix(unlist(strsplit(namesEB, ".", fixed=TRUE)), 3)[3,]
    names(model2.EB) <- namesEB
  
  entryEff <- ar1.h2.phi$entryEffects
  
  data2.entry <- cbind(entryEff, base.BLUP, model1.E.EB, model1.G.EB, model1.S.EB, model2.EB)
    

  #Correlation between baseBLUP and entry BLUPs from models
  
  
  entry.cor[s,1] <- phi[i]
  entry.cor[s,2] <- cor(entryEff, base.BLUP)
  entry.cor[s,3] <- cor(entryEff, model1.E.EB)
  entry.cor[s,4] <- cor(entryEff, model1.G.EB)
  entry.cor[s,5] <- cor(entryEff, model1.S.EB)
  entry.cor[s,6] <- cor(entryEff, model2.EB)
  entry.cor[s,7] <- phi.1.E
  entry.cor[s,8] <- phi.1.G
  entry.cor[s,9] <- phi.1.S  
  entry.cor[s,10] <- phi.2a  

  
  entry.mse[s,1] <- phi[i]
  entry.mse[s,2] <- mean(abs(entryEff - base.BLUP))
  entry.mse[s,3] <- mean(abs(entryEff - model1.E.EB))  
  entry.mse[s,4] <- mean(abs(entryEff - model1.G.EB))
  entry.mse[s,5] <- mean(abs(entryEff - model1.S.EB))
  entry.mse[s,6] <- mean(abs(entryEff - model2.EB))
  


  #Correlation between simulated and model spatial effects
 
  sp.cor[s,1] <- phi[i]
  sp.cor[s,2] <- cor(data2$spatialEff, model1.E.spBLUP)
  sp.cor[s,3] <- cor(data2$spatialEff, model1.G.spBLUP)
  sp.cor[s,4] <- cor(data2$spatialEff, model1.S.spBLUP)
  sp.cor[s,5] <- cor(data2$spatialEff, model2.spBLUP)
  sp.cor[s,6] <- phi.1.E
  sp.cor[s,7] <- phi.1.G
  sp.cor[s,8] <- phi.1.S   
  sp.cor[s,9] <- phi.2a  

  
  sp.mse[s,1] <- phi[i]
  sp.mse[s,2] <- mean(abs(data2$spatialEff - model1.E.spBLUP))
  sp.mse[s,3] <- mean(abs(data2$spatialEff - model1.G.spBLUP))
  sp.mse[s,4] <- mean(abs(data2$spatialEff - model1.S.spBLUP))
  sp.mse[s,5] <- mean(abs(data2$spatialEff - model2.spBLUP))  
  

  #Correlation between simulated and model residual error
  
  resid.cor[s,1] <- phi[i]
  resid.cor[s,2] <- cor(data2$nuggetEff, (data2[,"phenoVal"] - model1.E$predicted))
  resid.cor[s,3] <- cor(data2$nuggetEff, (data2[,"phenoVal"] - model1.G$predicted))
  resid.cor[s,4] <- cor(data2$nuggetEff, (data2[,"phenoVal"] - model1.S$predicted))
  resid.cor[s,5] <- cor(data2$nuggetEff, (data2[,"phenoVal"] - model2$predicted))
  
  
  resid.mse[s,1] <- phi[i]
  resid.mse[s,2] <- mean(abs(data2$nuggetEff - (data2[,"phenoVal"] - model1.E$predicted)))
  resid.mse[s,3] <- mean(abs(data2$nuggetEff - (data2[,"phenoVal"] - model1.G$predicted)))
  resid.mse[s,4] <- mean(abs(data2$nuggetEff - (data2[,"phenoVal"] - model1.S$predicted)))
  resid.mse[s,5] <- mean(abs(data2$nuggetEff - (data2[,"phenoVal"] - model2$predicted)))
  

  #Fraction of spatial error

  fra.Sp[s,1] <- phi[i]
  fra.Sp[s,2] <- fracSpatial
  fra.Sp[s,3] <- var(model1.E.spBLUP) /(var(model1.E.spBLUP) + var(data2[,"phenoVal"] - model1.E$predicted))
  fra.Sp[s,4] <- var(model1.G.spBLUP) /(var(model1.G.spBLUP) + var(data2[,"phenoVal"] - model1.G$predicted))
  fra.Sp[s,5] <- var(model1.S.spBLUP) /(var(model1.S.spBLUP) + var(data2[,"phenoVal"] - model1.S$predicted))
  fra.Sp[s,6] <- var(model2.spBLUP) /(var(model2.spBLUP) + var(data2[,"phenoVal"] - model2$predicted))
  
  #Heritability

  herit[s,1] <- phi[i]
  herit[s,2] <- heritability
  herit[s,3] <- var(model1.E.EB) /(var(model1.E.EB) + var(model1.E.spBLUP) + var(data2[,"phenoVal"] - model1.E$predicted))
  herit[s,4] <- var(model1.G.EB) /(var(model1.G.EB) + var(model1.G.spBLUP) + var(data2[,"phenoVal"] - model1.G$predicted))
  herit[s,5] <- var(model1.S.EB) /(var(model1.S.EB) + var(model1.S.spBLUP) + var(data2[,"phenoVal"] - model1.S$predicted))
  herit[s,6] <- var(model2.EB) /(var(model2.EB) + var(model2.spBLUP) + var(data2[,"phenoVal"] - model2$predicted)) 

  #Failure

  model.fail[s,1] <- model1.E.failed
  model.fail[s,2] <- model1.G.failed
  model.fail[s,3] <- model1.S.failed
  model.fail[s,4] <- model2.failed 

 } #end of for w.r.t simulation for one phi value 

 #for visualization of cor and mse value per phi value
 entry.cor.melt <- melt(entry.cor)
 entry.mse.melt <- melt(entry.mse)
 sp.cor.melt <- melt(sp.cor)
 sp.mse.melt <- melt(sp.mse)
 resid.cor.melt <- melt(resid.cor)
 resid.mse.melt <- melt(resid.mse)
 fra.Sp.melt <- melt(fra.Sp)
 herit.melt <- melt(herit)

 #output the files simulation for one phi
 write.csv(entry.cor.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_entry_Effect_cor_melt.csv"))
 write.csv(entry.mse.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_entry_Effect_rmse_melt.csv"))
 write.csv(sp.cor.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_spatial_Effect_cor_melt.csv"))
 write.csv(sp.mse.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_spatial_Effect_rmse_melt.csv"))
 write.csv(resid.cor.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_residual_cor_melt.csv"))
 write.csv(resid.mse.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_residual_rmse_melt.csv"))
 write.csv(fra.Sp.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_fration_Spatial_melt.csv"))
 write.csv(herit.melt, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_heritability_melt.csv"))
 write.csv(model.fail, paste0("fracSp_",fracSpatial,"-phi_",phi[i],"-rep_",repl,"-h2",heritability,"_failure_melt.csv"))

  #mean of correlation and mse 
  #entry    
  entry.cor.2[i,] <- apply(entry.cor,2,mean, na.rm = TRUE) 
  entry.mse.2[i,] <- apply(entry.mse,2,mean, na.rm = TRUE)  

  #spatial
  sp.cor.2[i,] <- apply(sp.cor,2,mean, na.rm = TRUE)
  sp.mse.2[i,] <- apply(sp.mse,2,mean, na.rm = TRUE)

  #residual
  resid.cor.2[i,] <- apply(resid.cor, 2, mean, na.rm = TRUE)
  resid.mse.2[i,] <- apply(resid.mse, 2, mean, na.rm = TRUE)

  #fra.spatial
  fra.Sp.2[i,] <- apply(fra.Sp, 2, mean, na.rm = TRUE)

  #heritability
  herit.2[i,] <- apply(herit, 2, mean, na.rm = TRUE) 

  #failure
  model.fail.2[i,] <- apply(model.fail, 2, sum, na.ram = TRUE) 

} # end of for w.r.t. phi

  write.csv(entry.cor.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_entry_Effect_cor.csv"))
  write.csv(entry.mse.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_entry_Effect_error.csv"))
  write.csv(sp.cor.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_spatial_Effect_cor.csv"))
  write.csv(sp.mse.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_spatial_Effect_error.csv"))
  write.csv(resid.cor.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_residual_cor.csv"))
  write.csv(resid.mse.2, paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_residual_error.csv"))
  write.csv(fra.Sp.2,  paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_fraction_Spatial.csv"))
  write.csv(herit.2,  paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_heritability.csv"))
  write.csv(model.fail.2,  paste0("fracSp_",fracSpatial,"-rep_",repl,"-h2",heritability,"_failure.csv"))

} # end of function





