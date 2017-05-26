###Spatial analysis using K matrix and various spatial variance-covariance structures##
###Ani A. Elias #### March - April, 2016 ####

library(rrBLUP) # to calculate K matrix
library(nlme) # to fit the first base model
library(LMERConvenienceFunctions) # remove outliers 
library(regress) #to fit models 

#this functon assumes that field coordinates are used for spatial analysis. It is required to change the 
  #method of distance calculation etc. if geocoordinates are used. 
#data is expected to have no rows with all NA values 
#genotype is the name of column containing genotypic predictor variable
#trait is a vector of names of traits that needs to be tested - response variable
#K is the relationship matrix from 'rrBLUP' OR snp is the snp data 
 
      #K matrix from snp data if snp data is to be provided and not the K.mat matrix (see below).
  #NOTE: the function is written in such a way that the user provide the K matrix (= K.mat).
        #but if snp data is provided instead, remove K.mat from function, load snp data and run the following:
        #eg. of loading snp data 
          #load(".../V6_IITA_BeagleDosages_V100715_MAF01.Rdata")
        #subsetting snp data - helps to remove clones that are not genotyped but phenotyped 
          #table(data1a[,genotype] %in% rownames(snps))  
          #snp2.names <- intersect(data1a[,genotype], rownames(snps))
          #snps2 <- subset(snps, rownames(snps) %in% snp2.names)

        #calculate K matrix
          #K.mat <- A.mat(snps2-1)
  #NOTE:- A.mat normalizes the K matrix using a normalization constant,c,
          #the mean of the diagonal elements is 1 + f. So, when using A.mat,
          #snps to K has to be done inside the loop because outliers are also 
          #removed from the original data which can change the number of genotypes
          #in the K matrix. 
           #when no normalization is done to K matrix, the K matrix can be supplied directly
           #to the function and the values in K matrix are not affected by 
           #outlier removal. However, the final model selected can be different because
           #K matrix influence the genotypic variance as g ~ N(0,sigma^2 * K), and therefore,
           # the leftover error variance which is partitined into spatial and residual error variances. 
           #One more reason is because the prmse values from CV are not drastically different from each other.


#assuming that range is on x axis and column on y axix as field coordinates; 
  #Range is the row (to be named as 'Range') and column (to be named as 'Column') is the column
  #range in E-W direction and column in N-S direction
  
  #plotSize as a vector of plotdimensions, plotSize[1] is the width of plot 
          #plots share longer edges across the range
          #row distance is the distance between ranges which is the shorter distance, 
          #plotSize[2] is length of plot
          #column distance is the longer distnace = distance between columns
          #in IITA, plots are arranged with the longer side shared across the ranges. 
          #That is, range distance is always shorter than that of column
        #assuming unit distance between adjascent plots in any direction. But, this can be ch
#loc is the name of location  - for easiness in identifying files in output
#kfold is the number of folds for each k-fold CV
#mkfold is the total number of k-fold CVs
#K.mat is the relationship square matrix based on snp or pedigree  whose rownames and colnames are the names of genotypes

#remove outliers (2.5 times the sd) before proceeding to analysis

#modify phi.2 and phi values based on prior information in order to speed up analysis


geno_spatial_snp <- function(data1, trait, genotype, plotSize, loc, kfold, mkfold, K.mat){ 
  
  phi.2 <- seq(0.1,0.9, 0.1)  
  
  for(i in c(1:length(trait))){  

  test.geno <- vector("list",mkfold)  

  base.prmse <- base.pcor <- base.R.prmse <- base.R.pcor <- base.C.prmse <- matrix(NA,mkfold,kfold)  
  base.C.pcor <- base.RC.prmse <- base.RC.pcor <- matrix(NA,mkfold,kfold)
  base.fail <- base.R.fail <- base.C.fail <- base.RC.fail <- matrix(NA, mkfold,kfold)  
  
  data1a <- na.omit(subset(data1,select = c("Range","Column",genotype,trait[i])))
  
  fixed.effect <- formula(paste(trait[i],"1", sep="~"))
  base <- lme(fixed = fixed.effect, random = as.formula(paste0("~ 1|", genotype)), data = na.omit(data1a), method = "ML")
  data1a <- (romr.fnc(base, na.omit(data1a), trim = 2.5))$data #outlier removal

  #K matrix aligning with data1a
  #subsetting K matrix  - helps to remove clones that are not genotyped but phenotyped   
  K.names <- intersect(data1a[,genotype], rownames(K.mat))
  K <- K.mat[K.names, K.names]
  

  data2 <- subset(data1a,data1a[,genotype] %in% K.names)
    write.csv(data2, "data2.csv")
    data2 <- read.csv("data2.csv", h=T, sep=",")   


  #computing the distance matrix
  rowVec<- data2$Range
  colVec <- data2$Column

  rowVec2 <- rowVec * plotSize[1] 
  rowDistMat2 <- sapply(rowVec2, function(rowNum) abs(rowNum - rowVec2))

  colVec2 <- colVec * plotSize[2]
  colDistMat2 <- sapply(colVec2, function(colNum) abs(colNum - colVec2))

  distMat <- sqrt(rowDistMat2^2 + colDistMat2^2)

  data2$Slno <- c(1:nrow(data2))
  rownames(distMat) <- colnames(distMat) <- data2$Slno

#phi and matrices for prmse and pcor
  phi <- seq(0.5,max(distMat), 0.5) 

  model1.E.prmse.1 <- model1.E.pcor.1 <- model1.E.fail.1 <- model1.Er.prmse.1 <- model1.Er.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Er.fail.1 <- model1.Ec.prmse.1 <- model1.Ec.pcor.1 <- model1.Ec.fail.1 <-  matrix(NA,mkfold,length(phi))

  model1.E.R.prmse.1 <- model1.E.R.pcor.1 <- model1.E.R.fail.1 <- model1.Er.R.prmse.1 <- model1.Er.R.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Er.R.fail.1 <- model1.Ec.R.prmse.1 <- model1.Ec.R.pcor.1 <- model1.Ec.R.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.E.C.prmse.1 <- model1.E.C.pcor.1 <- model1.E.C.fail.1 <- model1.Er.C.prmse.1 <- model1.Er.C.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Er.C.fail.1 <- model1.Ec.C.prmse.1 <- model1.Ec.C.pcor.1 <- model1.Ec.C.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.E.RC.prmse.1 <- model1.E.RC.pcor.1 <- model1.E.RC.fail.1 <- model1.Er.RC.prmse.1 <- model1.Er.RC.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Er.RC.fail.1 <- model1.Ec.RC.prmse.1 <- model1.Ec.RC.pcor.1 <- model1.Ec.RC.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.G.prmse.1 <- model1.G.pcor.1 <- model1.G.fail.1 <- model1.Gr.prmse.1 <- model1.Gr.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Gr.fail.1 <- model1.Gc.prmse.1 <- model1.Gc.pcor.1 <- model1.Gc.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.G.R.prmse.1 <- model1.G.R.pcor.1 <- model1.G.R.fail.1 <- model1.Gr.R.prmse.1 <- model1.Gr.R.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Gr.R.fail.1 <- model1.Gc.R.prmse.1 <- model1.Gc.R.pcor.1 <- model1.Gc.R.fail.1 <-  matrix(NA,mkfold,length(phi))

  model1.G.C.prmse.1 <- model1.G.C.pcor.1 <- model1.G.C.fail.1 <- model1.Gr.C.prmse.1 <- model1.Gr.C.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Gr.C.fail.1 <- model1.Gc.C.prmse.1 <- model1.Gc.C.pcor.1 <- model1.Gc.C.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.G.RC.prmse.1 <- model1.G.RC.pcor.1 <- model1.G.RC.fail.1 <- model1.Gr.RC.prmse.1 <- model1.Gr.RC.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Gr.RC.fail.1 <- model1.Gc.RC.prmse.1 <- model1.Gc.RC.pcor.1 <- model1.Gc.RC.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.S.prmse.1 <- model1.S.pcor.1 <- model1.S.fail.1 <- model1.Sr.prmse.1 <- model1.Sr.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Sr.fail.1 <- model1.Sc.prmse.1 <- model1.Sc.pcor.1 <- model1.Sc.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.S.R.prmse.1 <- model1.S.R.pcor.1 <- model1.S.R.fail.1 <- model1.Sr.R.prmse.1 <- model1.Sr.R.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Sr.R.fail.1 <- model1.Sc.R.prmse.1 <- model1.Sc.R.pcor.1 <-  model1.Sc.R.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.S.C.prmse.1 <- model1.S.C.pcor.1 <- model1.S.C.fail.1 <- model1.Sr.C.prmse.1 <- model1.Sr.C.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Sr.C.fail.1 <- model1.Sc.C.prmse.1 <- model1.Sc.C.pcor.1 <- model1.Sc.C.fail.1 <- matrix(NA,mkfold,length(phi))

  model1.S.RC.prmse.1 <- model1.S.RC.pcor.1 <- model1.S.RC.fail.1 <- model1.Sr.RC.prmse.1 <- model1.Sr.RC.pcor.1 <- matrix(NA,mkfold,length(phi))
  model1.Sr.RC.fail.1 <- model1.Sc.RC.prmse.1 <- model1.Sc.RC.pcor.1 <- model1.Sc.RC.fail.1 <- matrix(NA,mkfold,length(phi))

  model2.AR.prmse.1 <- model2.AR.pcor.1 <- model2.AR.fail.1 <- model2.ARr.prmse.1 <- model2.ARr.pcor.1 <- matrix(NA,mkfold,length(phi.2))
  model2.ARr.fail.1 <- model2.ARc.prmse.1 <- model2.ARc.pcor.1 <- model2.ARc.fail.1 <- matrix(NA,mkfold,length(phi.2)) 

  model2.AR.R.prmse.1 <- model2.AR.R.pcor.1 <- model2.AR.R.fail.1 <- model2.ARr.R.prmse.1 <- model2.ARr.R.pcor.1 <- matrix(NA,mkfold,length(phi.2))
  model2.ARr.R.fail.1 <- model2.ARc.R.prmse.1 <- model2.ARc.R.pcor.1 <- model2.ARc.R.fail.1 <- matrix(NA,mkfold,length(phi.2)) 

  model2.AR.C.prmse.1 <- model2.AR.C.pcor.1 <- model2.AR.C.fail.1 <- model2.ARr.C.prmse.1 <- model2.ARr.C.pcor.1 <- matrix(NA,mkfold,length(phi.2))
  model2.ARr.C.fail.1 <- model2.ARc.C.prmse.1 <- model2.ARc.C.pcor.1 <- model2.ARc.C.fail.1 <- matrix(NA,mkfold,length(phi.2)) 

  model2.AR.RC.prmse.1 <- model2.AR.RC.pcor.1 <- model2.AR.RC.fail.1 <- model2.ARr.RC.prmse.1 <- model2.ARr.RC.pcor.1 <- matrix(NA,mkfold,length(phi.2))
  model2.ARr.RC.fail.1 <- model2.ARc.RC.prmse.1 <- model2.ARc.RC.pcor.1 <-  model2.ARc.RC.fail.1 <- matrix(NA,mkfold,length(phi.2)) 


  model1.E.prmse.2 <- model1.E.pcor.2 <- model1.Er.prmse.2 <- model1.Er.pcor.2 <- model1.Ec.prmse.2 <- vector()
  model1.Ec.pcor.2 <-  model1.E.R.prmse.2 <- model1.E.R.pcor.2 <- vector()
  model1.Er.R.prmse.2 <- model1.Er.R.pcor.2 <- model1.Ec.R.prmse.2 <- model1.Ec.R.pcor.2 <-  vector()
  model1.E.C.prmse.2 <- model1.E.C.pcor.2 <- model1.Er.C.prmse.2 <- model1.Er.C.pcor.2 <- vector()
  model1.Ec.C.prmse.2 <- model1.Ec.C.pcor.2 <-  model1.E.RC.prmse.2 <- vector()
  model1.E.RC.pcor.2 <- model1.Er.RC.prmse.2 <- model1.Er.RC.pcor.2 <- model1.Ec.RC.prmse.2 <- model1.Ec.RC.pcor.2 <- vector()
  model1.E.fail.2 <- model1.Er.fail.2 <- model1.Ec.fail.2 <- model1.E.R.fail.2 <- model1.Er.R.fail.2 <- vector()
  model1.Ec.R.fail.2 <- model1.E.C.fail.2 <- model1.Ec.C.fail.2 <- model1.Er.C.fail.2 <- model1.E.RC.fail.2<- vector()
  model1.Er.RC.fail.2 <- model1.Ec.RC.fail.2 <- vector()
  model1.G.prmse.2 <- model1.G.pcor.2 <- model1.Gr.prmse.2 <- model1.G.fail.2 <- model1.Gr.fail.2 <- vector()
  model1.Gr.pcor.2 <- model1.Gc.prmse.2 <- model1.Gc.pcor.2 <-  model1.Gc.fail.2 <- model1.G.R.fail.2 <- vector()
  model1.G.R.prmse.2 <- model1.G.R.pcor.2 <- model1.Gr.R.prmse.2 <- model1.Gr.R.pcor.2 <- model1.Gr.R.fail.2 <- model1.Gc.R.prmse.2 <- vector()
  model1.Gc.R.pcor.2 <- model1.Gc.R.fail.2 <- model1.G.C.prmse.2 <- model1.G.C.pcor.2 <- model1.G.C.fail.2 <- vector()
  model1.Gr.C.prmse.2 <- model1.Gr.C.pcor.2 <- model1.Gr.C.fail.2 <- model1.Gc.C.prmse.2 <- model1.Gc.C.pcor.2 <- model1.Gc.C.fail.2 <- vector()
  model1.G.RC.prmse.2 <- model1.G.RC.pcor.2 <- model1.G.RC.fail.2 <-  model1.Gr.RC.prmse.2 <- model1.Gr.RC.pcor.2 <- model1.Gr.RC.fail.2 <- vector()
  model1.Gc.RC.prmse.2 <- model1.Gc.RC.pcor.2 <- model1.Gc.RC.fail.2 <- model1.S.prmse.2 <- vector()
  model1.S.pcor.2 <- model1.Sr.prmse.2 <- model1.Sr.pcor.2 <- model1.Sc.prmse.2 <- model1.Sc.pcor.2 <- vector()
  model1.S.R.prmse.2 <- model1.S.R.pcor.2 <- model1.Sr.R.prmse.2 <- vector()
  model1.Sr.R.pcor.2 <- model1.Sc.R.prmse.2 <- model1.Sc.R.pcor.2 <-  vector()
  model1.S.C.prmse.2 <- model1.S.C.pcor.2 <- model1.Sr.C.prmse.2 <- model1.Sr.C.pcor.2 <- model1.Sc.C.prmse.2 <- vector()
  model1.Sc.C.pcor.2 <-  model1.S.RC.prmse.2 <- model1.S.RC.pcor.2 <- vector()
  model1.Sr.RC.prmse.2 <- model1.Sr.RC.pcor.2 <- model1.Sc.RC.prmse.2 <- model1.Sc.RC.pcor.2 <-  vector()
  model1.S.fail.2 <- model1.Sr.fail.2 <- model1.Sc.fail.2 <- model1.S.R.fail.2 <- model1.Sr.R.fail.2 <- vector()
  model1.Sc.R.fail.2 <- model1.S.C.fail.2 <- model1.Sr.C.fail.2 <- model1.Sc.C.fail.2 <- vector()
  model1.S.RC.fail.2 <- model1.Sr.RC.fail.2<- model1.Sc.RC.fail.2 <- vector()
  model2.AR.prmse.2 <- model2.AR.pcor.2 <- model2.ARr.prmse.2 <- model2.ARr.prmse.2 <- vector()
  model2.ARc.prmse.2 <- model2.ARc.pcor.2 <-  model2.AR.R.prmse.2 <- vector()
  model2.AR.R.pcor.2 <- model2.AR.C.prmse.2 <- model2.AR.C.pcor.2 <- model2.AR.RC.prmse.2 <- model2.AR.RC.pcor.2 <- vector()
  model2.ARr.R.prmse.2 <- model2.ARr.R.pcor.2 <- model2.ARr.C.prmse.2 <- model2.ARr.C.pcor.2 <- model2.ARr.RC.prmse.2 <- vector()
  model2.ARr.RC.pcor.2 <- model2.ARc.R.prmse.2 <- model2.ARc.R.pcor.2 <- model2.ARc.C.prmse.2 <- model2.ARc.C.pcor.2 <- vector()
  model2.ARc.RC.prmse.2 <- model2.ARc.RC.pcor.2 <-  model2.AR.fail.2 <- model2.ARr.fail.2 <- vector()
  model2.ARc.fail.2 <- model2.AR.R.fail.2 <- model2.ARr.R.fail.2 <- model2.ARc.R.fail.2 <- vector()
  model2.AR.C.fail.2 <- model2.ARr.C.fail.2 <- model2.ARc.C.fail.2 <- model2.AR.RC.fail.2 <- vector()
  model2.ARr.RC.fail.2 <- model2.ARc.RC.fail.2 <- model.phi <- vector()  


#k-fold cross-validation  
  #splitting data - same genotypes not present in training and test data
  data.geno <- unique(data2[,genotype])

  for(m in c(1:mkfold)){
  idx <- sample(rep(1:kfold, length.out=length(data.geno)))
  test.geno[[m]] <- split(data.geno, idx) 
  } 


# CV for base model
for(m in c(1:length(test.geno))){
for (s in c(1:length(test.geno[[m]]))){ 
  test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
  train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]

  y <- train.data[,trait[i]]

  #X design matrix 
  X.train <- model.matrix(fixed.effect,data=train.data)
  X.test <- model.matrix(fixed.effect,data=test.data)

  #Identity matrix for training data
  Identity <- diag(nrow(train.data))

  idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
    Z.geno <- model.matrix(~idg - 1)

  idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
    Ztest.geno <- model.matrix(~idgt - 1)

  G.geno <- Z.geno%*%K%*%t(Z.geno) 

#base model
  base <- try(regress(y ~ X.train, ~G.geno, pos= rep(TRUE,2), tol = 1e-4, data = train.data),silent = TRUE)

  if(class(base) != "try-error"){

  R.residual <- Identity * base$sigma[[2]]

  Khat <- K * base$sigma[[1]] 
  Vhat <- G.geno * base$sigma[[1]] + R.residual

  gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - base$fitted)

  gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% base$beta)
    yhat.test <- gamma.test[,1] + gamma.test[,2]


base.prmse[m,s] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
base.pcor[m,s] <- cor(yhat.test, test.data[,trait[i]])
}else{
  base.fail[m,s] <- 1
}#end of try 
} #end of each CV
} #end of multiple CVs

#for comparison and plotting
base.prmse.2 <- mean(base.prmse, na.rm = TRUE)
base.pcor.2 <- mean(base.pcor, na.rm = TRUE)
base.fail.2 <- sum(base.fail, na.rm = TRUE)


    
#model 1 - add correlation structures to selected base model
  #updating base model to model1   
  #Exponential, Gaussian, Spherical - isotropic, row, column, row + column directional 
  for(m in c(1:length(test.geno))){ 
    model1.E.prmse <- model1.E.pcor <- model1.E.fail <- matrix(NA,kfold,length(phi))
    model1.G.prmse <- model1.G.pcor <- model1.G.fail <- matrix(NA,kfold,length(phi))
    model1.S.prmse <- model1.S.pcor <- model1.S.fail <- matrix(NA,kfold,length(phi))
    for (p in c(1:length(phi))){
      spR.E <- exp(-(distMat/phi[p])) 

      spR.G <- exp(-(distMat/phi[p])^2)

      distMat.s <- distMat
      distMat.s[distMat.s > phi[p]] <- NA     
      spR.S <- 1-(1.5*(distMat.s/phi[p])) + 0.5*((distMat.s/phi[p])^3)
      spR.S[is.na(spR.S)] <- 0 

      for (s in c(1:length(test.geno[[m]]))){ 
        test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
        train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]
        
        y <- train.data[,trait[i]]

        #X design matrix 
        X.train <- model.matrix(fixed.effect,data=train.data)
        X.test <- model.matrix(fixed.effect,data=test.data)

        #Identity matrix for training data
        Identity <- diag(nrow(train.data))

        #genotypic design matrix
        idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
          Z.geno <- model.matrix(~idg - 1)

        idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
          Ztest.geno <- model.matrix(~idgt - 1)

        G.geno <- Z.geno%*%K%*%t(Z.geno) 

        #distance design matrix
        idDist <- factor(as.character(train.data$Slno),levels =rownames(distMat))
        Z.spDist <- model.matrix(~idDist-1)

        idDist.test <- factor(as.character(test.data$Slno),levels =rownames(distMat))
        Ztest.spDist <- model.matrix(~idDist.test-1)        


        if (all(is.nan(diag(spR.E)))){ 
        model1.E <- "Pure.Nugget"
        model1.E.prmse [s,p] <- "PN"
        model1.E.pcor [s,p] <- "PN"
        } else{
          G.spR.E <- Z.spDist%*%spR.E%*%t(Z.spDist)
          #iso - exponential
          model1.E <- try(regress(y ~ X.train, ~ G.geno + G.spR.E, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.E) != "try-error"){

          R.residual <- Identity * model1.E$sigma[[3]]

          Khat <- K * model1.E$sigma[[1]] 
          Shat <- spR.E * model1.E$sigma[[2]]
          Vhat <- G.geno * model1.E$sigma[[1]] + G.spR.E * model1.E$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.E$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.E$fitted)

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.E$beta)  
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.E.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.E.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.E.fail[s,p] <- 1
          }  # end of try 
        } # end of else 
print("iso-Exp")

        if (all(is.nan(diag(spR.G)))){ 
        model1.G <- "Pure.Nugget"
        model1.G.prmse [s,p] <- "PN"
        model1.G.pcor [s,p] <- "PN"
        } else{
          G.spR.G <- Z.spDist%*%spR.G%*%t(Z.spDist)
          model1.G <- try(regress(y ~ X.train, ~ G.geno + G.spR.G, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.G) != "try-error"){

          R.residual <- Identity * model1.G$sigma[[3]]

          Khat <- K * model1.G$sigma[[1]] 
          Shat <- spR.G * model1.G$sigma[[2]]
          Vhat <- G.geno * model1.G$sigma[[1]] + G.spR.G * model1.G$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.G$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.G$fitted)
          
          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.G$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
          yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.G.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.G.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.G.fail[s,p] <- 1
          }         # end of try     
        } # end of else 
print("iso-Gaus")

        if (all(is.nan(diag(spR.S)))){ 
        model1.S <- "Pure.Nugget"
        model1.S.prmse [s,p] <- "PN"
        model1.S.pcor [s,p] <- "PN"
        } else{
          G.spR.S <- Z.spDist%*%spR.S%*%t(Z.spDist)
          model1.S <- try(regress(y ~ X.train, ~ G.geno + G.spR.S, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.S) != "try-error"){

          R.residual <- Identity * model1.S$sigma[[3]]

          Khat <- K * model1.S$sigma[[1]] 
          Shat <- spR.S * model1.S$sigma[[2]]
          Vhat <- G.geno * model1.S$sigma[[1]] + G.spR.S * model1.S$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.S$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.S$fitted)        

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.S$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.S.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.S.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } # end of try       
        } # end of else 
print("iso-sph")
        } # end of each CV 
      } # end of phi
      model1.E.prmse.1[m,] <-  apply(model1.E.prmse,2,mean, na.rm = TRUE) 
      model1.E.pcor.1[m,] <- apply(model1.E.pcor,2,mean, na.rm = TRUE) 
      model1.E.fail.1[m,] <- apply(model1.E.fail,2,sum, na.rm = TRUE)       

      model1.G.prmse.1[m,] <- apply(model1.G.prmse,2,mean, na.rm = TRUE)
      model1.G.pcor.1[m,] <- apply(model1.G.pcor,2,mean, na.rm = TRUE)
      model1.G.fail.1[m,] <- apply(model1.G.fail,2,sum, na.rm = TRUE)      

      model1.S.prmse.1[m,] <- apply(model1.S.prmse,2,mean, na.rm = TRUE)
      model1.S.pcor.1[m,] <- apply(model1.S.pcor,2,mean, na.rm = TRUE) 
      model1.S.fail.1[m,] <- apply(model1.S.fail, 2, sum, na.rm = TRUE)                  
    } # end of multiple CVs

    model1.E.prmse.2 <- apply(model1.E.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.E.1 <- which(model1.E.prmse.2 == min(model1.E.prmse.2))
    phi.E <- phi[phi.E.1]
    model1.E.prmse.3 <- model1.E.prmse.2[phi.E.1] # for comparing various models

    model1.E.pcor.2 <- apply(model1.E.pcor.1,2,mean) #for plotting CV pcor
    model1.E.pcor.3 <- model1.E.pcor.2[which(model1.E.pcor.2 == max(model1.E.pcor.2))] #for comparing various models

    model1.E.fail.2 <- apply(model1.E.fail.1,2,sum)    

    model1.G.prmse.2 <- apply(model1.G.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.G.1 <- which(model1.G.prmse.2 == min(model1.G.prmse.2))
    phi.G <- phi[phi.G.1]
    model1.G.prmse.3 <- model1.G.prmse.2[phi.G.1] # for comparing various models

    model1.G.pcor.2 <- apply(model1.G.pcor.1,2,mean) #for plotting CV pcor
    model1.G.pcor.3 <- model1.G.pcor.2[which(model1.G.pcor.2 == max(model1.G.pcor.2))] #for comparing various models

    model1.G.fail.2 <- apply(model1.G.fail.1,2,sum)    

    model1.S.prmse.2 <- apply(model1.S.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.S.1 <- which(model1.S.prmse.2 == min(model1.S.prmse.2))
    phi.S <- phi[phi.S.1]
    model1.S.prmse.3 <- model1.S.prmse.2[phi.S.1] # for comparing various models

    model1.S.pcor.2 <- apply(model1.S.pcor.1,2,mean) #for plotting CV pcor
    model1.S.pcor.3 <- model1.S.pcor.2[which(model1.S.pcor.2 == max(model1.S.pcor.2))] #for comparing various models

    model1.S.fail.2 <- apply(model1.S.fail.1,2,sum)
    
print("end of iso")  

  #Exponential,Gaussian,Spherical - range direction  
  for(m in c(1:length(test.geno))){ 
    model1.Er.prmse <- model1.Er.pcor <- model1.Er.fail <- matrix(NA,kfold,length(phi))    
    model1.Gr.prmse <- model1.Gr.pcor <- model1.Gr.fail <- matrix(NA,kfold,length(phi))
    model1.Sr.prmse <- model1.Sr.pcor <- model1.Sr.fail <- matrix(NA,kfold,length(phi))
   
    for (p in c(1:length(phi))){
      spR.E <-  exp(-(rowDistMat2/phi[p]))  

      spR.G <-  exp(-(rowDistMat2/phi[p])^2)

      rowDistMat2.s <- rowDistMat2
      rowDistMat2.s[rowDistMat2.s > phi[p]] <- NA
      spR.S <- 1-(1.5*(rowDistMat2.s/phi[p])) + 0.5*((rowDistMat2.s/phi[p])^3)
      spR.S[is.na(spR.S)] <- 0      

      for (s in c(1:length(test.geno[[m]]))){ 
        test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
        train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]
        
        y <- train.data[,trait[i]]

        #X design matrix 
        X.train <- model.matrix(fixed.effect, data=train.data)
        X.test <- model.matrix(fixed.effect, data = test.data)
        
        #Identity matrix for training data
        Identity <- diag(nrow(train.data))

        #genotypic design matrix
        idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
          Z.geno <- model.matrix(~idg - 1)

        idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
          Ztest.geno <- model.matrix(~idgt - 1)

        G.geno <- Z.geno%*%K%*%t(Z.geno) 

        #distance design matrix
        idDist <- factor(as.character(train.data$Slno),levels =rownames(distMat))
        Z.spDist <- model.matrix(~idDist-1)

        idDist.test <- factor(as.character(test.data$Slno),levels =rownames(distMat))
        Ztest.spDist <- model.matrix(~idDist.test-1)
        
        if (all(is.nan(diag(spR.E)))){ 
        model1.Er <- "Pure.Nugget"
        model1.Er.prmse [s,p] <- "PN"
        model1.Er.pcor [s,p] <- "PN"
        } else{
          G.spR.E <- Z.spDist%*%spR.E%*%t(Z.spDist)
          model1.Er <- try(regress(y ~ X.train, ~ G.geno + G.spR.E, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Er) != "try-error"){

          R.residual <- Identity * model1.Er$sigma[[3]]

          Khat <- K * model1.Er$sigma[[1]] 
          Shat <- spR.E * model1.Er$sigma[[2]]
          Vhat <- G.geno * model1.Er$sigma[[1]] + G.spR.E * model1.Er$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Er$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Er$fitted)         

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Er$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Er.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Er.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.Er.fail[s,p] <- 1
          }        # end of try 
        } # end of else     

print("range -exp")
      if (all(is.nan(diag(spR.G)))){ 
        model1.Gr <- "Pure.Nugget"
        model1.Gr.prmse [s,p] <- "PN"
        model1.Gr.pcor [s,p] <- "PN"
        } else{
          G.spR.G <- Z.spDist%*%spR.G%*%t(Z.spDist)
          model1.Gr <- try(regress(y ~ X.train, ~ G.geno + G.spR.G, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Gr) != "try-error"){

          R.residual <- Identity * model1.Gr$sigma[[3]]

          Khat <- K * model1.Gr$sigma[[1]] 
          Shat <- spR.G * model1.Gr$sigma[[2]]
          Vhat <- G.geno * model1.Gr$sigma[[1]] + G.spR.G * model1.Gr$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Gr$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Gr$fitted)
        
          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Gr$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Gr.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Gr.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.Gr.fail[s,p] <- 1
          }         # end of try       
        } # end of else  
print("range - Gaus")

        if (all(is.nan(diag(spR.S)))){ 
        model1.Sr <- "Pure.Nugget"
        model1.Sr.prmse [s,p] <- "PN"
        model1.Sr.pcor [s,p] <- "PN"
        } else{
          G.spR.S <- Z.spDist%*%spR.S%*%t(Z.spDist)
          model1.Sr <- try(regress(y ~ X.train, ~ G.geno + G.spR.S, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Sr) != "try-error"){

          R.residual <- Identity * model1.Sr$sigma[[3]]

          Khat <- K * model1.Sr$sigma[[1]] 
          Shat <- spR.S * model1.Sr$sigma[[2]]
          Vhat <- G.geno * model1.Sr$sigma[[1]] + G.spR.S * model1.Sr$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Sr$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Sr$fitted)          

          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Sr$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Sr.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Sr.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.Sr.fail[s,p] <- 1
          }  # end of try 
        } # end of else 
        } # end of each CV 
      } # end of phi
print("range - Sph")
      model1.Er.prmse.1[m,] <-  apply(model1.Er.prmse,2,mean, na.rm = TRUE) 
      model1.Er.pcor.1[m,] <- apply(model1.Er.pcor,2,mean, na.rm = TRUE) 
      model1.Er.fail.1[m,] <- apply(model1.Er.fail,2,sum, na.rm = TRUE)
      
      model1.Gr.prmse.1[m,] <- apply(model1.Gr.prmse,2,mean,na.rm = TRUE)
      model1.Gr.pcor.1[m,] <- apply(model1.Gr.pcor,2,mean,na.rm = TRUE)
      model1.Gr.fail.1[m,] <- apply(model1.Gr.fail,2,sum,na.rm = TRUE)  
     
      model1.Sr.prmse.1[m,] <- apply(model1.Sr.prmse,2,mean,na.rm = TRUE)
      model1.Sr.pcor.1[m,] <- apply(model1.Sr.pcor,2,mean,na.rm = TRUE) 
      model1.Sr.fail.1[m,] <- apply(model1.Sr.fail,2,sum, na.rm = TRUE)                  
    } # end of multiple CVs
print("range - .1")

    model1.Er.prmse.2 <- apply(model1.Er.prmse.1,2,mean) # for plotting CV prmse
    phi.Er.1 <- which(model1.Er.prmse.2 == min(model1.Er.prmse.2))
    phi.Er <- phi[phi.Er.1]
    model1.Er.prmse.3 <- model1.Er.prmse.2[phi.Er.1] # for comparing various models

    model1.Er.pcor.2 <- apply(model1.Er.pcor.1,2,mean) #for plotting CV pcor
    model1.Er.pcor.3 <- model1.Er.pcor.2[which(model1.Er.pcor.2 == max(model1.Er.pcor.2))] #for comparing various models  

    model1.Er.fail.2 <- apply(model1.Er.fail.1,2,sum)
    
    model1.Gr.prmse.2 <- apply(model1.Gr.prmse.1,2,mean) # for plotting CV prmse
    phi.Gr.1 <- which(model1.Gr.prmse.2 == min(model1.Gr.prmse.2))
    phi.Gr <- phi[phi.Gr.1]
    model1.Gr.prmse.3 <- model1.Gr.prmse.2[phi.Gr.1] # for comparing various models

    model1.Gr.pcor.2 <- apply(model1.Gr.pcor.1,2,mean) #for plotting CV pcor
    model1.Gr.pcor.3 <- model1.Gr.pcor.2[which(model1.Gr.pcor.2 == max(model1.Gr.pcor.2))] #for comparing various models 

    model1.Gr.fail.2 <- apply(model1.Gr.fail.1,2,sum)     

    model1.Sr.prmse.2 <- apply(model1.Sr.prmse.1,2,mean) # for plotting CV prmse
    phi.Sr.1 <- which(model1.Sr.prmse.2 == min(model1.Sr.prmse.2))
    phi.Sr <- phi[phi.Sr.1]
    model1.Sr.prmse.3 <- model1.Sr.prmse.2[phi.Sr.1] # for comparing various models

    model1.Sr.pcor.2 <- apply(model1.Sr.pcor.1,2,mean) #for plotting CV pcor
    model1.Sr.pcor.3 <- model1.Sr.pcor.2[which(model1.Sr.pcor.2 == max(model1.Sr.pcor.2))] #for comparing various models 

    model1.Sr.fail.2 <- apply(model1.Sr.fail.1,2,sum)
    
print("end of range direction")

  #Exponential, Gaussian, Spherical - column direction  
  for(m in c(1:length(test.geno))){ 
    model1.Ec.prmse <- model1.Ec.pcor <- model1.Ec.fail <- matrix(NA,kfold,length(phi))
    model1.Gc.prmse <- model1.Gc.pcor <- model1.Gc.fail <- matrix(NA,kfold,length(phi))
    model1.Sc.prmse <- model1.Sc.pcor <- model1.Sc.fail <- matrix(NA,kfold,length(phi))
    for (p in c(1:length(phi))){
      spR.E <-  exp(-(colDistMat2/phi[p])) 

      spR.G <-  exp(-(colDistMat2/phi[p])^2)

      colDistMat2.s <- colDistMat2
      colDistMat2.s[colDistMat2.s > phi[p]] <- NA
      spR.S <- 1-(1.5*(colDistMat2.s/phi[p])) + 0.5*((colDistMat2.s/phi[p])^3)
      spR.S[is.na(spR.S)] <- 0

      for (s in c(1:length(test.geno[[m]]))){ 
        test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
        train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]
        
        y <- train.data[,trait[i]]

        #X design matrix 
        X.train <- model.matrix(fixed.effect,data=train.data)
        X.test <- model.matrix(fixed.effect, data=test.data)
        
        #Identity matrix for training data
        Identity <- diag(nrow(train.data))

        #genotypic design matrix
        idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
          Z.geno <- model.matrix(~idg - 1)

        idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
          Ztest.geno <- model.matrix(~idgt - 1)

        G.geno <- Z.geno%*%K%*%t(Z.geno) 

        #distance design matrix
        idDist <- factor(as.character(train.data$Slno),levels =rownames(distMat))
        Z.spDist <- model.matrix(~idDist-1)

        idDist.test <- factor(as.character(test.data$Slno),levels =rownames(distMat))
        Ztest.spDist <- model.matrix(~idDist.test-1)     


        if (all(is.nan(diag(spR.E)))){ 
        model1.Ec <- "Pure.Nugget"
        model1.Ec.prmse [s,p] <- "PN"
        model1.Ec.pcor [s,p] <- "PN"
        } else{
          G.spR.E <- Z.spDist%*%spR.E%*%t(Z.spDist)
          model1.Ec <- try(regress(y ~ X.train, ~ G.geno + G.spR.E, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Ec) != "try-error"){

          R.residual <- Identity * model1.Ec$sigma[[3]]

          Khat <- K * model1.Ec$sigma[[1]] 
          Shat <- spR.E * model1.Ec$sigma[[2]]
          Vhat <- G.geno * model1.Ec$sigma[[1]] + G.spR.E * model1.Ec$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Ec$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Ec$fitted)
         
          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Ec$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Ec.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Ec.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          }else{
            model1.Ec.fail[s,p] <- 1
          }       # end of try     
        } # end of else  
print("column - exp")   

        if (all(is.nan(diag(spR.G)))){ 
        model1.Gc <- "Pure.Nugget"
        model1.Gc.prmse [s,p] <- "PN"
        model1.Gc.pcor [s,p] <- "PN"
        } else{
          G.spR.G <- Z.spDist%*%spR.G%*%t(Z.spDist)
          model1.Gc <- try(regress(y ~ X.train, ~ G.geno + G.spR.G, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Gc) != "try-error"){

          R.residual <- Identity * model1.Gc$sigma[[3]]

          Khat <- K * model1.Gc$sigma[[1]] 
          Shat <- spR.G * model1.Gc$sigma[[2]]
          Vhat <- G.geno * model1.Gc$sigma[[1]] + G.spR.G * model1.Gc$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Gc$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Gc$fitted)
         
          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Gc$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Gc.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Gc.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.Gc.fail[s,p] <- 1
          }      # end of try 
        } # end of else 
print("column - Gaus")

        if (all(is.nan(diag(spR.S)))){ 
        model1.Sc <- "Pure.Nugget"
        model1.Sc.prmse [s,p] <- "PN"
        model1.Sc.pcor [s,p] <- "PN"
        } else{
          G.spR.S <- Z.spDist%*%spR.S%*%t(Z.spDist)
          model1.Sc <- try(regress(y ~ X.train, ~ G.geno + G.spR.S, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

          if(class(model1.Sc) != "try-error"){

          R.residual <- Identity * model1.Sc$sigma[[3]]

          Khat <- K * model1.Sc$sigma[[1]] 
          Shat <- spR.S * model1.Sc$sigma[[2]]
          Vhat <- G.geno * model1.Sc$sigma[[1]] + G.spR.S * model1.Sc$sigma[[2]] + R.residual

          gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model1.Sc$fitted)
          delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model1.Sc$fitted)
          
          gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model1.Sc$beta)
          delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
            yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

          model1.Sc.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
          model1.Sc.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
          } else{
            model1.Sc.fail[s,p] <- 1
          }          # end of try 
        } # end of else 
print("column - Sph")
        } # end of each CV 
      } # end of phi
      model1.Ec.prmse.1[m,] <-  apply(model1.Ec.prmse,2,mean, na.rm = TRUE) 
      model1.Ec.pcor.1[m,] <- apply(model1.Ec.pcor,2,mean,na.rm = TRUE)  
      model1.Ec.fail.1[m,] <- apply(model1.Ec.fail,2,sum,na.rm = TRUE)
 
      model1.Gc.prmse.1[m,] <-  apply(model1.Gc.prmse,2,mean,na.rm = TRUE) 
      model1.Gc.pcor.1[m,] <- apply(model1.Gc.pcor,2,mean,na.rm = TRUE)
      model1.Gc.fail.1[m,] <- apply(model1.Gc.fail,2,sum, na.rm = TRUE) 

      model1.Sc.prmse.1[m,] <-  apply(model1.Sc.prmse,2,mean,na.rm = TRUE) 
      model1.Sc.pcor.1[m,] <- apply(model1.Sc.pcor,2,mean,na.rm = TRUE) 
      model1.Sc.fail.1[m,] <- apply(model1.Sc.fail,2,sum, na.rm = TRUE) 
    } # end of multiple CVs
print("column - CVs")

    model1.Ec.prmse.2 <- apply(model1.Ec.prmse.1,2,mean) # for plotting CV prmse
    phi.Ec.1 <- which(model1.Ec.prmse.2 == min(model1.Ec.prmse.2))
    phi.Ec <- phi[phi.Ec.1]
    model1.Ec.prmse.3 <- model1.Ec.prmse.2[phi.Ec.1] # for comparing various models

    model1.Ec.pcor.2 <- apply(model1.Ec.pcor.1,2,mean) #for plotting CV pcor
    model1.Ec.pcor.3 <- model1.Ec.pcor.2[which(model1.Ec.pcor.2 == max(model1.Ec.pcor.2))] #for comparing various models 

    model1.Ec.fail.2 <- apply(model1.Ec.fail.1,2,sum)

    model1.Gc.prmse.2 <- apply(model1.Gc.prmse.1,2,mean) # for plotting CV prmse
    phi.Gc.1 <- which(model1.Gc.prmse.2 == min(model1.Gc.prmse.2))
    phi.Gc <- phi[phi.Gc.1]
    model1.Gc.prmse.3 <- model1.Gc.prmse.2[phi.Gc.1] # for comparing various models

    model1.Gc.pcor.2 <- apply(model1.Gc.pcor.1,2,mean) #for plotting CV pcor
    model1.Gc.pcor.3 <- model1.Gc.pcor.2[which(model1.Gc.pcor.2 == max(model1.Gc.pcor.2))] #for comparing various models 

    model1.Gc.fail.2 <- apply(model1.Gc.fail.1,2,sum)

    model1.Sc.prmse.2 <- apply(model1.Sc.prmse.1,2,mean) # for plotting CV prmse
    phi.Sc.1 <- which(model1.Sc.prmse.2 == min(model1.Sc.prmse.2))
    phi.Sc <- phi[phi.Sc.1]
    model1.Sc.prmse.3 <- model1.Sc.prmse.2[phi.Sc.1] # for comparing various models

    model1.Sc.pcor.2 <- apply(model1.Sc.pcor.1,2,mean) #for plotting CV pcor
    model1.Sc.pcor.3 <- model1.Sc.pcor.2[which(model1.Sc.pcor.2 == max(model1.Sc.pcor.2))] #for comparing various models 

    model1.Sc.fail.2 <- apply(model1.Sc.fail.1,2,sum)
print("end of column dir")  

#model 2 - AR1 models 
  #updating base model to model2 
  #AR1 - isotropic, row, colum, and row + colum direction 
  for(m in c(1:length(test.geno))){ 
    model2.AR.prmse <- model2.AR.pcor <- model2.AR.fail <- matrix(NA,kfold,length(phi.2))
    model2.ARr.prmse <- model2.ARr.pcor <- model2.ARr.fail <- matrix(NA,kfold,length(phi.2))
    model2.ARc.prmse <- model2.ARc.pcor <- model2.ARc.fail <- matrix(NA,kfold,length(phi.2))
 
    for (p in c(1:length(phi.2))){
      rcAR <-  phi.2[p]^distMat
      rAR <-  phi.2[p]^rowDistMat2
      cAR <- phi.2[p]^colDistMat2      
      for (s in c(1:length(test.geno[[m]]))){ 
        test.data <- data2[(data2[,genotype] %in% test.geno[[m]][[s]]),]
        train.data <-data2[!(data2[,genotype] %in% test.geno[[m]][[s]]),]
        
        y <- train.data[,trait[i]]

        #X design matrix 
        X.train <- model.matrix(fixed.effect, data = train.data)
        X.test <- model.matrix(fixed.effect, data = test.data)
        
        #Identity matrix for training data
        Identity <- diag(nrow(train.data))

        #genotypic design matrix
        idg <- factor(as.character(train.data[,genotype]), levels = rownames(K))
          Z.geno <- model.matrix(~idg - 1)

        idgt <- factor(as.character(test.data[,genotype]), levels = rownames(K))
          Ztest.geno <- model.matrix(~idgt - 1)

        G.geno <- Z.geno%*%K%*%t(Z.geno) 

        #distance design matrix
        idDist <- factor(as.character(train.data$Slno),levels =rownames(distMat))
        Z.spDist <- model.matrix(~idDist-1)

        idDist.test <- factor(as.character(test.data$Slno),levels =rownames(distMat))
        Ztest.spDist <- model.matrix(~idDist.test-1)

 
        G.rcAR <- Z.spDist%*%rcAR%*%t(Z.spDist)
        model2.AR <- try(regress(y ~ X.train, ~ G.geno + G.rcAR, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

        if(class(model2.AR) != "try-error"){

        R.residual <- Identity * model2.AR$sigma[[3]]

        Khat <- K * model2.AR$sigma[[1]] 
        Shat <- rcAR * model2.AR$sigma[[2]]
        Vhat <- G.geno * model2.AR$sigma[[1]] + G.rcAR * model2.AR$sigma[[2]] + R.residual

        gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model2.AR$fitted)
        delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model2.AR$fitted)
        
        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model2.AR$beta)
        delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
          yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

        model2.AR.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model2.AR.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
        } else{
          model2.AR.fail[s,p] <- 1
        }        # end of try  
print("iso - AR")
                                                

        G.rAR <- Z.spDist%*%rAR%*%t(Z.spDist)
        model2.ARr <- try(regress(y ~ X.train, ~ G.geno + G.rAR, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

        if(class(model2.ARr) != "try-error"){

        R.residual <- Identity * model2.ARr$sigma[[3]]

        Khat <- K * model2.ARr$sigma[[1]] 
        Shat <- rAR * model2.ARr$sigma[[2]]
        Vhat <- G.geno * model2.ARr$sigma[[1]] + G.rAR * model2.ARr$sigma[[2]] + R.residual

        gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model2.ARr$fitted)
        delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model2.ARr$fitted)      

        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model2.ARr$beta)
        delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
          yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

        model2.ARr.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model2.ARr.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
        } else{
          model2.ARr.fail[s,p] <- 1
        }        # end of try    
print("range - AR")

        G.cAR <- Z.spDist%*%cAR%*%t(Z.spDist)
        model2.ARc <- try(regress(y ~ X.train, ~ G.geno + G.cAR, pos= rep(TRUE,3), tol = 1e-4,data = train.data),silent = TRUE)

        if(class(model2.ARc) != "try-error"){

        R.residual <- Identity * model2.ARc$sigma[[3]]

        Khat <- K * model2.ARc$sigma[[1]] 
        Shat <- cAR * model2.ARc$sigma[[2]]
        Vhat <- G.geno * model2.ARc$sigma[[1]] + G.cAR * model2.ARc$sigma[[2]] + R.residual

        gamma<- Khat %*% t(Z.geno) %*% solve(Vhat) %*% (y - model2.ARc$fitted)
        delta <- Shat %*% t(Z.spDist) %*% solve(Vhat) %*% (y - model2.ARc$fitted)
        
        gamma.test <- cbind(Ztest.geno %*% gamma, X.test %*% model2.ARc$beta)
        delta.test <- cbind(Ztest.spDist %*% delta, gamma.test)
          yhat.test <- delta.test[,1] + delta.test[,2] + delta.test[,3]

        model2.ARc.prmse[s,p] <- sqrt(mean((yhat.test - test.data[,trait[i]])^2))
        model2.ARc.pcor[s,p] <- cor(yhat.test, test.data[,trait[i]])  
        } else{
          model2.ARc.fail[s,p] <- 1
        }         # end of try  

print("column - AR")    
        } # end of each CV 
      } # end of phi
      model2.AR.prmse.1[m,] <-  apply(model2.AR.prmse,2,mean, na.rm = TRUE) 
      model2.AR.pcor.1[m,] <- apply(model2.AR.pcor,2,mean, na.rm = TRUE) 
      model2.AR.fail.1[m,] <- apply(model2.AR.fail,2,sum, na.rm = TRUE)

      model2.ARr.prmse.1[m,] <-  apply(model2.ARr.prmse,2,mean, na.rm = TRUE) 
      model2.ARr.pcor.1[m,] <- apply(model2.ARr.pcor,2,mean, na.rm = TRUE) 
      model2.ARr.fail.1[m,] <- apply(model2.ARr.fail,2,sum, na.rm = TRUE)

      model2.ARc.prmse.1[m,] <-  apply(model2.ARc.prmse,2,mean, na.rm = TRUE) 
      model2.ARc.pcor.1[m,] <- apply(model2.ARc.pcor,2,mean, na.rm = TRUE)
      model2.ARc.fail.1[m,] <- apply(model2.ARc.fail,2,sum, na.rm = TRUE)
     } # end of multiple CVs
print("AR - CVs")

    model2.AR.prmse.2 <- apply(model2.AR.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.AR.1 <- which(model2.AR.prmse.2 == min(model2.AR.prmse.2))
    phi.AR <- phi.2[phi.AR.1]
    model2.AR.prmse.3 <- model2.AR.prmse.2[phi.AR.1] # for comparing various models

    model2.AR.pcor.2 <- apply(model2.AR.pcor.1,2,mean) #for plotting CV pcor
    model2.AR.pcor.3 <- model2.AR.pcor.2[which(model2.AR.pcor.2 == max(model2.AR.pcor.2))] #for comparing various models

    model2.AR.fail.2 <- apply(model2.AR.fail.1,2, sum)
    
    model2.ARr.prmse.2 <- apply(model2.ARr.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.ARr.1 <- which(model2.ARr.prmse.2 == min(model2.ARr.prmse.2))
    phi.ARr <- phi.2[phi.ARr.1]
    model2.ARr.prmse.3 <- model2.ARr.prmse.2[phi.ARr.1] # for comparing various models

    model2.ARr.pcor.2 <- apply(model2.ARr.pcor.1,2,mean) #for plotting CV pcor
    model2.ARr.pcor.3 <- model2.ARr.pcor.2[which(model2.ARr.pcor.2 == max(model2.ARr.pcor.2))] #for comparing various models

    model2.ARr.fail.2 <- apply(model2.ARr.fail.1,2,sum)
    
    model2.ARc.prmse.2 <- apply(model2.ARc.prmse.1,2,mean) # for plotting CV prmse - plot this as bxpt and overlap with lines connecting min including the one value from base
    phi.ARc.1 <- which(model2.ARc.prmse.2 == min(model2.ARc.prmse.2))
    phi.ARc <- phi.2[phi.ARc.1]
    model2.ARc.prmse.3 <- model2.ARc.prmse.2[phi.ARc.1] # for comparing various models

    model2.ARc.pcor.2 <- apply(model2.ARc.pcor.1,2,mean) #for plotting CV pcor
    model2.ARc.pcor.3 <- model2.ARc.pcor.2[which(model2.ARc.pcor.2 == max(model2.ARc.pcor.2))] #for comparing various models

    model2.ARc.fail.2 <- apply(model2.ARc.fail.1,2,sum)
    
print("end of AR")

#best model from each cross-validation set from all spatial correlation structures
model1.prmse.2 <- model2.prmse.2 <- model1.pcor.2 <- model2.pcor.2 <- model.prmse <- model.pcor <- model.fail <- NULL

  model1.prmse.2 <- cbind(model1.E.prmse.2, model1.Er.prmse.2, model1.Ec.prmse.2, 
                      model1.G.prmse.2, model1.Gr.prmse.2, model1.Gc.prmse.2, model1.S.prmse.2, 
                      model1.Sr.prmse.2, model1.Sc.prmse.2)                   
                      
  model2.prmse.2 <- cbind(model2.AR.prmse.2, model2.ARr.prmse.2, model2.ARc.prmse.2)
                     
  model1.pcor.2 <- cbind(model1.E.pcor.2, model1.Er.pcor.2, model1.Ec.pcor.2, model1.G.pcor.2, 
                      model1.Gr.pcor.2, model1.Gc.pcor.2, model1.S.pcor.2, model1.Sr.pcor.2, 
                      model1.Sc.pcor.2)                   
                      
  model2.pcor.2 <- cbind(model2.AR.pcor.2, model2.ARr.pcor.2, model2.ARc.pcor.2)
                     
  model.prmse <- cbind(base.prmse.2, model1.E.prmse.3, model1.Er.prmse.3, model1.Ec.prmse.3, 
                      model1.G.prmse.3, model1.Gr.prmse.3, model1.Gc.prmse.3, model1.S.prmse.3, 
                      model1.Sr.prmse.3, model1.Sc.prmse.3, model2.AR.prmse.3, model2.ARr.prmse.3, 
                      model2.ARc.prmse.3)
                     

  model.pcor <- cbind(base.pcor.2, model1.E.pcor.3, model1.Er.pcor.3, model1.Ec.pcor.3, 
                      model1.G.pcor.3, model1.Gr.pcor.3, model1.Gc.pcor.3, model1.S.pcor.3, 
                      model1.Sr.pcor.3, model1.Sc.pcor.3, model2.AR.pcor.3, model2.ARr.pcor.3, 
                      model2.ARc.pcor.3)

  model.fail <- rbind(base.fail.2, model1.E.fail.2, model1.Er.fail.2, model1.Ec.fail.2, model1.G.fail.2, 
                      model1.Gr.fail.2, model1.Gc.fail.2, model1.S.fail.2, model1.Sr.fail.2, 
                      model1.Sc.fail.2, model2.AR.fail.2, model2.ARr.fail.2, model2.ARc.fail.2)
                        
  model.phi <- cbind(0, phi.E, phi.Er, phi.Ec, phi.G, phi.Gr, phi.Gc, phi.S, phi.Sr, phi.Sc, phi.AR, phi.ARr, phi.ARc) 
  
  #patch work
  if (nrow(model.prmse) >1){
    model.prmse <- model.prmse[1,]    
  }

  if (nrow(model.pcor) >1){
    model.pcor <- model.pcor[1,]    
  }

  if (nrow(model.phi) >1){
    model.phi <- model.phi[1,]    
  }


  phi.1 <- which(model.prmse == min(model.prmse))
  phi.win <- model.phi[phi.1]


#output tables from cross-validation
 write.csv(model1.prmse.2, paste0(trait[i],"-",loc,"_K_non-AR1_models_prmse.csv")) 
 write.csv(model1.pcor.2, paste0(trait[i],"-",loc,"_K_non-AR1_models_pcor.csv")) 
 write.csv(model2.prmse.2, paste0(trait[i],"-",loc,"_K_AR1_models_prmse.csv")) 
 write.csv(model2.pcor.2, paste0(trait[i],"-",loc,"_K_AR1_models_pcor.csv")) 
 write.csv(model.fail, paste0(trait[i],"-",loc,"_K_models_final_failed.csv")) 
 write.csv(model.prmse, paste0(trait[i],"-",loc,"_K_models_final_prmse.csv")) 
 write.csv(model.pcor, paste0(trait[i],"-",loc,"_K_models_final_pcor.csv")) 
 write.csv(model.phi, paste0(trait[i],"-",loc,"_K_models_final_phi.csv")) 


#correlation and design matrices for final model

  y <- data2[,trait[i]]

  #X design matrix 
  X.mat <- model.matrix(fixed.effect,data=data2)
 
  #genotypic design matrix
  idg <- factor(as.character(data2[,genotype]), levels = rownames(K))
    Z.geno <- model.matrix(~idg - 1)
    geno <- Z.geno%*%K%*%t(Z.geno) 

  #distance design matrix
  idDist <- factor(as.character(data2$Slno),levels =rownames(distMat))
    Z.spDist <- model.matrix(~idDist-1)
  

 #spatial relationship matrices

  spR.E <- exp(-(distMat/phi.win)) 
  spR.G <- exp(-(distMat/phi.win)^2)
  distMat.s <- distMat
  distMat.s[distMat.s > phi.win] <- NA     
  spR.S <- 1-(1.5*(distMat.s/phi.win)) + 0.5*((distMat.s/phi.win)^3)
  spR.S[is.na(spR.S)] <- 0 

  spR.Er <-  exp(-(rowDistMat2/phi.win))
  spR.Gr <-  exp(-(rowDistMat2/phi.win)^2)
  rowDistMat2.sr <- rowDistMat2
  rowDistMat2.sr[rowDistMat2.sr > phi.win] <- NA
  spR.Sr <- 1-(1.5*(rowDistMat2.sr/phi.win)) + 0.5*((rowDistMat2.sr/phi.win)^3)
  spR.Sr[is.na(spR.Sr)] <- 0 

  spR.Ec <-  exp(-(colDistMat2/phi.win))
  spR.Gc <-  exp(-(colDistMat2/phi.win)^2)
  rowDistMat2.sc <- rowDistMat2
  rowDistMat2.sc[rowDistMat2.sc > phi.win] <- NA
  spR.Sc <- 1-(1.5*(rowDistMat2.sc/phi.win)) + 0.5*((rowDistMat2.sc/phi.win)^3)
  spR.Sc[is.na(spR.Sc)] <- 0 

  rcAR <-  phi.win^distMat
  rAR <-  phi.win^rowDistMat2
  cAR <- phi.win^colDistMat2  


  G.spR.E <- Z.spDist %*% spR.E %*% t(Z.spDist)
  G.spR.G <- Z.spDist%*%spR.G%*%t(Z.spDist)
  G.spR.S <- Z.spDist%*%spR.S%*%t(Z.spDist)

  G.spR.Er <- Z.spDist%*%spR.Er%*%t(Z.spDist)
  G.spR.Gr <- Z.spDist%*%spR.Gr%*%t(Z.spDist)
  G.spR.Sr <- Z.spDist%*%spR.Sr%*%t(Z.spDist)

  G.spR.Ec <- Z.spDist%*%spR.Ec%*%t(Z.spDist)
  G.spR.Gc <- Z.spDist%*%spR.Gc%*%t(Z.spDist)
  G.spR.Sc <- Z.spDist%*%spR.Sc%*%t(Z.spDist)

  G.spR.AR <- Z.spDist%*%rcAR%*%t(Z.spDist)
  G.spR.ARr <- Z.spDist%*%rAR%*%t(Z.spDist)
  G.spR.ARc <- Z.spDist%*%cAR%*%t(Z.spDist)

  

#choosing the final best model  - default is base model

model <- switch(phi.1,
     "1" <- regress(y ~ X.mat, ~geno, pos = rep(TRUE,2), tol = 1e-4, data = data2), 
     
     "2" <- regress(y ~ X.mat, ~geno + G.spR.E, pos= rep(TRUE,3), tol = 1e-4, data = data2),
     "3" <- regress(y ~ X.mat, ~geno + G.spR.Er, pos= rep(TRUE,3), tol = 1e-4,data = data2),
     "4" <- regress(y ~ X.mat, ~geno + G.spR.Ec, pos= rep(TRUE,3), tol = 1e-4,data = data2),
          
     "5" <- regress(y ~ X.mat, ~geno + G.spR.G, pos= rep(TRUE,3), tol = 1e-4, data = data2),
     "6" <- regress(y ~ X.mat, ~geno + G.spR.Gr, pos= rep(TRUE,3), tol = 1e-4,data = data2),
     "7" <- regress(y ~ X.mat, ~geno + G.spR.Gc, pos= rep(TRUE,3), tol = 1e-4,data = data2),
          
     "8" <- regress(y ~ X.mat, ~geno + G.spR.S, pos= rep(TRUE,3), tol = 1e-4, data = data2),
     "9" <- regress(y ~ X.mat, ~geno + G.spR.Sr, pos= rep(TRUE,3), tol = 1e-4,data = data2),
     "10" <- regress(y ~ X.mat, ~geno + G.spR.Sc, pos= rep(TRUE,3), tol = 1e-4,data = data2),
          
     "11" <- regress(y ~ X.mat, ~geno + G.spR.AR, pos= rep(TRUE,3), tol = 1e-4,data = data2),
     "12" <- regress(y ~ X.mat, ~geno + G.spR.ARr, pos= rep(TRUE,3), tol = 1e-4,data = data2),
     "13" <- regress(y ~ X.mat, ~geno + G.spR.ARc, pos= rep(TRUE,3), tol = 1e-4,data = data2),
          
    regress(y ~ X.mat, ~geno, pos = rep(TRUE,2), tol = 1e-4, data = data2)
    )

  
  base <- regress(y ~ X.mat, ~geno, pos = rep(TRUE,2), tol = 1e-4, data = data2)


#summary of models
  base.formula <- capture.output(base$formula)
  base.Vformula <- capture.output(base$Vformula)    
  base.out <- capture.output(summary(base))
  base.cov <- capture.output(base$sigma.cov)

  model.formula <- capture.output(model$formula)
  model.Vformula <- capture.output(model$Vformula)
  model.out <- capture.output(summary(model))  
  model.cov <- capture.output(model$sigma.cov)


#output summary of final model
 cat(base.formula, file = paste0(trait[i],"-",loc,"_K_base_summary.txt"), sep = "\n", append=TRUE)
 cat(base.Vformula, file = paste0(trait[i],"-",loc,"_K_base_summary.txt"), sep = "\n", append=TRUE )
 cat(base.out, file= paste0(trait[i],"-",loc,"_K_base_summary.txt"), sep="\n", append=TRUE)
 cat(base.cov, file= paste0(trait[i],"-",loc,"_K_base_summary.txt"), sep="\n", append=TRUE)

 cat(model.formula, file = paste0(trait[i],"-",loc,"_K_model_summary.txt"), sep = "\n", append=TRUE)
 cat(model.Vformula, file = paste0(trait[i],"-",loc,"_K_model_summary.txt"), sep = "\n", append=TRUE )
 cat(paste0("phi = ",phi.win), file = paste0(trait[i],"-",loc,"_K_model_summary.txt"), sep = "\n", append=TRUE)
 cat(model.out, file= paste0(trait[i],"-",loc,"_K_model_summary.txt"), sep="\n", append=TRUE)
 cat(model.cov, file= paste0(trait[i],"-",loc,"_K_model_summary.txt"), sep="\n", append=TRUE)



#blup from base and models
  baseBLUP <- BLUP(base)$Mean
  modelBLUP <- BLUP(model)$Mean
  

  #separating spatial and entry blups for models
   model.spBLUP <- modelBLUP[grep("spR", names(modelBLUP), fixed=TRUE)]
      data2.sp <- cbind(data2, model.spBLUP)
    
  base.EntryBLUP <- baseBLUP[grep("geno", names(baseBLUP), fixed=TRUE)]
    
  model.EntryBLUP <- modelBLUP[grep("geno", names(modelBLUP), fixed=TRUE)]
       
  data2.entry <- cbind(data2.sp,base.EntryBLUP, model.EntryBLUP) 

#output BLUP
  write.csv(data2.sp, paste0(trait[i],"-",loc,"_K_spatial_BLUP.csv"))
  write.csv(data2.entry, paste0(trait[i],"-",loc,"_K_genotype_spatial_BLUP.csv"))

} # end of for w.r.t. trait  
} # end of function
