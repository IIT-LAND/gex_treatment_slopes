# lasso_functions.R
# Functions to do the bulk of the heavy lifting for the LASSO models
#

#==============================================================================
.preproc <- function(dataTraining, labelsTraining, dataTest, labelsTest, dataType) {
  #--------------------------------------------------------------------------
  if (dataType=="Microarray"){
    full_model = model.matrix(~ 1 + rfx_slope + batch + sex + RIN, data=labelsTraining)
    cov_column_names = c("batch2","batchWG6","sexM","RIN")
    # cov_column_names = c("batch1","batch2","batchWG6","sexM","RIN")
    # colnames(full_model) = c("rfx_slope",cov_column_names)
    colnames(full_model) = c("intercept","rfx_slope",cov_column_names)
  } else if (dataType=="RNAseq"){
    full_model = model.matrix(~ 0 + rfx_slope + batch + sex + rin, data = labelsTraining)
    cov_column_names = c("batchp2157","batchp2321","batchp2327","batchp2341","batchp2350","batchp2354","sexM","rin")
    colnames(full_model) = c("rfx_slope",cov_column_names)
  }
  
  # fit model
  fit = lmFit(dataTraining,full_model)
  
  # remove batch, sex, and RIN
  beta1 = fit$coefficients[, cov_column_names, drop = FALSE]
  beta1[is.na(beta1)] = 0
  dataTraining = dataTraining - beta1 %*% t(full_model[,cov_column_names])
  
  if (dataType=="Microarray"){
    full_model = model.matrix(~ 1 + rfx_slope + batch + sex + RIN, data=labelsTest)
    # cov_column_names = c("batch1","batch2","batchWG6","sexM","RIN")
    cov_column_names = c("batch2","batchWG6","sexM","RIN")
    # colnames(full_model) = c("rfx_slope",cov_column_names)
    colnames(full_model) = c("intercept","rfx_slope",cov_column_names)
  } else if (dataType=="RNAseq"){
    full_model = model.matrix(~ 0 + rfx_slope + batch + sex + rin, data = labelsTest)
    cov_column_names = c("batchp2157","batchp2321","batchp2327","batchp2341","batchp2350","batchp2354","sexM","rin")
    colnames(full_model) = c("rfx_slope",cov_column_names)
  }
  
  if (dim(dataTest)[2]>1){
    dataTest = dataTest - beta1 %*% t(full_model[,cov_column_names])
  } else{
    dataTest = dataTest - beta1 %*% t(t(full_model[,cov_column_names]))
  }
  #--------------------------------------------------------------------------
  res = list(dataTraining=dataTraining, dataTest=dataTest)
  return(res)
} # function .preproc


#==============================================================================
# run LASSO regression with Leave one out (LOO)
.myLasso.LOO.new <- function(inputData, inputTrait, dataType, seed2use=123, ncores=4, PLOT=TRUE, PERM=FALSE) {
  
  # load required package
  require(glmnet)
  require(doMC)
  require(ggplot2)
  
  registerDoMC(cores = ncores)
  
  # define a range for Lambda
  range.L = 10^seq(10, -2, length.out = 100) 
  
  # define the number of samples
  nsample = ncol(inputData) 
  
  # creat empty list of the same length of the number of samples
  loocRes = rep(NA, nsample) 
  predVal = rep(NA, nsample)
  realLabels = rep(NA, nsample)
  
  # realLabels = inputTrait[,"rfx_slope"]
  
  # repeat for each sample
  for (i in 1:nsample) { 
    
    # remove 1 sample
    dataTraining = inputData[,-i] 
    dataTest = matrix(inputData[,i])
    
    # remove the trait for the same sample
    classTraining = as.matrix(inputTrait[-i,"rfx_slope"]) 
    classTest = inputTrait[i,"rfx_slope"]
    labelsTraining = inputTrait[-i,]
    labelsTest = inputTrait[i,]
    
    # #--------------------------------------------------------------------------
    # full_model = model.matrix(~ 0 + rfx_slope + batch + sex + RIN, data=labelsTraining)
    # cov_column_names = c("batch1","batch2","batchWG6","sexM","RIN")
    # colnames(full_model) = c("rfx_slope",cov_column_names)
    # 
    # # fit model
    # fit = lmFit(dataTraining,full_model)
    # 
    # # remove batch, sex, and RIN
    # beta1 = fit$coefficients[, cov_column_names, drop = FALSE]
    # beta1[is.na(beta1)] = 0
    # dataTraining = dataTraining - beta1 %*% t(full_model[,cov_column_names])
    # 
    # full_model = model.matrix(~ 0 + rfx_slope + batch + sex + RIN, data=labelsTest)
    # cov_column_names = c("batch1","batch2","batchWG6","sexM","RIN")
    # colnames(full_model) = c("rfx_slope",cov_column_names)
    # 
    # dataTest = dataTest - beta1 %*% t(t(full_model[,cov_column_names]))
    # #--------------------------------------------------------------------------
    preproc_res = .preproc(dataTraining = dataTraining, labelsTraining = labelsTraining, 
                dataTest = dataTest, labelsTest = labelsTest, dataType = dataType) 
    dataTraining = preproc_res$dataTraining
    dataTest = preproc_res$dataTest
    
    if (PERM==TRUE){
      classTraining = sample(classTraining)
    }
    
    # set seed for reproducibility
    set.seed(seed2use) 
    
    # find min Lambda with a cv of 10
    min.L = cv.glmnet(x = t(dataTraining), 
                      y = classTraining, 
                      family = "gaussian", 
                      type.measure = "mse", 
                      alpha = 1, 
                      lambda = range.L,
                      parallel=TRUE)
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                       newx = matrix(dataTest,nrow=1), 
                       s = min.L$lambda.min, 
                       type = "response") 
    
    # find the mean squared error for the left out samples
    loocRes[i] = mean((classTest-traitHat)^2) 
    predVal[i] = traitHat
    realLabels[i] = classTest
  } # for (i in 1:nsample) {
  
  # find the mean squared error for the left out samples
  loocRes = mean(loocRes) 
  tmp = data.frame(predVal, inputTrait=realLabels)
  sst = sum((realLabels-mean(realLabels))^2)
  sse = sum((predVal-realLabels)^2)
  rsq = 1-(sse/sst)
  
  # make plot
  if (PLOT){
    minVal = min(tmp$inputTrait)
    maxVal = max(tmp$inputTrait)
    gg = ggplot(data = tmp, aes(x = predVal,y = inputTrait)) 
    gg = gg + geom_point(col="steelblue", size=3) +
      geom_smooth(method = "lm", col="firebrick", size=1.5) + 
      ylim(minVal-1,maxVal+1) + xlim(minVal-1, maxVal+1) + 
      ggtitle(sprintf("R squared = %f",rsq)) + xlab("Predicted") + ylab("Observed")
    print(gg)
  } # if (PLOT)
  
  res = list(mse=loocRes, data2plot = tmp)
  # return the funtion result
  # return(loocRes) 
  return(res) 
  
} # .myLasso.LOO.new function
#==============================================================================





#==============================================================================
# run LASSO regression with Leave one out (LOO) and find features in the model
.myLasso.LOO.genes.new <- function(inputData, inputTrait, dataType, seed2use=123, ncores=4) {
  
  # load required package
  require(doMC)
  require(glmnet) 
  
  registerDoMC(cores = ncores)
  
  # define a range for Lambda
  range.L = 10^seq(10, -2, length.out = 100)
  
  # define the number of samples
  nsample = ncol(inputData) 
  
  # create empty list of the same length of the number of samples
  loocRes = rep(NA, nsample) 
  gene.list = list()
  
  # repeat for each sample
  for (i in 1:nsample) { 
    
    # remove 1 sample
    dataTraining = inputData[,-i] 
    dataTest = matrix(inputData[,i])
    
    # remove the trait for the same sample
    classTraining = as.matrix(inputTrait[-i,"rfx_slope"]) 
    classTest = inputTrait[i,"rfx_slope"]
    labelsTraining = inputTrait[-i,]
    labelsTest = inputTrait[i,]
    
    # remove covariates
    preproc_res = .preproc(dataTraining = dataTraining, labelsTraining = labelsTraining, 
                              dataTest = dataTest, labelsTest = labelsTest, dataType = dataType) 
    dataTraining = preproc_res$dataTraining
    dataTest = preproc_res$dataTest
    
    # set seed for reproducibility
    set.seed(seed2use) 
    
    # find min Lambda with a cv of 10
    min.L = cv.glmnet(x = t(dataTraining), 
                      y = classTraining, 
                      family = "gaussian", 
                      type.measure = "mse", 
                      alpha = 1, 
                      lambda = range.L,
                      parallel=TRUE)
    
    # find the genes with min lambda
    genes = as.matrix(coef(min.L, s = "lambda.min")) 
    
    # remove 1 row
    genes = genes[-c(1),]
    
    # keep genes that are not zero
    tmp = genes!=0
    
    # create a matrix of the kept genes
    tmp2 = as.matrix(subset(tmp, tmp==TRUE))
    
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                              newx = matrix(dataTest,nrow=1), 
                              s = min.L$lambda.min, 
                              type = "response") 
    
    loocRes[i] = mean((classTest-traitHat)^2)
    # make a list with genes and MSE
    gene.list = cbind(gene.list, list(genes = rownames(tmp2), score = loocRes)) 
  } # for (i in 1:nsample) {
  
  #loocRes=median(loocRes) # in the 10 fold analysi it would take the median MSE from the 10 runs, but here is just the mean of the LOO
  
  # return the funtion result
  return(gene.list) 
  
} # .myLasso.LOO.genes.new function
#==============================================================================









#==============================================================================
# run LASSO regression with cross validation in cvfold # see above for comments 
.myLasso.CV.new <- function(inputData, inputTrait, cvfold = 10, dataType, 
                            seed2use=123, ncores=4, PLOT=TRUE, PERM=FALSE) {
  
  # load required package
  require(glmnet) 
  require(doMC)
  
  registerDoMC(cores = ncores)
  
  range.L = 10^seq(10, -2, length.out = 100)
  
  nsample = ncol(inputData)
  
  # set seed for reproducibility
  set.seed(seed2use)
  
  # split the samples in 10 and train only on 9
  inputRange = split(sample(1:nsample), cut(1:nsample,cvfold)) 
  
  mse_res = numeric(length=0)
  predVal = numeric(length=0)
  trueVal = numeric(length=0)
  
  for (i in 1:length(inputRange)) {
    
    dataTraining = inputData[,-inputRange[[i]]] 
    dataTest = inputData[,inputRange[[i]]]
    
    classTraining = as.matrix(inputTrait[-inputRange[[i]],"rfx_slope"]) 
    classTest = inputTrait[inputRange[[i]],"rfx_slope"]

    labelsTraining = inputTrait[-inputRange[[i]],]
    labelsTest = inputTrait[inputRange[[i]],]
    
    preproc_res = .preproc(dataTraining = dataTraining, labelsTraining = labelsTraining, 
                              dataTest = dataTest, labelsTest = labelsTest, dataType = dataType) 
    dataTraining = preproc_res$dataTraining
    dataTest = preproc_res$dataTest
    
    if (PERM==TRUE){
      classTraining = sample(classTraining)
    } # if (PERM==TRUE){
    
    if (nsample > 30) {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        parallel=TRUE)

    } else {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        nfolds = 9,
                        parallel=TRUE)

    } # if (nsample > 30) {
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                              newx = t(dataTest), 
                              s = min.L$lambda.min, 
                              type = "response") 
    
    # compute mean squared error
    mse_res = c(mse_res, mean((classTest-traitHat)^2))
    
    # save true data and predicted values over folds
    trueVal = c(trueVal, classTest)
    predVal = c(predVal, traitHat)
    
  } # for (i in 1:length(inputRange)) {
  
  # compute R^2
  tmp = data.frame(predVal,inputTrait = trueVal)
  sst = sum((trueVal-mean(trueVal))^2)
  sse = sum((predVal-trueVal)^2)
  rsq = 1-sse/sst
  
  # make plot
  if (PLOT){
    minVal = min(tmp$inputTrait)
    maxVal = max(tmp$inputTrait)
    gg = ggplot(data = tmp, aes(x = predVal,y = inputTrait)) 
    gg = gg + geom_point(col="steelblue", size=3) +
      geom_smooth(method = "lm", col="firebrick", size=1.5) + 
      ylim(minVal-1,maxVal+1) + xlim(minVal-1, maxVal+1) +  
      ggtitle(sprintf("R squared = %f",rsq)) + xlab("Predicted") + ylab("Observed")
    print(gg)
  } # if (PLOT)
  
  # compute average mse across folds
  mean_mse_res = mean(mse_res)
  
  res = list(mse=mean_mse_res, data2plot = tmp)
  
  return(res)
  
} # .myLasso.CV.new function
#==============================================================================

#==============================================================================
# run LASSO regression with 10-fold Cross validation and find the genes in the model
.myLasso.CV.genes.new <- function(inputData, inputTrait, cvfold=10, dataType, 
                                  seed2use=123, ncores=4, PERM=FALSE) { 
  
  require(glmnet)
  require(doMC)
  
  registerDoMC(cores = ncores)
  
  range.L = 10^seq(10, -2, length.out = 100)
  
  nsample = ncol(inputData)
  
  set.seed(seed2use)
  
  inputRange = split(sample(1:nsample), cut(1:nsample,cvfold))
  
  loocRes = rep(NA,cvfold)
  
  gene.list=list()
  
  for (i in 1:length(inputRange)) {
    
    dataTraining = inputData[,-inputRange[[i]]] 
    dataTest = inputData[,inputRange[[i]]]
    
    classTraining = as.matrix(inputTrait[-inputRange[[i]],"rfx_slope"]) 
    classTest = inputTrait[inputRange[[i]],"rfx_slope"]
    
    labelsTraining = inputTrait[-inputRange[[i]],]
    labelsTest = inputTrait[inputRange[[i]],]

    if (PERM==TRUE){
      classTraining = sample(classTraining)
    } # if (PERM==TRUE){
    
    # remove covariates
    preproc_res = .preproc(dataTraining = dataTraining, labelsTraining = labelsTraining, 
                           dataTest = dataTest, labelsTest = labelsTest, dataType = dataType) 
    dataTraining = preproc_res$dataTraining
    dataTest = preproc_res$dataTest
    
    if (nsample > 30) {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        parallel=TRUE)
      
    } else {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        nfolds = 9,
                        parallel=TRUE)
      
    } # if (nsample > 30) {
    
    # find the genes with min lambda
    genes = as.matrix(coef(min.L, s = "lambda.min")) 
    
    # remove 1 row
    genes = genes[-c(1),] 
    
    # keep genes that are not zero
    tmp = genes!=0
    
    # create a matrix of the kept genes
    tmp2 = as.matrix(subset(tmp, tmp==TRUE)) 
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                              newx = t(dataTest), 
                              s = min.L$lambda.min, 
                              type = "response") 
    
    loocRes = mean((classTest-traitHat)^2)
    
    # make a list with genes and MSE
    gene.list = cbind(gene.list, list(genes = rownames(tmp2), score = loocRes)) 
  } # for (i in 1:length(inputRange)) {
  
  #loocRes=median(loocRes)
  
  return(gene.list)
  
} # .myLasso.CV.genes.new function
#==============================================================================

#==============================================================================
# define the function for the shuffled permutation using multiple cores for LOO
.myRandFn.LOO.new <- function(inputData, 
                              inputTrait, 
                              dataType,
                              nperm = 10000, 
                              seed2use=123, 
                              verbose=FALSE) {
  
  randRes = rep(NA, nperm)
  
  for (i in 1:nperm) {
    set.seed(i)
    res = .myLasso.LOO.new(inputData, 
                           inputTrait, 
                           dataType=dataType,
                           seed2use=i, 
                           PLOT=FALSE, 
                           PERM=TRUE)
    randRes[i] = res$mse
    if (verbose){
      print(sprintf("Permutation %d: MSE = %f",i,randRes[i]))
    } # if (verbose)
  } # for (i in 1:repNo) {
  
  return(randRes)
  
} # .myRandFn.LOO,new function
#==============================================================================


#==============================================================================
# define the function for the shuffled permutation using multiple cores for CV
.myRandFn.CV.new <- function(inputData, 
                             inputTrait,
                             nperm = 10000,
                             cvfold = 10,
                             dataType,
                             seed2use=123,
                             verbose=FALSE){
  
  randRes = rep(NA, nperm)
  
  for (i in 1:nperm) {
    set.seed(i)
    res = .myLasso.CV.new(inputData = inputData,
                          inputTrait = inputTrait,
                          cvfold = cvfold,
                          dataType = dataType,
                          seed2use = i,
                          PLOT = FALSE,
                          PERM=TRUE)
    randRes[i] = res$mse
    if (verbose){
      print(sprintf("Permutation %d: MSE = %f",i,randRes[i]))
    } # if (verbose)
  } # for (i in 1:repNo) {
  
  return(randRes)
  
} # .myRandFn.CV function
#==============================================================================


#==============================================================================
# run LASSO regression with cross validation in cvfold USING ONLY THE CLINCIAL DATA AND NO GENE EXPRESSION 
.myLasso.CV.new.clin2 <- function(inputData, inputTrait, cvfold = 10, dataType, 
                                 seed2use=123, ncores=4, PLOT=TRUE, PERM=FALSE, clin_vars) {
  
  # load required package
  require(glmnet) 
  require(doMC)
  
  registerDoMC(cores = ncores)
  
  range.L = 10^seq(10, -2, length.out = 100)
  
  nsample = ncol(inputData)
  
  # set seed for reproducibility
  set.seed(seed2use)
  
  # split the samples in 10 and train only on 9
  inputRange = split(sample(1:nsample), cut(1:nsample,cvfold)) 
  
  mse_res = numeric(length=0)
  predVal = numeric(length=0)
  trueVal = numeric(length=0)
  
  for (i in 1:length(inputRange)) {
    
    dataTraining = inputData[,-inputRange[[i]]] 
    dataTest = inputData[,inputRange[[i]]]
    
    classTraining = as.matrix(inputTrait[-inputRange[[i]],"rfx_slope"]) 
    classTest = inputTrait[inputRange[[i]],"rfx_slope"]
    
    labelsTraining = inputTrait[-inputRange[[i]],]
    labelsTest = inputTrait[inputRange[[i]],]
    
    # preproc_res = .preproc(dataTraining = dataTraining, labelsTraining = labelsTraining, 
    #                        dataTest = dataTest, labelsTest = labelsTest, dataType = dataType) 
    # dataTraining = preproc_res$dataTraining
    # dataTest = preproc_res$dataTest
    
    if (PERM==TRUE){
      classTraining = sample(classTraining)
    } # if (PERM==TRUE){
    
    if (nsample > 30) {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        parallel=TRUE)
      
    } else {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        nfolds = 9,
                        parallel=TRUE)
      
    } # if (nsample > 30) {
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                              newx = t(dataTest), 
                              s = min.L$lambda.min, 
                              type = "response") 
    
    # compute mean squared error
    mse_res = c(mse_res, mean((classTest-traitHat)^2))
    
    # save true data and predicted values over folds
    trueVal = c(trueVal, classTest)
    predVal = c(predVal, traitHat)
    
  } # for (i in 1:length(inputRange)) {
  
  # compute R^2
  tmp = data.frame(predVal,inputTrait = trueVal)
  sst = sum((trueVal-mean(trueVal))^2)
  sse = sum((predVal-trueVal)^2)
  rsq = 1-sse/sst
  
  # make plot
  if (PLOT){
    minVal = min(tmp$inputTrait)
    maxVal = max(tmp$inputTrait)
    gg = ggplot(data = tmp, aes(x = predVal,y = inputTrait)) 
    gg = gg + geom_point(col="steelblue", size=3) +
      geom_smooth(method = "lm", col="firebrick", size=1.5) + 
      ylim(minVal-1,maxVal+1) + xlim(minVal-1, maxVal+1) +  
      ggtitle(sprintf("R squared = %f",rsq)) + xlab("Predicted") + ylab("Observed")
    print(gg)
  } # if (PLOT)
  
  # compute average mse across folds
  mean_mse_res = mean(mse_res)
  
  res = list(mse=mean_mse_res, data2plot = tmp)
  
  return(res)
  
} # .myLasso.CV.new.clin2 function
#==============================================================================



#==============================================================================
# run LASSO regression with 10-fold Cross validation and find the genes in the model
.myLasso.CV.genes.new.clin2 <- function(inputData, inputTrait, cvfold=10, dataType, 
                                       seed2use=123, ncores=4, PERM=FALSE, clin_vars) { 
  
  require(glmnet)
  require(doMC)
  
  registerDoMC(cores = ncores)
  
  range.L = 10^seq(10, -2, length.out = 100)
  
  nsample = ncol(inputData)
  
  set.seed(seed2use)
  
  inputRange = split(sample(1:nsample), cut(1:nsample,cvfold))
  
  loocRes = rep(NA,cvfold)
  
  gene.list=list()
  
  for (i in 1:length(inputRange)) {
    
    dataTraining = inputData[,-inputRange[[i]]] 
    dataTest = inputData[,inputRange[[i]]]
    
    classTraining = as.matrix(inputTrait[-inputRange[[i]],"rfx_slope"]) 
    classTest = inputTrait[inputRange[[i]],"rfx_slope"]
    
    labelsTraining = inputTrait[-inputRange[[i]],]
    labelsTest = inputTrait[inputRange[[i]],]
    
    if (PERM==TRUE){
      classTraining = sample(classTraining)
    } # if (PERM==TRUE){
    
    # # remove covariates
    # preproc_res = .preproc(dataTraining = dataTraining_gex, labelsTraining = labelsTraining, 
    #                        dataTest = dataTest_gex, labelsTest = labelsTest, dataType = dataType) 
    # dataTraining = cbind(preproc_res$dataTraining, dataTraining_clin)
    # dataTest = cbind(preproc_res$dataTest, dataTest_clin)
    
    if (nsample > 30) {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        parallel=TRUE)
      
    } else {
      
      # find min Lambda with a cv of 10
      min.L = cv.glmnet(x = t(dataTraining), 
                        y = classTraining, 
                        family = "gaussian", 
                        type.measure = "mse", 
                        alpha = 1, 
                        lambda = range.L,
                        nfolds = 9,
                        parallel=TRUE)
      
    } # if (nsample > 30) {
    
    # find the genes with min lambda
    genes = as.matrix(coef(min.L, s = "lambda.min")) 
    
    # remove 1 row
    genes = genes[-c(1),] 
    
    # keep genes that are not zero
    tmp = genes!=0
    
    # create a matrix of the kept genes
    tmp2 = as.matrix(subset(tmp, tmp==TRUE)) 
    
    # train the model with the min lambda on the training set (all samples minus one)
    lasso.mod = glmnet(t(dataTraining), 
                       classTraining, 
                       alpha = 1, 
                       lambda = min.L$lambda.min, 
                       thresh = 1e-12) 
    
    # predict the outcome based on the trained model
    traitHat = predict.glmnet(lasso.mod, 
                              newx = t(dataTest), 
                              s = min.L$lambda.min, 
                              type = "response") 
    
    loocRes = mean((classTest-traitHat)^2)
    
    # make a list with genes and MSE
    gene.list = cbind(gene.list, list(genes = rownames(tmp2), score = loocRes)) 
  } # for (i in 1:length(inputRange)) {
  
  #loocRes=median(loocRes)
  
  return(gene.list)
  
} # .myLasso.CV.genes.new.clin2 function
#==============================================================================


#==============================================================================
# define the function for the shuffled permutation using multiple cores for CV
.myRandFn.CV.new.clin2 <- function(inputData, 
                                  inputTrait,
                                  nperm = 10000,
                                  cvfold = 10,
                                  dataType,
                                  seed2use=123,
                                  verbose=FALSE, 
                                  clin_vars){
  
  randRes = rep(NA, nperm)
  
  for (i in 1:nperm) {
    set.seed(i)
    res = .myLasso.CV.new.clin2(inputData = inputData,
                               inputTrait = inputTrait,
                               cvfold = cvfold,
                               dataType = dataType,
                               seed2use = i,
                               PLOT = FALSE,
                               PERM=TRUE, 
                               clin_vars = clin_vars)
    randRes[i] = res$mse
    if (verbose){
      print(sprintf("Permutation %d: MSE = %f",i,randRes[i]))
    } # if (verbose)
  } # for (i in 1:repNo) {
  
  return(randRes)
  
} # .myRandFn.CV.new.clin2 function
#==============================================================================
