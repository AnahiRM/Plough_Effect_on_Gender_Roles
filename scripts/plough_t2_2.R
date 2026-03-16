  # DML: Casual Machine Learning
  # References: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017     
  #             "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016 
  
  # This empirical example uses the data from Alesina, Alberto, Paola Giuliano, and Nathan Nunn. "On the origins of gender roles: Women and the plough." The Quarterly Journal of Economics 128.2 (2013): 469-530.
  #######################################################################################################################################################
  
  ###################### Loading packages ###########################
  library(foreign);
  library(quantreg);
  library(mnormt);
  library(gbm);
  library(glmnet);
  library(MASS);
  library(rpart);
  library(foreach)
  library(iterators)
  library(parallel)
  library(doParallel)
  library(sandwich);
  library(hdm);
  library(randomForest);
  library(nnet)
  library(neuralnet)
  library(matrixStats)
  library(quadprog)
  library(xtable)
  library(ivmodel)
  library(Hmisc)
  library(dummy)
  library(DoubleML)
  rm(list = ls()) 
  trim <- FALSE
  
  ###################### Loading functions and Data ##############################
  
  #setwd(".../research/plough/t4")
  #source(".../research/functions/Moment_Functions_02.R")
  #source(".../research/functions/ML_Functions_01.R")
  
#setwd("The effect of plough agriculture on gender roles A machine learning approach (replication data)/Replication files/plough")
  # Uncomment moment function below and trim change if using trimmed propensity score, comment out the other moment function
#  source("../functions/Moment_Functions_05_v01_trim_adjusted.R")
#  trim <- TRUE
  source("functions/Moment_Functions_05_v01_adjusted2.R")
  source("functions/ML_Functions_01.R")
  
  options(warn=-1)
  set.seed(1211);
  #cl <- makeCluster(4)
  #registerDoParallel(cl)

  # Method names: Boosting, Nnet, RLasso, PostRLasso, Forest, Trees, Ridge, Lasso, Elnet, Ensemble
  
  Boosting     <- list(n.minobsinnode = 1, bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=2, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
  Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
  Lasso       <- list(s = "lambda.min",intercept = TRUE)
  RLasso       <- list(intercept = TRUE)
  Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)
  Trees        <- list(reg_method="anova", clas_method="class")
  
  arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Lasso=Lasso, Nnet=Nnet, Trees=Trees)
  
  methods      <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet")
  
  split        <- 100
  nfold        <- 2
  start_time <- Sys.time()
  
  # Column 1 ----------------------------------------------------------------
  
  
  data <- read.dta("data/crosscountry_dataset.dta")
  
  #Create Dummy Variables for Continent Fixed Effect
  data["North_America"] <- ifelse(data$continent == 'North America', 1, 0)
  data["Asia"] <- ifelse(data$continent == 'Asia', 1, 0)
  data["Europe"] <- ifelse(data$continent == 'Europe', 1, 0)
  data["Africa"] <- ifelse(data$continent == 'Africa', 1, 0)
  data["South_America"] <- ifelse(data$continent == 'South America', 1, 0)
  data["Oceania"] <- ifelse(data$continent == 'Oceania', 1, 0)
  
  data_lasso <- data[ ,c("flfp2000","plow","plow_positive_crops", "plow_negative_crops", "agricultural_suitability",
                    "tropical_climate", "large_animals", "political_hierarchies", "economic_complexity", "ln_income",
                    "ln_income_squared", "terrain_slope", "soil_depth", "avg_temperature", "avg_precipitation",
                    "terrslope2", "soil2", "avg_temp2", "avg_precip2", "North_America", "Asia", "Africa", "Europe", "South_America", "Oceania")]
  data_lasso <- na.omit(data_lasso) 
  
  data <- data[ ,c("flfp2000","plow","plow_positive_crops", "plow_negative_crops", "agricultural_suitability", 
                   "tropical_climate", "large_animals", "political_hierarchies", "economic_complexity", "ln_income", 
                   "ln_income_squared", "terrain_slope", "soil_depth", "avg_temperature", "avg_precipitation", 
                    "North_America", "Asia", "Africa", "Europe", "South_America", "Oceania")] 
  data<- na.omit(data)
  
  
  # Outcome Variable
  y  <- "flfp2000"
  
  # Treatment Variable
  d  <- "plow"
  
  #Read Instruments
  z <- c("plow_positive_crops", "plow_negative_crops")
 # z <- c("plow_positive_crops")
  # Controls
  x <- paste(colnames(data[,!colnames(data)%in% c(y,d,z, "country")]), collapse=" + ") #use this for tree-based methods like forests and boosted trees
  x_lasso <- paste(colnames(data_lasso[,!colnames(data_lasso)%in% c(y,d,z, "country")]), collapse=" + ")
  xl <- paste("(", x, ")^2", sep="") #use this for rlasso etc.
  xl_lasso <- paste("(", x_lasso, ")^2", sep="")
  ############## Arguments for DoubleML function:
  
  r_1 <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet')) %dopar% { 
    set.seed(k)
    dml <- DoubleML(data=data, data_alt=data_lasso, y=y, d=d, z=z, xx=x, xx_alt=x_lasso, xL=xl, xL_alt=xl_lasso, methods=methods, DML="DML2", nfold=nfold, est="IV", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=c(0.01,0.99))
    
    #data.frame(t(dml[1,]), t(dml[2,]), t(dml[3,]), t(dml[4,]))
  }
  
#  source("../functions/Moment_Functions_05_v01.R")
#  ensemble     <- list(methods=c("Boosting","Forest", "Nnet"))                       # methods for the ensemble estimation
#  methods      <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet","Ensemble")
  
#  r_2 <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet')) %dopar% { 
#    set.seed(k)
#    dml <- DoubleML(data=data, y=y, d=d, z=z, xx=x, xL=xl, methods=methods, DML="DML2", nfold=nfold, est="IV", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=c(0.01,0.99))
    
    #data.frame(t(dml[1,]), t(dml[2,]), t(dml[3,]), t(dml[4,]))
#  }
  
  end_time <- Sys.time()
  time <- end_time - start_time
  
  data_info <- list("n"=nrow(data), "p_col" = ncol(data), "y"= y, "d"=d, "x"=x, "xl"=xl, "methods"=methods, "time"= time, "split"=split, "nfold"=nfold)
  
  save(r_1, data_info, file = file.path("outputs", "plough_t4_c1_output_1instrument_trim.RData"))
  
  #stopCluster(cl)
  
  #### Extract median ATE and se to table
 # result           <- matrix(0,2, length(methods)+1)
#  colnames(result) <- cbind(t(methods), "best")
 # rownames(result) <- cbind("Median ATE", "se")
  
  result           <- matrix(0, 2, 7)
  colnames(result) <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet", "Ensemble", "Best")
  rownames(result) <- c("Median ATE", "se")
  
  rvLasso <- vector()
  rvLasso1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['Lasso']
    rvLasso <- append(rvLasso, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['Lasso.1']
    rvLasso1 <- append(rvLasso1, row2[1,1])
  }
  result[1,1] = median(rvLasso)
  result[2,1] = median(rvLasso1)
  
  rvTrees <- vector()
  rvTrees1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['Trees']
    rvTrees <- append(rvTrees, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['Trees.1']
    rvTrees1 <- append(rvTrees1, row2[1,1])
  }
  result[1,2] = median(rvTrees)
  result[2,2] = median(rvTrees1)
  
  rvBoosting <- vector()
  rvBoosting1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['Boosting']
    rvBoosting <- append(rvBoosting, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['Boosting.1']
    rvBoosting1 <- append(rvBoosting1, row2[1,1])
  }
  result[1,3] = median(rvBoosting)
  result[2,3] = median(rvBoosting1)
  
  rvForest <- vector()
  rvForest1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['Forest']
    rvForest <- append(rvForest, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['Forest.1']
    rvForest1 <- append(rvForest1, row2[1,1])
  }
  result[1,4] = median(rvForest)
  result[2,4] = median(rvForest1)
  
  rvNNet <- vector()
  rvNNet1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['Nnet']
    rvNNet <- append(rvNNet, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['Nnet.1']
    rvNNet1 <- append(rvNNet1, row2[1,1])
  }
  result[1,5] = median(rvNNet)
  result[2,5] = median(rvNNet1)
  
#  rvEnsemble <- vector()
#  rvEnsemble1 <- vector()
#  for (i in 1:split){
#    row1 <- r_2[i,1][[1]]$plow['Ensemble']
#    rvEnsemble <- append(rvEnsemble, row1[1,1])
#    row2 <- r_2[i,1][[1]]$plow['Ensemble.1']
#    rvEnsemble1 <- append(rvEnsemble1, row2[1,1])
#  }
#  result[1,6] = median(rvEnsemble)
#  result[2,6] = median(rvEnsemble1)
  
  # Manually set Ensemble values to NA since we skipped r_2
  result[1, 6] <- NA  # Median ATE
  result[2, 6] <- NA  # Standard Error
  
  rvBest <- vector()
  rvBest1 <- vector()
  for (i in 1:split){
    row1 <- r_1[i,1][[1]]$plow['best']
    rvBest <- append(rvBest, row1[1,1])
    row2 <- r_1[i,1][[1]]$plow['best.1']
    rvBest1 <- append(rvBest1, row2[1,1])
  }
  result[1,7] = median(rvBest)
  result[2,7] = median(rvBest1)
  print(xtable(result, type="latex", digits=3), file = file.path("outputs", "PL_T2_PB.txt"))
  
  if (trim==FALSE && length(z)==2 && nfold==2) {
    print(xtable(result, type="latex", digits=3), file=file.path("outputs", "PL_T2_PB.txt"))
  }
  
  
