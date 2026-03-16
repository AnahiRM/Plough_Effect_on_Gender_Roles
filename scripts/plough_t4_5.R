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
rm(list = ls()) 

###################### Loading functions and Data ##############################

#setwd("The effect of plough agriculture on gender roles A machine learning approach (replication data)/Replication files/plough")

source("functions/Moment_Functions_01.R")
source("functions/ML_Functions_01.R")


options(warn=-1)
set.seed(1211);
cl <- makeCluster(4)
registerDoParallel(cl)

# Method names: Boosting, Nnet, RLasso, PostRLasso, Forest, Trees, Ridge, Lasso, Elnet, Ensemble

Boosting     <- list(n.minobsinnode = 1, bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=2, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
Lasso       <- list(s = "lambda.min",intercept = TRUE)
RLasso       <- list(intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")

arguments    <- list(Boosting=Boosting, Forest=Forest, RLasso=RLasso, Lasso=Lasso, Nnet=Nnet, Trees=Trees)

ensemble     <- list(methods=c("Lasso", "Boosting", "Forest" ,"Nnet"))                       # methods for the ensemble estimation
methods      <- c("Lasso", "Trees", "Boosting", "Forest", "Nnet","Ensemble") 

split        <- 100
nfold <- 2
start_time <- Sys.time()

# Column 1 ----------------------------------------------------------------


data <- read.dta("data/crosscountry_dataset.dta")

data <- data[ ,c("women_politics","plow", "agricultural_suitability", "tropical_climate", "large_animals", "political_hierarchies", "economic_complexity", "ln_income", "ln_income_squared")] 
data <- na.omit(data)

# Outcome Variable
y  <- "women_politics"

# Treatment Variable
d  <- "plow"

# Controls
x <- paste(colnames(data[,!colnames(data)%in% c(y,d, "country")]), collapse=" + ") #use this for tree-based methods like forests and boosted trees

xl <- paste("(", x, ")^2", sep="") #use this for rlasso etc.

############## Arguments for DoubleML function:


r_5 <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','neuralnet','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet')) %dopar% { 
  set.seed(k)
  dml <- DoubleML(data=data, y=y, d=d, z=NULL, xx=x, xL=xl, methods=methods, DML="DML2", nfold=nfold, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=NULL) 
  
  data.frame(t(dml[1,]), t(dml[2,]), t(dml[3,]), t(dml[4,]))
  
}

end_time <- Sys.time()
time <- end_time - start_time

data_info <- list("n"=nrow(data), "p_col" = ncol(data), "y"= y, "d"=d, "x"=x, "xl"=xl, "methods"=methods, "time"= time, "split"=split, "nfold"=nfold)

save(r_5, time, data_info, file = file.path( "outputs" ,"plough_t4_c5_output.RData"))


#stopCluster(cl)

#### Extract median ATE and se to table
result           <- matrix(0,2, length(methods)+1)
colnames(result) <- cbind(t(methods), "best")
rownames(result) <- cbind("Median ATE", "se")

for (i in 1:(length(methods)+1)) {
  result[1,i] = median(r_5[,i])
  result[2,i] = median(r_5[,i+7])
}
print(xtable(result, type="latex", digits=3), file = file.path("outputs", "PL_T1_PC.txt"))

