##########################################################################################################
#  Program:     Functions for estimating moments for using Machine Learning Methods.                    #
#  Modified: footnote page 1, allowing for vector D and Z for partial linear estimation method          #
#  Based on: "Double/Debiased Machine Learning of Treatment and Causal Parameters",  AER P&P 2017       #
#              "Double Machine Learning for Treatment and Causal Parameters",  Arxiv 2016               #
#  by V.Chernozhukov, D. Chetverikov, M. Demirer, E. Duflo, C. Hansen, W. Newey                         #           
#########################################################################################################


#source("../research/functions/ML_Functions_01.R")


DoubleML <- function(data, data_alt, y, d,z, xx, xx_alt, xL, xL_alt, methods, DML, nfold, est, arguments, ensemble, silent=FALSE, trim){
  K         <- nfold
  TE        <- rep(list(matrix(0,1,(length(methods)+1))), length(d))
  STE       <- rep(list(matrix(0,1,(length(methods)+1))), length(d))
  result    <-rep(list(matrix(0,2,(length(methods)+1))), length(d))
  result2   <- rep(list(matrix(0,2,(length(methods)+1))), length(d))
  MSE1      <- matrix(0,length(methods)+1,K)
  MSE2      <- rep(list(matrix(0,length(methods)+1,K)),length(d))
  MSE3      <- rep(list(matrix(0,length(methods)+1,K)), length(z))
  MSE4      <- rep(list(matrix(0,length(methods)+1,K)), length(d))
  MSE5      <- rep(list(matrix(0,length(methods)+1,K)), length(d))
  cond.comp <- matrix(list(),length(methods),K)
  dpool     <- rep(list(vector("list", (length(methods)+1))), length(d))
  ypool     <- vector("list", (length(methods)+1))
  zopool    <- rep(list(vector("list", (length(methods)+1))), length(d))
  if(is.null(z)){zpool     <- rep(list(vector("list", (length(methods)+1))), length(d))}else{zpool     <- rep(list(vector("list", (length(methods)+1))), length(z))}
  #zpool     <- rep(list(vector("list", (length(methods)+1))), length(z))###### rename rest to d 
  z1pool    <- rep(list(vector("list", (length(methods)+1))), length(d))
  z0pool    <- rep(list(vector("list", (length(methods)+1))), length(d))
  dz1pool   <- rep(list(vector("list", (length(methods)+1))), length(d))
  dz0pool   <- rep(list(vector("list", (length(methods)+1))), length(d))
  
  binary <- rep(0, length(d))
  for(i in 1:length(d)){
    binary[i] <- as.numeric(checkBinary(data[,d[[i]]]))
  }
  
  flag      <- 0
  
  if(est=="LATE"){
    
    binary.z <- as.numeric(checkBinary(data[,z])) 
    if(!(binary.z==1)){
      print("instrument is not binary")
      stop()
    } 
    
    if(sum(!(data[data[,z]==0,d]==0))==0){
      
      flag <- 1
      
    }
  }
  
  split     <- runif(nrow(data))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
  
  for(k in 1:length(methods)){   
    
    if(silent==FALSE){
      cat(methods[k],'\n')
    }
    
    if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==methods[k])){
      x=xL
      x_alt = xL_alt
    } else {
      x=xx
      x_alt = xx_alt
    }
    
    for(j in 1:K){   
      
      if(silent==FALSE){
        cat('  fold',j,'\n')
      }
      
      ii  <- cvgroup == j
      nii <- cvgroup != j
      
      if(K==1){
        
        ii  <- cvgroup == j
        nii <- cvgroup == j
        
      }
      
      datause <- as.data.frame(data[nii,])
      dataout <- as.data.frame(data[ii,]) 
      
      datause_alt <- as.data.frame(data_alt[nii,])
      dataout_alt <- as.data.frame(data_alt[ii,]) 
      
      if(est=="LATE" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, x=x, z=z, method=methods[k], plinear=3, xL=xL, binary=binary, flag=flag, arguments=arguments, ensemble=ensemble)}
        else{                       cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, z=z, method=methods[k], plinear=3, xL=xL, binary=binary, flag=flag, arguments=arguments);  }
        
        MSE1[k,j]               <- cond.comp[[k,j]]$err.yz0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yz1
        if(flag==1){ MSE3[k,j]  <- 0}
        else{  MSE3[k,j]        <- cond.comp[[k,j]]$err.dz0 }
        MSE4[k,j]               <- cond.comp[[k,j]]$err.dz1
        MSE5[k,j]               <- cond.comp[[k,j]]$err.z2
        
        drop                   <- which(cond.comp[[k,j]]$mz2_x>trim[1] & cond.comp[[k,j]]$mz2_x<trim[2])      
        mz_x                   <- cond.comp[[k,j]]$mz2_x[drop]
        my_z1x                 <- cond.comp[[k,j]]$my_z1x[drop]
        my_z0x                 <- cond.comp[[k,j]]$my_z0x[drop]
        md_z1x                 <- cond.comp[[k,j]]$md_z1x[drop]
        if(flag==1){ md_z0x    <- matrix(0,1,length(my_z0x))}
        else{  md_z0x          <- cond.comp[[k,j]]$md_z0x[drop] }
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]
        zout                   <- dataout[drop,z]
        
        
        TE[1,k]                <- LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)/K + TE[1,k];
        STE[1,k]               <- (1/(K^2))*((SE.LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x))^2) + STE[1,k];
        
        ypool[[k]]             <- c(ypool[[k]], yout)
        dpool[[k]]             <- c(dpool[[k]], dout)
        zopool[[k]]            <- c(zopool[[k]], zout)
        zpool[[k]]             <- c(zpool[[k]], mz_x)
        z1pool[[k]]            <- c(z1pool[[k]], my_z1x)
        z0pool[[k]]            <- c(z0pool[[k]], my_z0x)
        dz1pool[[k]]           <- c(dz1pool[[k]], md_z1x)
        dz0pool[[k]]           <- c(dz0pool[[k]], md_z0x)
        
        MSE1[(length(methods)+1),j] <- error(mean(datause[datause[,z]==0,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==0,y]),y])$err
        MSE2[(length(methods)+1),j] <- error(mean(datause[datause[,z]==1,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==1,y]),y])$err
        if(flag==1){ MSE3[(length(methods)+1),j]=0    }
        else{  MSE3[(length(methods)+1),j] <- error(mean(datause[datause[,z]==0,d], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==0,d]),d])$err }
        MSE4[(length(methods)+1),j] <- error(mean(datause[datause[,z]==1,d], na.rm = TRUE), dataout[!is.na(dataout[dataout[,z]==1,d]),d])$err
        MSE5[(length(methods)+1),j] <- error(mean(datause[,z], na.rm = TRUE), dataout[!is.na(dataout[,z]),z])$err
        
        
      }
      
      
      if(est=="interactive" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, x=x, method=methods[k], plinear=0, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)}
        else{                       cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, method=methods[k], plinear=0, xL=xL, binary=binary, arguments=arguments);  }
        
        MSE1[k,j]               <- cond.comp[[k,j]]$err.yz0
        MSE2[k,j]               <- cond.comp[[k,j]]$err.yz1
        MSE3[k,j]               <- cond.comp[[k,j]]$err.z
        
        drop                   <- which(cond.comp[[k,j]]$mz_x>trim[1] & cond.comp[[k,j]]$mz_x<trim[2])      
        mz_x                   <- cond.comp[[k,j]]$mz_x[drop]
        my_z1x                 <- cond.comp[[k,j]]$my_z1x[drop]
        my_z0x                 <- cond.comp[[k,j]]$my_z0x[drop]
        yout                   <- dataout[drop,y]
        dout                   <- dataout[drop,d]
        
        TE[1,k]                <- ATE(yout, dout, my_z1x, my_z0x, mz_x)/K + TE[1,k];
        STE[1,k]               <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,k];
        
        ypool[[k]]             <- c(ypool[[k]], yout)
        dpool[[k]]             <- c(dpool[[k]], dout)
        zpool[[k]]             <- c(zpool[[k]], mz_x)
        z1pool[[k]]            <- c(z1pool[[k]], my_z1x)
        z0pool[[k]]            <- c(z0pool[[k]], my_z0x)
        
        
        MSE1[(length(methods)+1),j] <- error(mean(datause[datause[,d]==0,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,d]==0,y]),y])$err
        MSE2[(length(methods)+1),j] <- error(mean(datause[datause[,d]==1,y], na.rm = TRUE), dataout[!is.na(dataout[dataout[,d]==1,y]),y])$err
        MSE3[(length(methods)+1),j] <- error(mean(datause[,d], na.rm = TRUE), dataout[!is.na(dataout[,d]),d])$err
        
      }
      
      if(est=="plinear" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, x=x, method=methods[k], plinear=1, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)}
        if(methods[k]=="Lasso") {
          cond.comp[[k,j]] <- cond_comp(datause=datause_alt, dataout=dataout_alt, y=y, d=d, z=z, x=x_alt, method=methods[k], plinear=1, xL=xL_alt, binary=binary, arguments=arguments)
        }
        else{                        cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, method=methods[k], plinear=1, xL=xL, binary=binary, arguments=arguments)}
        
        
        MSE1[k,j]              <- cond.comp[[k,j]]$err.y
        
        for(i in 1:length(d)){
          MSE2[[i]][k,j]              <- cond.comp[[k,j]]$err.z[[i]]
        }
        
        
        df <- data.frame(matrix(unlist(cond.comp[[k,j]]$rz), ncol=length(cond.comp[[k,j]]$rz), byrow=F))
        
        # lm.fit.ry              <- lm(as.matrix(cond.comp[[k,j]]$ry) ~ as.matrix(cond.comp[[k,j]]$rz)-1);
        #lm.fit.ry              <- lm(as.matrix(cond.comp[[k,j]]$ry) ~ as.matrix(cond.comp[[k,j]]$rz[[1]]) + as.matrix(cond.comp[[k,j]]$rz[[2]]) + as.matrix(cond.comp[[k,j]]$rz[[3]]) -1);
        lm.fit.ry              <- lm(as.matrix(cond.comp[[k,j]]$ry) ~.-1,df);
        ate                    <- lm.fit.ry$coef;
        HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
        
        ypool[[k]]             <- c(ypool[[k]], cond.comp[[k,j]]$ry)
        MSE1[(length(methods)+1),j] <- error(rep(mean(datause[,y], na.rm = TRUE), length(dataout[!is.na(dataout[,y]),y])), dataout[!is.na(dataout[,y]),y])$err
        
        for(i in 1:length(d)){
          STE[[i]][1,k]               <- (1/(K^2))*(diag(HCV.coefs)[i]) +  STE[[i]][1,k]  
          TE[[i]][1,k]                <- ate[i]/K + TE[[i]][1,k] ;
          zpool[[i]][[k]]             <- c(zpool[[i]][[k]], cond.comp[[k,j]]$rz[[i]])
          MSE2[[i]][(length(methods)+1),j] <- error(rep(mean(datause[,d[[i]]], na.rm = TRUE), length(dataout[!is.na(dataout[,d[[i]]]),d[[i]]])), dataout[!is.na(dataout[,d[[i]]]),d[[i]]])$err
          
        }
        
        
      }
      
      if(est=="IV" && (length(methods)>0)){
        
        if(methods[k]=="Ensemble") { 
          cond.comp[[k,j]] <- ensembleF(datause=datause, dataout=dataout, y=y, d=d, z=z, xx=x, method=methods[k], plinear=2, xL=xL, binary=binary, arguments=arguments, ensemble=ensemble)
        }
        if(methods[k]=="Lasso") {
          cond.comp[[k,j]] <- cond_comp(datause=datause_alt, dataout=dataout_alt, y=y, d=d, z=z, x=x_alt, method=methods[k], plinear=2, xL=xL_alt, binary=binary, arguments=arguments)
        }
        else{       
          cond.comp[[k,j]] <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, z=z, x=x, method=methods[k], plinear=2, xL=xL, binary=binary, arguments=arguments)
        }
        
        
        MSE1[k,j]              <- cond.comp[[k,j]]$err.y
        
        for(i in 1:length(d)){
          MSE2[[i]][k,j]              <- cond.comp[[k,j]]$err.z[[i]]
        }
        
        for(v in 1:length(z)){
          MSE3[[v]][k,j]              <- cond.comp[[k,j]]$err.z2[[v]]
        }
        
        df <- data.frame(matrix(unlist(cond.comp[[k,j]]$rz), ncol=length(cond.comp[[k,j]]$rz), byrow=F))
        
        dfz2 <- data.frame(matrix(unlist(cond.comp[[k,j]]$rz2), ncol=length(cond.comp[[k,j]]$rz2), byrow=F))
        
        
        lm.fit.ry              <- tsls(y=cond.comp[[k,j]]$ry,d=as.matrix(df), x=NULL, z=as.matrix(dfz2), intercept = FALSE)
        ate                    <- lm.fit.ry$coef
        
        HCV.coefs              <- lm.fit.ry$vcov
        
        ypool[[k]]             <- c(ypool[[k]], cond.comp[[k,j]]$ry)
        MSE1[(length(methods)+1),j] <- NA
        
        for(i in 1:length(d)){
          STE[[i]][1,k]               <- (1/(K^2))*(diag(HCV.coefs)[i]) +  STE[[i]][1,k]  
          TE[[i]][1,k]                <- ate[i]/K + TE[[i]][1,k] ;
          dpool[[i]][[k]]             <- c(dpool[[i]][[k]], cond.comp[[k,j]]$rz[[i]])
          MSE2[[i]][(length(methods)+1),j] <- NA
        }
        
        for(v in 1:length(z)){
          zpool[[v]][[k]]             <- c(zpool[[v]][[k]], cond.comp[[k,j]]$rz2[[v]])
          MSE3[[v]][(length(methods)+1),j] <- NA
          
        }
        
        
      }
    }  
  }
  
  
  
  if(est=="LATE"){
    
    if(length(methods)>1){
      
      min1 <- which.min(rowMeans(MSE1[1:length(methods),]))
      min2 <- which.min(rowMeans(MSE2[1:length(methods),]))
      min3 <- which.min(rowMeans(MSE3[1:length(methods),]))
      min4 <- which.min(rowMeans(MSE4[1:length(methods),]))
      min5 <- which.min(rowMeans(MSE5[1:length(methods),]))
      
    }
    
    if(length(methods)==1){
      
      min1 <- which.min(mean(MSE1[1:length(methods),]))
      min2 <- which.min(mean(MSE2[1:length(methods),]))
      min3 <- which.min(mean(MSE3[1:length(methods),]))
      min4 <- which.min(mean(MSE4[1:length(methods),]))
      min5 <- which.min(mean(MSE5[1:length(methods),]))
      
    }
    
    if(silent==FALSE){
      cat('  best methods for E[Y|X, D=0]:',methods[min1],'\n')
      cat('  best methods for E[Y|X, D=1]:',methods[min2],'\n')
      cat('  best methods for E[D|X, Z=0]:',methods[min3],'\n')
      cat('  best methods for E[D|X, Z=1]:',methods[min4],'\n')
      cat('  best methods for E[Z|X]:',methods[min5],'\n')
    }
  }
  
  
  if(est=="interactive"){
    
    if(length(methods)>1){
      
      min1 <- which.min(rowMeans(MSE1[1:length(methods),]))
      min2 <- which.min(rowMeans(MSE2[1:length(methods),]))
      min3 <- which.min(rowMeans(MSE3[1:length(methods),]))
      
    }
    
    if(length(methods)>1){
      
      min1 <- which.min(mean(MSE1[1:length(methods),]))
      min2 <- which.min(mean(MSE2[1:length(methods),]))
      min3 <- which.min(mean(MSE3[1:length(methods),]))
      
    }
    
    if(silent==FALSE){
      cat('  best methods for E[Y|X, D=0]:',methods[min1],'\n')
      cat('  best methods for E[Y|X, D=1]:',methods[min2],'\n')
      cat('  best methods for E[D|X]:',methods[min3],'\n')
    }
  }
  
  if(est=="plinear"){
    
    if(length(methods)>1){
      
      min1 <- which.min(rowMeans(MSE1[1:length(methods),]))
      min2 <- rep(0, length(d))
      for(i in 1: length(d)){
        min2[i] <- which.min(rowMeans(MSE2[[i]][1:length(methods),]))
      }
      
      
    }
    
    if(length(methods)==1){
      min1 <- which.min(mean(MSE1[1:length(methods),]))
      min2 <- rep(0, length(d))
      for(i in 1:length(d)){
        min2[i] <- which.min(mean(MSE2[[i]][1:length(methods),]))
      }
      
      
    }
    
    if(silent==FALSE){   
      cat('  best methods for E[Y|X]:',methods[min1],'\n')
      for(i in 1:length(d)){
        cat('  best methods for E[D|X]:',i, methods[min2[i]],'\n')
      }
    }    
  }
  
  if(est=="IV"){
    
    if(length(methods)>1){
      
      min1 <- which.min(rowMeans(MSE1[1:length(methods),]))
      min2 <- rep(0, length(d))
      for(i in 1: length(d)){
        min2[i] <- which.min(rowMeans(MSE2[[i]][1:length(methods),]))
      }
      min3 <- rep(0, length(z))
      for(v in 1:length(z)){
        min3[v] <- which.min(rowMeans(MSE3[[v]][1:length(methods),]))
      }
      
    }
    
    if(length(methods)==1){
      
      min1 <- which.min(mean(MSE1[1:length(methods),]))
      min2 <- rep(0, length(d))
      for(i in 1:length(d)){
        min2[i] <- which.min(mean(MSE2[[i]][1:length(methods),]))
      }
      min3 <- rep(0, length(z))
      for(v in 1:length(z)){
        min3[v] <- which.min(mean(MSE3[[v]][1:length(methods),]))
      }
      
    }
    
    if(silent==FALSE){   
      cat('  best methods for E[Y|X]:',methods[min1],'\n')
      for(i in 1:length(d)){
        cat('  best methods for E[D|X]:',i, methods[min2[i]],'\n')
      }
      for(v in 1:length(z)){
        cat('  best methods for E[Z|X]:', v, methods[min3[[v]]],'\n')
      }
    }    
  }
  
  for(j in 1:K){  
    
    ii = cvgroup == j
    nii = cvgroup != j
    
    datause = as.data.frame(data[nii,])
    dataout = as.data.frame(data[ii,])  
    
    
    if(est=="LATE"){
      
      drop                   <- which(cond.comp[[min5,j]]$mz2_x>trim[1] & cond.comp[[min5,j]]$mz2_x<trim[2])      
      mz_x                   <- cond.comp[[min1,j]]$mz2_x[drop]
      my_z1x                 <- cond.comp[[min2,j]]$my_z1x[drop]
      my_z0x                 <- cond.comp[[min3,j]]$my_z0x[drop]
      md_z1x                 <- cond.comp[[min4,j]]$md_z1x[drop]
      if(flag==1){ md_z0x    <- matrix(0,1,length(my_z0x))}
      else{  md_z0x          <- cond.comp[[min5,j]]$md_z0x[drop] }
      
      yout                   <- dataout[drop,y]
      dout                   <- dataout[drop,d]
      zout                   <- dataout[drop,z]
      
      TE[1,(k+1)]            <- LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)/K + TE[1,(k+1)];
      STE[1,(k+1)]           <- (1/(K^2))*((SE.LATE(yout, dout, zout, my_z1x, my_z0x, mz_x, md_z1x, md_z0x))^2) + STE[1,(k+1)];
      
      ypool[[k+1]]             <- c(ypool[[k+1]], yout)
      dpool[[k+1]]             <- c(dpool[[k+1]], dout)
      zopool[[k+1]]            <- c(zopool[[k+1]], zout)
      zpool[[k+1]]             <- c(zpool[[k+1]], mz_x)
      z1pool[[k+1]]            <- c(z1pool[[k+1]], my_z1x)
      z0pool[[k+1]]            <- c(z0pool[[k+1]], my_z0x)
      dz1pool[[k+1]]           <- c(dz1pool[[k+1]], md_z1x)
      dz0pool[[k+1]]           <- c(dz0pool[[k+1]], md_z0x)
      
    }
    
    
    if(est=="interactive"){
      
      drop                   <- which(cond.comp[[min3,j]]$mz_x>trim[1] & cond.comp[[min3,j]]$mz_x<trim[2])      
      mz_x                   <- cond.comp[[min1,j]]$mz_x[drop]
      my_z1x                 <- cond.comp[[min2,j]]$my_z1x[drop]
      my_z0x                 <- cond.comp[[min3,j]]$my_z0x[drop]
      yout                   <- dataout[drop,y]
      dout                   <- dataout[drop,d]
      
      TE[1,(k+1)]            <- ATE(yout, dout, my_z1x, my_z0x, mz_x)/K + TE[1,(k+1)];
      STE[1,(k+1)]           <- (1/(K^2))*((SE.ATE(yout, dout, my_z1x, my_z0x, mz_x))^2) + STE[1,(k+1)];
      
      ypool[[k+1]]             <- c(ypool[[k+1]], yout)
      dpool[[k+1]]             <- c(dpool[[k+1]], dout)
      zpool[[k+1]]             <- c(zpool[[k+1]], mz_x)
      z1pool[[k+1]]            <- c(z1pool[[k+1]], my_z1x)
      z0pool[[k+1]]            <- c(z0pool[[k+1]], my_z0x)
      
    }
    
    if(est=="plinear"){
      
      best_df <- rep(list(NULL), length(d))
      for(i in 1:length(d)){
        best_df[[i]] <- cond.comp[[min2[i],j]]$rz[[i]]
      }
      
      best_df <- data.frame(matrix(unlist(best_df), ncol=length(d), byrow=F))  
      
      lm.fit.ry              <- lm(as.matrix(cond.comp[[min1,j]]$ry)~.-1, best_df);
      lm.fit.ry              <- lm(as.matrix(cond.comp[[min1,j]]$ry)~as.matrix(best_df)-1);
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
      
      for(i in 1:length(d)){
        STE[[i]][1,(k+1)]           <- (1/(K^2))*(diag(HCV.coefs)[i]) +  STE[[i]][1,(k+1)] 
        TE[[i]][1,(k+1)]            <- ate[i]/K + TE[[i]][1,(k+1)] ;
        MSE2[[i]][(length(methods)+1),j] <- MSE2[[i]][min2[i], j]
        zpool[[i]][[k+1]]             <- c(zpool[[i]][[k+1]], cond.comp[[min2[i],j]]$rz[[i]])
      }
      
      
      #change 
      MSE1[(length(methods)+1),j] <- MSE1[min1, j]
      ypool[[k+1]]             <- c(ypool[[k+1]], cond.comp[[min1,j]]$ry)
      
      
    }
    
    if(est=="IV"){
      best_df <- rep(list(NULL), length(d))
      for(i in 1:length(d)){
        best_df[[i]] <- cond.comp[[min2[i],j]]$rz[[i]]
      }
      
      best_df <- data.frame(matrix(unlist(best_df), ncol=length(d), byrow=F))  
      
      best_dfz2 <- rep(list(NULL), length(z))
      for(v in 1:length(z)){
        best_dfz2[[v]] <- cond.comp[[min3[v],j]]$rz2[[v]]
      }
      
      best_dfz2 <- data.frame(matrix(unlist(best_dfz2), ncol=length(z), byrow=F))  
      
      
      lm.fit.ry              <- tsls(y=cond.comp[[min1,j]]$ry, d=as.matrix(best_df), x=NULL, z=as.matrix(best_dfz2), intercept = FALSE)
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- lm.fit.ry$vcov
      
      ypool[[k+1]]             <- c(ypool[[k+1]], cond.comp[[min1,j]]$ry)
      
      for(i in 1:length(d)){
        STE[[i]][1,(k+1)]           <- (1/(K^2))*(diag(HCV.coefs)[i]) +  STE[[i]][1,(k+1)] 
        TE[[i]][1,(k+1)]            <- ate[i]/K + TE[[i]][1,(k+1)] ;
        dpool[[i]][[k+1]]             <- c(dpool[[i]][[k+1]], cond.comp[[min2[i],j]]$rz[[i]])
      }
      for(v in 1:length(z)){
        zpool[[v]][[k+1]]             <- c(zpool[[v]][[k+1]], cond.comp[[min3[i],j]]$rz2[[v]])
      }
    }
    
  }  
  
  TE_pool        <- rep(list(matrix(0,1,(length(methods)+1))), length(d))
  STE_pool       <- rep(list(matrix(0,1,(length(methods)+1))), length(d))
  
  for(k in 1:(length(methods)+1)){ 
    
    if(est=="LATE"){
      
      TE_pool[1,(k)]         <- LATE(ypool[[k]], dpool[[k]], zopool[[k]],  z1pool[[k]], z0pool[[k]], zpool[[k]], dz1pool[[k]], dz0pool[[k]])
      STE_pool[1,(k)]        <- ((SE.LATE(ypool[[k]], dpool[[k]], zopool[[k]],  z1pool[[k]], z0pool[[k]], zpool[[k]], dz1pool[[k]], dz0pool[[k]]))^2) 
      
    }
    
    if(est=="interactive"){
      
      TE_pool[1,(k)]         <- ATE(ypool[[k]], dpool[[k]], z1pool[[k]], z0pool[[k]], zpool[[k]])
      STE_pool[1,(k)]        <- ((SE.ATE(ypool[[k]], dpool[[k]], z1pool[[k]], z0pool[[k]], zpool[[k]]))^2)
      
    }
    
    if(est=="plinear"){
      
      pool_df <- rep(list(NULL), length(d))
      for(i in 1:length(d)){
        pool_df[[i]] <- zpool[[i]][[k]]
      }
      pool_df <- data.frame(matrix(unlist(pool_df), ncol=length(d), byrow=F))  
      
      lm.fit.ry              <- lm(as.matrix(ypool[[k]]) ~ .-1, pool_df);
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- vcovHC(lm.fit.ry, type = 'HC');
      for(i in 1:length(d)){
        STE_pool[[i]][1,k]          <- (diag(HCV.coefs))[i]
        TE_pool[[i]][1,k]           <- ate[i]
      }
      
      
    }
    
    if(est=="IV"){
      pool_df <- rep(list(NULL), length(d))
      for(i in 1:length(d)){
        pool_df[[i]] <- dpool[[i]][[k]]
      }
      pool_df <- data.frame(matrix(unlist(pool_df), ncol=length(d), byrow=F))  
      
      pool_dfz2 <- rep(list(NULL), length(z))
      for(v in 1:length(z)){
        pool_dfz2[[v]] <- zpool[[v]][[k]]
      }
      pool_dfz2 <- data.frame(matrix(unlist(pool_dfz2), ncol=length(z), byrow=F)) 
      
      lm.fit.ry              <- tsls(y=ypool[[k]],d=as.matrix(pool_df), x=NULL, z=as.matrix(pool_dfz2), intercept = FALSE)
      ate                    <- lm.fit.ry$coef;
      HCV.coefs              <- lm.fit.ry$vcov
      
      for(i in 1:length(d)){
        STE_pool[[i]][1,k]          <- (diag(HCV.coefs))[i]
        TE_pool[[i]][1,k]           <- ate[i]
      }
      
    }
  }
  
  if(length(methods)==1){
    for(i in 1:length(d)){
      TE_pool[[i]][1,(length(methods)+1)]  <- TE_pool[[i]][1,(length(methods))]
      STE_pool[[i]][1,(length(methods)+1)] <- STE_pool[[i]][1,(length(methods))]
    }
  }
  table <- rep(list(NULL), length(d))
  for(i in 1:length(d)){
    colnames(result[[i]])   <- c(methods, "best") 
    colnames(result2[[i]])  <- c(methods, "best") 
    rownames(MSE1)     <- c(methods, "best") 
    rownames(MSE2[[i]])     <- c(methods, "best") 
    rownames(result[[i]])   <- c(paste("ATE", d[[i]]), "se")
    rownames(result2[[i]])  <- c(paste("ATE", d[[i]]), "se")
    
    if(DML=="DML1"){
      result[[i]][1,]         <- colMeans(TE[[i]])
      result[[i]][2,]         <- sqrt((STE[[i]]))
    }
    
    if(DML=="DML2"){
      result[[i]][1,]         <- colMeans(TE_pool[[i]])
      result[[i]][2,]         <- sqrt((STE_pool[[i]]))
    }
    
    
    if(est=="plinear"){   
      table[[i]] <- rbind(result[[i]], rowMeans(MSE1), rowMeans(MSE2[[i]])) 
      rownames(table[[i]])[3:4]   <- c("MSE[Y|X]", "MSE[D|X]") 
    }
    if(est=="IV"){   
      MSE_bind <- matrix(unlist(lapply(MSE3, rowMeans)), nrow=length(z), byrow=T)
      table[[i]] <- rbind(result[[i]], rowMeans(MSE1), rowMeans(MSE2[[i]]), MSE_bind) 
      rownames(table[[i]])[3:4]   <- c("MSE[Y|X]", "MSE[D|X]") 
      for(v in 1:length(z)){
        rownames(MSE3[[v]])     <- c(methods, "best") 
        rownames(table[[i]])[4+v]   <- paste("MSE[Z|X]", z[[v]])
      }
    }
  }
  
  if(est=="interactive"){    
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3))   
    rownames(table)[3:5]   <- c("MSE[Y|X, D=0]", "MSE[Y|X, D=1]", "MSE[D|X]")
  }  
  
  if(est=="LATE"){    
    table <- rbind(result, rowMeans(MSE1), rowMeans(MSE2) , rowMeans(MSE3),rowMeans(MSE4) , rowMeans(MSE5))    
    rownames(table)[3:7]   <- c("MSE[Y|X, Z=0]", "MSE[Y|X, Z=1]", "MSE[D|X,Z=0]", "MSE[D|X,Z=1]" ,"MSE[Z|X]")
  }  
  for(i in 1:length(d)){
    colnames(table[[i]])[length(methods)+1] = "best"
    table[[i]] <- data.frame(t(table[[i]][1,]), t(table[[i]][2,]), t(table[[i]][3,]), t(table[[i]][4,]))
    names(table)[i] <- d[[i]]
  }
  
  
  
  
  if(TRUE %in% (methods %in% c("RLasso"))){
    lasso_index <- which(methods %in% c("RLasso"))
    for(i in 1:nfold){
      lasso_y_coefficients[[i]] <- cond.comp[[lasso_index,i]]$fit.y$model$coefficients
      if(!FALSE %in% c(lasso_y_coefficients[[i]] == 0)){lasso_y_coefficients[[i]]["intercept"] <- cond.comp[[lasso_index,i]]$fit.y$model$intercept} 
      
      lasso_d_coefficients[[i]] <- cond.comp[[lasso_index,i]]$fit.z$model$coefficients
      if(!FALSE %in% c(lasso_d_coefficients[[i]] == 0)){lasso_d_coefficients[[i]]["intercept"] <- cond.comp[[lasso_index,i]]$fit.z$model$intercept} 
      
      if(!is.null(z)){
        lasso_z_coefficients[[i]] <- cond.comp[[lasso_index,i]]$fit.z2$model$coefficients 
        if(!FALSE %in% c(lasso_z_coefficients[[i]] == 0)){lasso_z_coefficients[[i]]["intercept"] <- cond.comp[[lasso_index,i]]$fit.z2$model$intercept} 
      }
    }
  }
  
  if(TRUE %in% (methods %in% c("Lasso"))){
    lasso_index <- which(methods %in% c("Lasso"))
    temp_y <- 0
    for(i in 1:nfold){
      temp_y <- temp_y + coef(cond.comp[[lasso_index,i]]$fit.y$model$glmnet.fit, s= cond.comp[[lasso_index,i]]$fit.y$model$lambda.min)
    }
    lasso_y_coefficients <- list(temp_y/nfold)
    
    temp_d <- rep(list(0), length(d))
    for(u in 1:length(d)){
      for(i in 1:nfold){
        temp_d[[u]] <- temp_d[[u]] + coef(cond.comp[[lasso_index,i]]$fit.z[[u]]$model$glmnet.fit, s= cond.comp[[lasso_index,i]]$fit.z[[u]]$model$lambda.min)
      }
      temp_d[[u]] <- temp_d[[u]]/nfold 
    }
    lasso_d_coefficients <- temp_d
    
    temp_z <- rep(list(0), length(z))
    if(!is.null(z)){
      for(v in 1:length(z)){
        for(i in 1:nfold){
          temp_z[[v]] <- temp_z[[v]] + coef(cond.comp[[lasso_index,i]]$fit.z2[[v]]$model$glmnet.fit, s= cond.comp[[lasso_index,i]]$fit.z2[[v]]$model$lambda.min)
        }
        temp_z[[v]] <- temp_z[[v]]/nfold 
      }
      lasso_z_coefficients <- temp_z
    }
    else{
      lasso_z_coefficients <- NA
    }
  }
  return(list("table" = table, "lasso_y_coeff" = lasso_y_coefficients, "lasso_d_coeff" = lasso_d_coefficients , "lasso_z_coeff" = lasso_z_coefficients))
}  

cond_comp <- function(datause, dataout, y, d, z=NULL, x, method, plinear,xL, binary,flag=0, arguments){
  
  form_y   <- y
  form_d   <- d
  form_z   <- z
  form_x   <- x
  form_xL  <- xL
  ind_u <- rep(list(0), length(d))
  ind_o <- rep(list(0), length(d))
  for(i in 1: length(d)){
    ind_u[[i]] <- which(datause[,d[[i]]]==1)
    ind_o[[i]] <- which(dataout[,d[[i]]]==1)
  }
  if(plinear==3){
    ind_u    <- which(datause[,z]==1)
    ind_o    <- which(dataout[,z]==1)
  }
  err.yz1  <- NULL
  err.yz0  <- NULL
  my_z1x   <- NULL
  my_z0x   <- NULL
  fit.yz1  <- NULL
  fit.yz0  <- NULL
  fit.dz1  <- NULL
  err.dz1  <- NULL
  md_z1x   <- NULL
  fit.dz0  <- NULL
  err.dz0  <- NULL
  md_z0x   <- NULL
  fit.y <-NULL
  
  #D
  fit.z <- rep(list(NULL), length(d))
  rz <- rep(list(NULL), length(d))
  mz_x <- rep(list(NULL), length(d))
  err.z <- rep(list(NULL), length(d))
  mis.z <- rep(list(NULL), length(d))
  
  #Z
  fit.z2   <- rep(list(NULL), length(z))
  rz2      <- rep(list(NULL), length(z))
  mz2_x    <- rep(list(NULL), length(z))
  err.z2   <- rep(list(NULL), length(z))
  mis.z2 <- rep(list(NULL), length(z))
  
  ########################## Boosted  Trees ###################################################;
  
  if(method=="Boosting"){
    
    option <- arguments[[method]]
    arg    <- option
    arg[which(names(arg) %in% c("clas_dist","reg_dist"))] <-  NULL
    
    if(plinear==3){
      
      fit.yz1        <- boost(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, n.trees=fit.yz1$best, dataout, type="response") 
      
      fit.dz1        <- boost(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_d, distribution=option[['clas_dist']], option=arg)
      err.dz1        <- error(fit.dz1$yhatout, dataout[ind_o,y])$err
      md_z1x         <- predict(fit.dz1$model, n.trees=fit.dz1$best, dataout, type="response") 
      
      fit.yz0        <- boost(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,d])$err
      my_z0x         <- predict(fit.yz0$model, n.trees=fit.yz0$best, dataout, type="response") 
      
      if(flag==0){
        fit.dz0        <- boost(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_d, distribution=option[['clas_dist']], option=arg)
        err.dz0        <- error(fit.dz0$yhatout, dataout[-ind_o,d])$err
        md_z0x         <- predict(fit.dz0$model, n.trees=fit.dz0$best, dataout, type="response") 
      }
    }
    
    if(plinear==0){
      
      fit.yz1        <- boost(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, n.trees=fit.yz1$best, dataout, type="response") 
      
      fit.yz0        <- boost(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, n.trees=fit.yz0$best, dataout, type="response") 
      
      
    }
    for(i in 1:length(d)){
      if(binary[i]==1){
        fit.z[[i]]          <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]], distribution=option[['clas_dist']], option=arg)
        mis.z[[i]]          <- error(fit.z[[i]]$yhatout, dataout[,d[[i]]])$mis
      }
      
      
      if(binary[i]==0){
        fit.z[[i]]          <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]], distribution=option[['reg_dist']], option=arg)
        mis.z[[i]]          <- NA
      }
    }
    fit.y            <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, distribution=option[['reg_dist']], option=arg)
    
    if(plinear==2 | plinear==3){
      for(v in 1:length(z)){
        fit.z2[[v]]         <- boost(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_z[[v]], distribution=option[['reg_dist']], option=arg)
      }
    }
    
  }  
  
  
  ########################## Neural Network(Nnet Package) ###################################################;   
  
  
  if(method=="Nnet"){
    
    option <- arguments[[method]]
    arg    <- option
    
    if(plinear==3){
      
      fit.yz1        <- nnetF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz1$f] <- as.data.frame(scale(dataouts[,!fit.yz1$f], center = fit.yz1$min, scale = fit.yz1$max - fit.yz1$min))
      my_z1x         <- predict(fit.yz1$model, dataouts)*(fit.yz1$max[fit.yz1$k]-fit.yz1$min[fit.yz1$k])+fit.yz1$min[fit.yz1$k] 
      
      fit.yz0        <- nnetF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz0$f] <- as.data.frame(scale(dataouts[,!fit.yz0$f], center = fit.yz0$min, scale = fit.yz0$max - fit.yz0$min))
      my_z0x         <- predict(fit.yz0$model, dataouts)*(fit.yz0$max[fit.yz0$k]-fit.yz0$min[fit.yz0$k])+fit.yz0$min[fit.yz0$k] 
      
      fit.dz1        <- nnetF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_d, arg=arg)
      err.dz1        <- error(fit.dz1$yhatout, dataout[ind_o,d])$err
      dataouts       <- dataout
      dataouts[,!fit.dz1$f] <- as.data.frame(scale(dataouts[,!fit.dz1$f], center = fit.dz1$min, scale = fit.dz1$max - fit.dz1$min))
      md_z1x         <- predict(fit.dz1$model, dataouts)*(fit.dz1$max[fit.dz1$k]-fit.dz1$min[fit.dz1$k])+fit.yz1$min[fit.dz1$k] 
      
      if(flag==0){
        fit.dz0        <- nnetF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_d, arg=arg)
        err.dz0        <- error(fit.dz0$yhatout, dataout[-ind_o,d])$err
        dataouts       <- dataout
        dataouts[,!fit.dz0$f] <- as.data.frame(scale(dataouts[,!fit.dz0$f], center = fit.dz0$min, scale = fit.dz0$max - fit.dz0$min))
        md_z0x         <- predict(fit.dz0$model, dataouts)*(fit.dz0$max[fit.dz0$k]-fit.dz0$min[fit.dz0$k])+fit.dz0$min[fit.dz0$k] 
      }
    }
    
    if(plinear==0){
      
      fit.yz1        <- nnetF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz1$f] <- as.data.frame(scale(dataouts[,!fit.yz1$f], center = fit.yz1$min, scale = fit.yz1$max - fit.yz1$min))
      my_z1x         <- predict(fit.yz1$model, dataouts)*(fit.yz1$max[fit.yz1$k]-fit.yz1$min[fit.yz1$k])+fit.yz1$min[fit.yz1$k] 
      
      fit.yz0        <- nnetF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      dataouts       <- dataout
      dataouts[,!fit.yz0$f] <- as.data.frame(scale(dataouts[,!fit.yz0$f], center = fit.yz0$min, scale = fit.yz0$max - fit.yz0$min))
      my_z0x         <- predict(fit.yz0$model, dataouts)*(fit.yz0$max[fit.yz0$k]-fit.yz0$min[fit.yz0$k])+fit.yz0$min[fit.yz0$k] 
    }
    for(i in 1:length(d)){
      if(binary[i]==1){
        fit.z[[i]]          <- nnetF(datause=datause, dataout=dataout, form_x=form_x, form_y=form_d[[i]], clas=TRUE, arg=arg)
        mis.z[[i]]          <- error(fit.z$yhatout, dataout[,d[[i]]])$mis
      }
      
      if(binary[i]==0){
        fit.z[[i]]          <- nnetF(datause=datause, dataout=dataout, form_x=form_x, form_y=form_d[[i]], clas=FALSE, arg=arg)
        mis.z[[i]]          <- NA
      }
    }
    
    
    fit.y          <- nnetF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, arg=arg)
    
    if(plinear==2 | plinear==3){
      for(v in 1:length(z)){
        fit.z2[[v]]           <- nnetF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_z[[v]], arg=arg)
      }
    }
    
  } 
  
  ########################## Lasso and Post Lasso(Hdm Package) ###################################################;    
  
  if(method=="RLasso" || method=="PostRLasso"){
    
    post = FALSE
    if(method=="PostRLasso"){ post=TRUE }
    
    option    <- arguments[[method]]
    arg       <- option
    
    if(plinear==3){
      
      fit.yz1        <- rlassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, post, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, newdata=formC(form_y, form_x, dataout)$x , type="response") 
      
      fit.yz0        <- rlassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, post, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, newdata=formC(form_y, form_x, dataout)$x, type="response")   
      
      fit.dz1        <- rlassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_d, post, arg=arg)
      err.dz1        <- error(fit.dz1$yhatout, dataout[ind_o,d])$err
      md_z1x         <- predict(fit.dz1$model, newdata=formC(form_y, form_x, dataout)$x , type="response") 
      
      if(flag==0){
        fit.dz0        <- rlassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_d, post, arg=arg)
        err.dz0        <- error(fit.dz0$yhatout, dataout[-ind_o,d])$err
        md_z0x         <- predict(fit.dz0$model, newdata=formC(form_y, form_x, dataout)$x, type="response")   
      }
    }
    
    if(plinear==0){
      
      fit.yz1        <- rlassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, post, arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, newdata=formC(form_y, form_x, dataout)$x , type="response") 
      
      fit.yz0        <- rlassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, post, arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, newdata=formC(form_y, form_x, dataout)$x, type="response")   
      
    }
    for(i in 1:length(d))
    {
      if(binary[i]==1){
        fit.z[[i]]          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_d[[i]], post, logit=TRUE, arg=arg)
        mis.z[[i]]          <- error(fit.z[[i]]$yhatout, dataout[,d[[i]]])$mis
      }
      
      if(binary[i]==0){
        fit.z[[i]]          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_d[[i]], post, logit=FALSE, arg=arg)
        mis.z[[i]]          <- NA
      }
    }
    
    
    fit.y          <- rlassoF(datause=datause, dataout=dataout,  form_x, form_y, post, arg=arg)
    
    if(plinear==2 | plinear==3){
      for(v in 1:length(z)){
        fit.z2[[v]]         <- rlassoF(datause=datause, dataout=dataout,  form_x, form_z[[v]], post, arg=arg)
      }
    }
    
  }    
  
  
  ########################## Lasso and Post Lasso(Glmnet) Package) ###################################################;    
  
  if(method=="Ridge" || method=="Lasso" || method=="Elnet"){
    
    if(method=="Ridge"){ alp=0 }
    if(method=="Lasso"){ alp=1 }
    if(method=="Elnet"){ alp=0.5 }
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("s"))] <-  NULL
    s         <- option[['s']]
    
    if(plinear==3){
      
      fit.yz1        <- lassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z1x         <- predict(fit.yz1$model, newx=fit.p$x[,-1] ) 
      
      fit.yz0        <- lassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z0x         <- predict(fit.yz0$model,  newx=fit.p$x[,-1])   
      
      fit.dz1        <- lassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_d, alp=alp, arg=arg, s=s)
      err.dz1        <- error(fit.yz1$yhatout, dataout[ind_o,d])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      md_z1x         <- predict(fit.dz1$model, newx=fit.p$x[,-1] ) 
      
      if(flag==0){
        fit.dz0        <- lassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_d, alp=alp, arg=arg, s=s)
        err.dz0        <- error(fit.yz0$yhatout, dataout[-ind_o,d])$err
        fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
        md_z0x         <- predict(fit.dz0$model,  newx=fit.p$x[,-1])   
      }
    }
    
    if(plinear==0){
      
      fit.yz1        <- lassoF(datause=datause[ind_u,], dataout=dataout[ind_o,],  form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z1x         <- predict(fit.yz1$model, newx=fit.p$x[,-1] ) 
      
      fit.yz0        <- lassoF(datause=datause[-ind_u,], dataout=dataout[-ind_o,], form_x, form_y, alp=alp, arg=arg, s=s)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      fit.p          <- lm(as.formula(paste(form_y, "~", form_x)),  x = TRUE, y = TRUE, data=dataout);
      my_z0x         <- predict(fit.yz0$model,  newx=fit.p$x[,-1])   
      
    }
    
    for(i in 1:length(d)){
      if(binary[i]==1){
        fit.z[[i]]          <- lassoF(datause=datause, dataout=dataout,  form_x, form_d[[i]], logit=TRUE, alp=alp, arg=arg, s=s)
        mis.z[[i]]          <- error(fit.z$yhatout, dataout[,d[[i]]])$mis
      }
      
      if(binary[i]==0){
        fit.z[[i]]          <- lassoF(datause=datause, dataout=dataout,  form_x, form_d[[i]], logit=FALSE, alp=alp, arg=arg, s=s)
        mis.z[[i]]          <- NA
      }   
      
    }
    
    fit.y            <- lassoF(datause=datause, dataout=dataout,  form_x, form_y, alp=alp, arg=arg, s=s)
    
    if(plinear==2 | plinear==3){
      for(v in 1:length(z)){
        fit.z2[[v]]         <- lassoF(datause=datause, dataout=dataout,  form_x, form_z[[v]], alp=alp, arg=arg, s=s)
      }
      
    }
    
  }    
  
  ############# Random Forest ###################################################;
  
  if(method=="Forest" | method=="TForest"){
    
    tune = FALSE
    if(method=="TForest"){tune=TRUE}
    
    option    <- arguments[[method]]
    
    arg       <- option
    arg[which(names(arg) %in% c("clas_nodesize","reg_nodesize"))] <-  NULL
    
    if(plinear==3){
      
      fit.yz1        <- RF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout, type="response") 
      
      fit.yz0        <- RF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout, type="response")
      
      fit.dz1        <- RF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_d, nodesize=option[["clas_nodesize"]], arg=arg, tune=tune)
      err.dz1        <- error(fit.dz1$yhatout, dataout[ind_o,d])$err
      md_z1x         <- predict(fit.dz1$model, dataout, type="response") 
      
      if(flag==0){
        fit.dz0        <- RF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_d, nodesize=option[["clas_nodesize"]], arg=arg, tune=tune)
        err.dz0        <- error(fit.dz0$yhatout, dataout[-ind_o,d])$err
        md_z0x         <- predict(fit.dz0$model, dataout, type="response")
      }
    }
    
    if(plinear==0){
      
      fit.yz1        <- RF(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout, type="response") 
      
      fit.yz0        <- RF(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]], arg=arg, tune=tune)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model, dataout, type="response")
      
    }
    for(i in 1:length(d)){
      if(binary[i]==1){
        fit.z[[i]]          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]], nodesize=option[["clas_nodesize"]], arg=arg, reg=TRUE, tune=tune)
        mis.z[[i]]          <- error(as.numeric(fit.z[[i]]$yhatout), dataout[,d[[i]]])$mis
      }
      
      if(binary[i]==0){
        fit.z[[i]]          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]],nodesize=option[["reg_nodesize"]], arg=arg, reg=TRUE, tune=tune)
        mis.z[[i]]          <- NA
      }   
      
    }
    
    fit.y           <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, nodesize=option[["reg_nodesize"]],  arg=arg, tune=tune)
    
    if(plinear==2 | plinear==3){
      for(v in 1:length(z)){
        fit.z2[[v]]          <- RF(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_z[[v]], nodesize=option[["clas_nodesize"]],  arg=arg, tune=tune)
      }
    }   
    
  }
  
  ########################## Regression Trees ###################################################;     
  
  if(method=="Trees"){
    
    option    <- arguments[[method]]
    arg       <- option
    arg[which(names(arg) %in% c("reg_method","clas_method"))] <-  NULL
    
    if(plinear==3){
      
      fit.yz1        <- tree(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout) 
      
      fit.yz0        <- tree(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model,dataout)   
      
      fit.dz1        <- tree(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_d, method=option[["clas_method"]], arg=arg)
      err.dz1        <- error(fit.dz1$yhatout, dataout[ind_o,d])$err
      md_z1x         <- predict(fit.dz1$model, dataout) 
      
      if(flag==0){
        fit.dz0        <- tree(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_d, method=option[["clas_method"]], arg=arg)
        err.dz0        <- error(fit.dz0$yhatout, dataout[-ind_o,d])$err
        md_z0x         <- predict(fit.dz0$model,dataout)  
      }
      fit.z2         <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_z, method=option[["clas_method"]], arg=arg)
      
    }
    
    if(plinear==0){
      
      fit.yz1        <- tree(datause=datause[ind_u,], dataout=dataout[ind_o,], form_x=form_x,  form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz1        <- error(fit.yz1$yhatout, dataout[ind_o,y])$err
      my_z1x         <- predict(fit.yz1$model, dataout) 
      
      fit.yz0        <- tree(datause=datause[-ind_u,], dataout=dataout[-ind_o,],  form_x=form_x, form_y=form_y, method=option[["reg_method"]], arg=arg)
      err.yz0        <- error(fit.yz0$yhatout, dataout[-ind_o,y])$err
      my_z0x         <- predict(fit.yz0$model,dataout)   
      
    }
    
    for(i in 1:length(d)){
      if(binary[i]==1){
        
        fit.z[[i]]          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]], method=option[["clas_method"]], arg=arg)
        mis.z[[i]]          <- error(as.numeric(fit.z$yhatout[[i]]), dataout[,d[[i]]])$mis
      }
      
      if(binary[i]==0){
        
        fit.z[[i]]          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_d[[i]], method=option[["reg_method"]], arg=arg)
        mis.z[[i]]          <- NA
      }      
    }
    
    
    fit.y           <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_y, method=option[["cont_method"]], arg=arg)
    
    if(plinear==2){
      for(v in 1:length(z)){
        fit.z2[[v]]          <- tree(datause=datause, dataout=dataout,  form_x=form_x, form_y=form_z[[v]], method=option[["cont_method"]], arg=arg)
      }
    }   
  }
  
  for(i in 1: length(d)){
    err.z[[i]]          <- error(fit.z[[i]]$yhatout, dataout[,d[[i]]])$err
    mz_x[[i]]           <- fit.z[[i]]$yhatout       
    rz[[i]]             <- fit.z[[i]]$resout 
  }
  
  ry             <- fit.y$resout
  my_x           <- fit.y$yhatout
  err.y          <- error(fit.y$yhatout, dataout[,y])$err
  
  if(plinear==2 | plinear==3){
    for(v in 1:length(z)){
      rz2[[v]]            <- fit.z2[[v]]$resout
      mz2_x[[v]]          <- fit.z2[[v]]$yhatout
      err.z2[[v]]         <- error(fit.z2[[v]]$yhatout, dataout[,z[[v]]])$err
    }
  }   
  
  return(list(fit.dz0 = fit.dz0, err.dz0=err.dz0, md_z0x=md_z0x, fit.dz1 = fit.dz1, err.dz1=err.dz1, md_z1x=md_z1x, rz2=rz2, mz2_x=mz2_x, err.z2=err.z2, my_z1x=my_z1x, mz_x= mz_x, my_z0x=my_z0x, my_x = my_x, err.z = err.z,  err.yz0= err.yz0,  err.yz1=err.yz1, mis.z=mis.z, ry=ry , rz=rz, err.y=err.y,  fit.y=fit.y, fit.z=fit.z, fit.z2=fit.z2,  fit.yz1out= fit.yz1$yhatout,  fit.yz0out= fit.yz0$yhatout));
  
}  

ensembleF <- function(datause, dataout, y, d, z=NULL, xx, method, plinear, xL, binary, flag=flag, arguments, ensemble){
  
  K         <- 2
  k         <- length(ensemble[['methods']])
  fits      <- vector("list", k)
  method    <- ensemble[['methods']]
  
  fit.z <- rep(list(NULL), length(d)) 
  mis.z <- rep(list(NULL), length(d))
  err.z <- rep(list(NULL), length(d))
  mz_x <- rep(list(NULL), length(d))
  rz <- rep(list(NULL), length(d))
  rz2            <- rep(list(NULL), length(z))
  mz2_x          <- rep(list(NULL), length(z))
  err.z2         <- rep(list(NULL), length(z))
  fit.z2 <- rep(list(NULL), length(z))
  
  ind_u <- rep(list(0), length(d))
  ind_o <- rep(list(0), length(d))
  for(i in 1: length(d)){
    ind_u[[i]] <- which(datause[,d[[i]]]==1)
    ind_o[[i]] <- which(dataout[,d[[i]]]==1)
  }
  
  if(plinear==3){
    ind_u     <- which(datause[,z]==1)
    ind_o     <- which(dataout[,z]==1)
  }
  
  split     <- runif(nrow(datause))
  cvgroup   <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/K)),include.lowest = TRUE))  
  
  if(k<4)  {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.01))) }
  if(k==4) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.02))) }
  if(k==5) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.04))) }
  if(k==6) {  lst <- lapply(numeric(k), function(x) as.vector(seq(0,1,0.1)))  }
  
  gr      <- as.matrix(expand.grid(lst))
  weight  <- gr[rowSums(gr)==1,]
  
  length_error_pred <- length(d) + length(z) + 1
  
  if(plinear==1 | plinear==2){
    
    
    errorM  <- array(0,dim=c(nrow(weight),5,length_error_pred)) #weights, splits, parameters
    pred2   <- array(0,dim=c(nrow(dataout),length_error_pred,k)) #n, parameters, method 
    
    for(j in 1:K){
      
      ii   <- cvgroup == j
      nii  <- cvgroup != j
      
      datause1 <- as.data.frame(datause[nii,])
      datause2 <- as.data.frame(datause[ii,])  
      pred1    <- array(0,dim=c(nrow(datause2),length_error_pred,k))  #length prediction, parameters, methods
      
      for(i in 1:k){
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
          x=xL
        } else {
          x=xx
        }
        
        if(plinear==2){
          fits[[i]]   <- cond_comp(datause=datause1, dataout=datause2, y=y, d=d, z=z, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
        }
        
        if(plinear==1){
          fits[[i]]   <- cond_comp(datause=datause1, dataout=datause2, y=y, d=d, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
        }
        
        pred1[,1,i] <- fits[[i]][['my_x']]
        for(u in 1:length(d)){
          pred1[,1+u,i] <- fits[[i]][['mz_x']][[u]]
        }
        
        
        
        if(plinear==2){
          for(v in 1:length(z)){
            pred1[,1+length(d)+v,i] <- fits[[i]][['mz2_x']][[v]]
          }  
        }
      }
      
      for(p in 1:nrow(weight)){
        errorM[p,j,1] <- error(pred1[,1,] %*% (weight[p,]), datause2[,y])$err #all predcition by one method in a matrix
        for(u in 1:length(d)){
          errorM[p,j,1+u] <- error(pred1[,1+u,] %*% (weight[p,]), datause2[,d[[u]]])$err 
        }
        if(plinear==2){
          for(v in 1:length(z)){
            errorM[p,j,1+length(d)+v] <- error(pred1[,(1+length(d)+v),] %*% (weight[p,]), datause2[,z[[v]]])$err 
          }
          
        }
      }
    }
    
    min1 <- which.min(as.matrix(rowSums(errorM[,,1])))
    min2 <- rep(0, length(d))
    for(u in 1:length(d)){
      min2[u] <- which.min(as.matrix(rowSums(errorM[,,1+u])))
    }
    
    if(plinear==2){
      min3 <- rep(0, length(z))
      for(v in 1:length(z)){
        min3[v] <- which.min(as.matrix(rowSums(errorM[,,1+length(d)+v])))
      }
      
    }
    
    for(i in 1:k){
      
      if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
        x=xL
      } else {
        x=xx
      }
      
      if(plinear==2){
        fits[[i]]   <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, z=z, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
      }
      
      if(plinear==1){
        fits[[i]]   <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
      }
      
      pred2[,1,i] <- fits[[i]][['my_x']]
      for(u in 1:length(d)){
        pred2[,1+u,i] <- fits[[i]][['mz_x']][[u]]  
      }
      
      if(plinear==2){
        for(v in 1:length(z)){
          pred2[,1+length(d)+v,i] <- fits[[i]][['mz2_x']][[v]]  
        }
        
      }
      
    }
    
    fit.y  <- pred2[,1,] %*% (weight[min1,])
    ry     <- dataout[,y] - fit.y
    err.y  <- error(fit.y, dataout[,y])$err 
    
    for(u in 1:length(d)){
      fit.z[[u]]  <- pred2[,1+u,] %*% (weight[min2[u],])
      rz[[u]]     <- dataout[,d[[u]]] - fit.z[[u]]
      err.z[[u]]  <- error(fit.z[[u]], dataout[,d[[u]]])$err
    }
    
    if(plinear==2){
      for(v in 1:length(z)){
        fit.z2[[v]] <- pred2[,1+length(d)+v,] %*% (weight[min3[v],])
        rz2[[v]] <- dataout[,z[[v]]] - fit.z2[[v]]
        err.z2[[v]] <- error(fit.z2[[v]], dataout[,z[[v]]])$err  
      }
      
    }
    
    return(list(err.z=err.z, err.y=err.y, err.z2=err.z2, ry=ry, rz=rz, rz2=rz2));
    
  }
  
  
  if(plinear==3){
    
    errorM  <- array(0,dim=c(nrow(weight),5,5))
    pred2   <- array(0,dim=c(nrow(dataout),5,k))  
    
    for(j in 1:K){
      
      ii   <- cvgroup == j
      nii  <- cvgroup != j
      
      datause1 <- as.data.frame(datause[nii,])
      datause2 <- as.data.frame(datause[ii,])  
      pred1    <- array(0,dim=c(nrow(datause2),5,k))  
      
      ind_u1     <- which(datause1[,d]==1)
      ind_u2     <- which(datause2[,d]==1)
      if(plinear==3){
        ind_u1     <- which(datause1[,z]==1)
        ind_u2     <- which(datause2[,z]==1)
      }
      
      for(i in 1:k){
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
          x=xL
        } else {
          x=xx
        }
        
        fits[[i]]   <- cond_comp(datause=datause1, dataout=datause2, y=y, d=d, z=z, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, flag=flag, arguments=arguments);
        
        pred1[,1,i] <- fits[[i]][['my_z1x']]
        pred1[,2,i] <- fits[[i]][['my_z0x']]
        pred1[,3,i] <- fits[[i]][['md_z1x']]
        if(flag==1){ pred1[,4,i] <- matrix(0,1,length(fits[[i]][['md_z1x']]))}
        else{pred1[,4,i] <- fits[[i]][['md_z0x']]}
        pred1[,5,i] <- fits[[i]][['mz2_x']]
        
      }
      
      for(p in 1:nrow(weight)){
        
        errorM[p,j,1] <- error(pred1[ind_u2,1,] %*% (weight[p,]), datause2[ind_u2,y])$err 
        errorM[p,j,2] <- error(pred1[-ind_u2,2,] %*% (weight[p,]), datause2[-ind_u2,y])$err 
        errorM[p,j,3] <- error(pred1[ind_u2,3,] %*% (weight[p,]), datause2[ind_u2,d])$err 
        errorM[p,j,4] <- error(pred1[-ind_u2,4,] %*% (weight[p,]), datause2[-ind_u2,d])$err 
        errorM[p,j,5] <- error(pred1[,5,] %*% (weight[p,]), datause2[,z])$err 
        
      }
    }
    
    min1 <- which.min(as.matrix(rowSums(errorM[,,1])))
    min2 <- which.min(as.matrix(rowSums(errorM[,,2])))
    min3 <- which.min(as.matrix(rowSums(errorM[,,3])))
    min4 <- which.min(as.matrix(rowSums(errorM[,,4])))
    min5 <- which.min(as.matrix(rowSums(errorM[,,5])))
    
    for(i in 1:k){
      
      if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
        x=xL
      } else {
        x=xx
      }
      
      fits[[i]] <-cond_comp(datause=datause, dataout=dataout, y=y, d=d, z=z, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary,flag=flag, arguments=arguments);
      
      pred2[,1,i] <- fits[[i]][['my_z1x']]
      pred2[,2,i] <- fits[[i]][['my_z0x']]
      pred2[,3,i] <- fits[[i]][['md_z1x']]
      if(flag==1){ pred2[,4,i] <- matrix(0,1,length(fits[[i]][['md_z1x']]))}
      else{pred2[,4,i] <- matrix(0,1,fits[[i]][['md_z0x']])}
      pred2[,5,i] <- fits[[i]][['mz2_x']]
      
    }
    
    fit.yz1x   <- pred2[,1,] %*% (weight[min1,])
    fit.yz0x   <- pred2[,2,] %*% (weight[min2,])
    fit.dz1x   <- pred2[,3,] %*% (weight[min3,])
    fit.dz0x   <- pred2[,4,] %*% (weight[min4,])
    fit.z2x    <- pred2[,5,] %*% (weight[min5,])
    
    err.yz0  <- error(fit.yz1x[-ind_o], dataout[-ind_o,y])$err 
    err.yz1  <- error(fit.yz0x[ind_o], dataout[ind_o,y])$err 
    err.dz0  <- error(fit.dz0x[-ind_o], dataout[-ind_o,d])$err 
    err.dz1  <- error(fit.dz1x[ind_o], dataout[ind_o,d])$err 
    err.z2   <- error(fit.z2x, dataout[,z])$err 
    
    return(list(err.yz0=err.yz0, err.yz1=err.yz1, err.dz0=err.dz0, err.dz1=err.dz1, err.z2=err.z2, my_z0x=fit.yz0x, my_z1x=fit.yz1x, md_z0x=fit.dz0x, md_z1x=fit.dz1x, mz2_x=fit.z2x));
    
  }
  
  
  if(plinear==0){
    
    errorM  <- array(0,dim=c(nrow(weight),5,3))
    pred2   <- array(0,dim=c(nrow(dataout),3,k))  
    pred3   <- matrix(0,nrow(dataout[ind_o,]),k)
    pred4   <- matrix(0,nrow(dataout[-ind_o,]),k)
    
    for(j in 1:K){
      
      ii = cvgroup == j
      nii = cvgroup != j
      
      datause1 = as.data.frame(datause[nii,])
      datause2 = as.data.frame(datause[ii,])  
      
      ind2_u   <- which(datause1[,d]==1)
      ind2_o   <- which(datause2[,d]==1)
      
      pred1  <- array(0,dim=c(nrow(datause2),3,k))  
      
      for(i in 1:k){
        
        if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
          x=xL
        } else {
          x=xx
        }
        
        fits[[i]] <-cond_comp(datause=datause1, dataout=datause2, y=y, d=d, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
        
        pred1[,1,i] <- fits[[i]][['my_z1x']]
        pred1[,2,i] <- fits[[i]][['my_z0x']]
        pred1[,3,i] <- fits[[i]][['mz_x']]
      }
      
      for(p in 1:nrow(weight)){
        
        errorM[p,j,1] <- error(pred1[ind2_o,1,] %*% (weight[p,]), datause2[ind2_o,y])$err 
        errorM[p,j,2] <- error(pred1[-ind2_o,2,] %*% (weight[p,]), datause2[-ind2_o,d])$err 
        errorM[p,j,3] <- error(pred1[,3,] %*% (weight[p,]), datause2[,d])$err 
      }
    }
    
    min1 <- which.min(as.matrix(rowSums(errorM[,,1])))
    min2 <- which.min(as.matrix(rowSums(errorM[,,2])))
    min3 <- which.min(as.matrix(rowSums(errorM[,,3])))
    
    for(i in 1:k){
      
      if (any(c("RLasso", "PostRLasso", "Ridge", "Lasso", "Elnet")==method[i])){
        x=xL
      } else {
        x=xx
      }
      
      fits[[i]]   <- cond_comp(datause=datause, dataout=dataout, y=y, d=d, x=x, method=method[i], plinear=plinear, xL=xL, binary=binary, arguments=arguments);
      
      pred2[,1,i] <- fits[[i]][['my_z1x']]
      pred2[,2,i] <- fits[[i]][['my_z0x']]
      pred2[,3,i] <- fits[[i]][['mz_x']]
      pred3[,i]   <- fits[[i]][['fit.yz1out']]
      pred4[,i]   <- fits[[i]][['fit.yz0out']]
      
    }
    
    fit.y1x    <- pred2[,1,] %*% (weight[min1,])
    fit.y0x    <- pred2[,2,] %*% (weight[min2,])
    fit.zx     <- pred2[,3,] %*% (weight[min3,])
    fit.yz1out <- pred3 %*% (weight[min1,])
    fit.yz0out <- pred4 %*% (weight[min2,])
    
    err.y1x <- error(fit.yz1out, dataout[ind_o,y])$err 
    err.y0x <- error(fit.yz0out, dataout[-ind_o,y])$err 
    err.z   <- error(fit.zx, dataout[,d])$err 
    
    
    return(list(err.z=err.z, err.yz1=err.y1x, err.yz0=err.y0x,  my_z1x=fit.y1x, my_z0x=fit.y0x,mz_x=fit.zx ));
    
  }
}



