library(ExtremeRisks)
library(dplyr)
library(tidyr)

library(fGarch)
library(rugarch)
library(qrmtools)
library(WeightedPortTest)
library(hms)

####
coins<-c("Bitcoin") #indexes are not "traded" during weekends

folder<-"raw_kaggle_data"
df <-as.data.frame(matrix(nrow=0,ncol=6))
colnames(df)<-c("Date", "Name", "Close","logClose", "logReturn","NlogReturn")
setwd("/Users/andreateruzzi/Deskcrypto_tail_expectilestop/data/")


#start data cleaning
for(coin in coins){
  filename<-paste0("coin_",coin,".csv")
  y_raw <-read.csv(file.path(folder,filename))
  y_raw$Date<-as.Date(y_raw$Date)
  
  # set dates in "yyyy-mm-dd" format
  start_date<-as.Date("1950-12-31")
  end_date<-as.Date("2050-12-31")
  min_date<-min(y_raw$Date)
  max_date<-max(y_raw$Date)
  
  if (start_date>end_date){
    start_date<-min_date
    end_date<-max_date
    warning("Start date was later than end date; they have been set to the first and last available dates")
  } 
  else{
    if (start_date<min_date){
      start_date<-min_date
      warning(paste(coin,": start date was earlier than the first available date and has been set to the latter"))
    }
    if (end_date>max_date){
      end_date<-max_date
      warning(paste(coin,": end date was later than the last available date and has been set to the latter"))
    }
  }
  #neglog return
  y_raw$Name = coin
  y <-y_raw %>% select("Date", "Name", "Close") %>% filter(Date >= start_date) %>% filter(Date <= end_date)
  y$logClose <- log(y$Close)
  y$logReturn <- NA
  y$logReturn[-1]<-diff(y$logClose)
  y$NlogReturn <- -y$logReturn
  
  df<-rbind(df,y)
  rm(y)
}
  
  # DATA
  #  coins<-c("Bitcoin", "Ethereum", "Litecoin", "Nasdaq", "DowJones","SP") #indexes are not "traded" during weekends
  coin <- "Bitcoin"
  k_fix=100  #number of obs to consider for the dynamic estimation
  
  y  <- df %>% filter(Name == coin) %>%  select(NlogReturn, Date) 
  Date <- y$Date
  NlogReturn <- y$NlogReturn
  
  # ROLLING WINDOW SIZE
  nsdj<-nrow(y)
  block <- 1000
  
  
  ############ CHECK DIFFERENT TIME SERIES MODELS ############ 
  
  
  #model1: ARMA(1,1)+GARCH(1,1)
  armaOrder <- c(1,1) # ARMA order
  garchOrder <- c(1,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Arma11Garch11FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj]) # fit
  
  #model2: ARMA(1,0)+GARCH(1,1)
  armaOrder <- c(1,0) # ARMA order
  garchOrder <- c(1,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Arma10Garch11FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj]) # fit
  
  #model3: : GARCH(1,1)
  armaOrder <- c(0,0) # ARMA order
  garchOrder <- c(1,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Garch11FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj]) # fit
  
  #model4: ARMA(1,1)+GARCH(2,1)
  armaOrder <- c(1,1) # ARMA order
  garchOrder <- c(2,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Arma11Garch21FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj]) # fit
  
  #model5: ARMA(1,2)+GARCH(2,1)
  armaOrder <- c(1,2) # ARMA order
  garchOrder <- c(2,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Arma12Garch21FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj]) # fit
  
  #model6: ARMA(1,1)
  armaOrder <- c(1,1) # ARMA order
  garchOrder <- c(0,0) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  Arma11FitStud <- ugarchfit(spec, data = y$NlogReturn[(nsdj-block+1):nsdj], solver = "nlminb") # fit
  
  # MODEL SELECTION ON THE BASIS OF STANDARD STATISTICAL CRITERIONS
  infocriteria(Arma11Garch11FitStud)
  infocriteria(Arma10Garch11FitStud)
  infocriteria(Arma11Garch21FitStud)
  infocriteria(Arma12Garch21FitStud)
  infocriteria(Garch11FitStud)       #<- best one (lowest information loss)
  infocriteria(Arma11FitStud)
  
  Errors <- as.numeric(residuals(Garch11FitStud, standardize=TRUE))
  
  plot(Date[(nsdj-block+1):nsdj],
       Errors, 
       type="l", xlab="Time", ylab="", cex.lab=1.3, cex.axis=1.3, main=paste(coin, " - RESIDUALS"), ylim=c(-4,7))
  
  acf(Errors)
  pacf(Errors)
  
  # PRELIMINARY
  ny <- nsdj # <- ???????????
  bigBlock <- ceiling((log(ny))^2)
  smallBlock <- ceiling(2*log(ny))
  k <- seq(10, 300, by=1)                            
  nk <- length(k)
  pn <- 1/block
  alpha_n <- 1 - pn
  
  # INITIALIZE RESULTS
  Tail_index_Hill <- array(0, c(nk,3))  # <- tail indexes
  Tail_index_Mom <- rep(0, nk)
  Tail_index_ML <- array(0, c(nk,3))
  
  Extreme_tau1 <- rep(0, nk)            # <- extreme tau
  
  Exp_LAWS_1 <- array(0, c(nk,6))  # <- extreme expectiles
  Exp_QB_1 <- array(0, c(nk,6))
  Exp_LAWS_2 <- array(0, c(nk,6))
  Exp_QB_2 <- array(0, c(nk,6))
  
  Qua_WEIS <- array(0, c(nk,8))    # <- extreme quantile
  
  
  # estimation method for the extreme expectile (eg QB-based)
  estMethod <- "QB"
  
  ############### TAIL, EXPECTILES AND VaR FOR DIFFERENT K's  ############### 
  
  for(i in 1:nk){
    
    ############### TAIL INDEX  ###############
    
    # estimates the tail index by the Hill estimator
    temp <- HTailIndex(Errors,
                       k[i], 
                       TRUE, 
                       bigBlock=bigBlock,
                       smallBlock=smallBlock)
    Tail_index_Hill[i, ] <- c(temp$CIgamHat[1], temp$gammaHat, temp$CIgamHat[2])
    
    # estimates the tail index by the ML estimator
    temp <- MLTailIndex(Errors,
                        k[i],
                        TRUE,
                        bigBlock=bigBlock,
                        smallBlock=smallBlock)
    Tail_index_ML[i,] <- c(temp$CIgamHat[1], temp$gammaHat, temp$CIgamHat[2])
    
    # estimates the tail index by the moment estimator
    Tail_index_Mom[i] <- MomTailIndex(Errors,
                                      k[i])
    
    ###############    EXTREME TAU    ###############
    
    # estimates the extreme level tau_n' given alpha_n
    Extreme_tau1[i] <- estExtLevel(alpha_n,
                                   Errors,
                                   k=k[i])$tauHat
    
    ############### EXTREME EXPECTILE  ###############
    
    ############################  TAU'=f (alpha_n) 
    # 1.1] LAWS-D 
    #      TAU'=alpha_n
    #      D-95%-CI 
    temp <- predExpectiles(Errors, 
                           NULL, 
                           alpha_n,
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock,
                           k=k[i])
    
    Exp_LAWS_1[i,1] <- temp$CIExpct[1]
    Exp_LAWS_1[i,2] <- temp$ExpctHat
    Exp_LAWS_1[i,3] <- temp$CIExpct[2]
    Exp_LAWS_1[i,4] <- temp$VarExtHat
    
    # 1.2] LAWS-D 
    #      TAU'=alpha_n 
    #      IID-95%-CI 
    temp <- predExpectiles(Errors,
                           NULL, 
                           alpha_n, 
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock,
                           varType="asym-Ind-Log", 
                           k=k[i])
    Exp_LAWS_1[i,5] <- temp$CIExpct[1]
    Exp_LAWS_1[i,6] <- temp$CIExpct[2]
    
    # 1.3] QB-D       
    #      TAU'=alpha_n
    #      D-95%-CI 
    temp <- predExpectiles(Errors,
                           NULL, 
                           alpha_n,
                           estMethod, 
                           var=TRUE,
                           bigBlock=bigBlock, 
                           smallBlock=smallBlock,
                           k=k[i])
    Exp_QB_1[i,1] <- temp$CIExpct[1]
    Exp_QB_1[i,2] <- temp$ExpctHat
    Exp_QB_1[i,3] <- temp$CIExpct[2]
    Exp_QB_1[i,4] <- temp$VarExtHat
    
    # 1.4] QB-D 
    #      TAU'=alpha_n
    #      IID-95%-CI 
    temp <- predExpectiles(Errors, 
                           NULL,
                           alpha_n,
                           estMethod, 
                           var=TRUE,
                           bigBlock=bigBlock, smallBlock=smallBlock, 
                           varType="asym-Ind-Log",
                           k=k[i])
    Exp_QB_1[i,5] <- temp$CIExpct[1]
    Exp_QB_1[i,6] <- temp$CIExpct[2]
    
    ############################
    # 2.1] LAWS-D 
    #      TAU'= f(alpha_n)
    #      D-95%-CI 
    temp <- predExpectiles(Errors,
                           NULL,
                           NULL, 
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock, 
                           k=k[i], 
                           alpha_n=alpha_n)
    
    Exp_LAWS_2[i,1] <- temp$CIExpct[1]
    Exp_LAWS_2[i,2] <- temp$ExpctHat
    Exp_LAWS_2[i,3] <- temp$CIExpct[2]
    Exp_LAWS_2[i,4] <- temp$VarExtHat
    
    # 2.2] LAWS-D 
    #      TAU'= f(alpha_n) 
    #      IID-95%-CI 
    temp <- predExpectiles(Errors,
                           NULL,
                           NULL,
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock,
                           varType="asym-Ind-Log",
                           k=k[i],
                           alpha_n=alpha_n)
    Exp_LAWS_2[i,5] <- temp$CIExpct[1]
    Exp_LAWS_2[i,6] <- temp$CIExpct[2]
    
    # 2.3] QB-D 
    #      TAU'= f(alpha_n)
    #      D-95%-CI 
    temp <- predExpectiles(Errors,
                           NULL,
                           NULL, 
                           estMethod,
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock,
                           k=k[i],
                           alpha_n=alpha_n)
    Exp_QB_2[i,1] <- temp$CIExpct[1]
    Exp_QB_2[i,2] <- temp$ExpctHat
    Exp_QB_2[i,3] <- temp$CIExpct[2]
    Exp_QB_2[i,4] <- temp$VarExtHat
    
    # 2.4] QB-D 
    #      TAU'= f(alpha_n) 
    #      IID-95%-CI 
    temp <- predExpectiles(Errors, 
                           NULL, 
                           NULL, 
                           estMethod, 
                           var=TRUE,
                           bigBlock=bigBlock,
                           smallBlock=smallBlock,
                           varType="asym-Ind-Log",
                           k=k[i], 
                           alpha_n=alpha_n)
    Exp_QB_2[i,5] <- temp$CIExpct[1]
    Exp_QB_2[i,6] <- temp$CIExpct[2]
    
    
    
    ############### EXTREME QUANTILE  ###############
    
    #  WEISS with D-95%-CI 
    temp <- extQuantile(Errors, 
                        NULL,
                        alpha_n,
                        var=TRUE,
                        bigBlock=bigBlock,
                        smallBlock=smallBlock,
                        k=k[i])
    
    Qua_WEIS[i,1] <- temp$CIExtQ[1]
    Qua_WEIS[i,2] <- temp$ExtQHat
    Qua_WEIS[i,3] <- temp$CIExtQ[2]
    Qua_WEIS[i,4] <- temp$VarExQHat
    
    #  WEISS with IID-95%-CI 
    temp <- extQuantile(Errors, 
                        NULL,
                        alpha_n,
                        var=TRUE,
                        k=k[i],
                        bigBlock=bigBlock,
                        smallBlock=smallBlock,
                        varType="asym-Ind")
    Qua_WEIS[i,5] <- temp$CIExtQ[1]
    Qua_WEIS[i,6] <- temp$CIExtQ[2]
    
    #  WEISS with D-ADJ-95%-CI 
    temp <- extQuantile(Errors, 
                        NULL,
                        alpha_n,
                        var=TRUE,
                        k=k[i],
                        bigBlock=bigBlock,
                        smallBlock=smallBlock,
                        varType="asym-Dep-Adj")
    Qua_WEIS[i,7] <- temp$CIExtQ[1]
    Qua_WEIS[i,8] <- temp$CIExtQ[2]
  }

  
  ##  DYNAMIC ESTIMATION
  
  ```{r}
  # MODEL SELECTION 
  armaOrder <- c(0,0) # ARMA order
  garchOrder <- c(1,1) # GARCH order
  varModel <- list(model = "sGARCH", garchOrder = garchOrder)
  spec <- ugarchspec(varModel, mean.model = list(armaOrder = armaOrder),
                     distribution.model = "std") # without fixed parameters here
  
  
  # INITIALIZE RESULTS
  
  Expc_LAWS_DYN <- NULL
  Expc_QB_DYN <- NULL
  Expc_LAWS_tau_hat_DYN <- NULL
  Expc_QB_tau_hat_DYN <- NULL
  Quan_DYN <- NULL
  Quan_BIAS_DYN <- NULL
  start_time <- Sys.time()
  
  
  ############### DYNAMIC  ############### 
  
  for(i in 1:(nsdj-block+1)){
    
    # fit ARMA(0,0)-GARCH(1,1) 
    GarchDataBlock <- try(ugarchfit(spec, data = y$NlogReturn[i:(i+block-1)], solver = "hybrid"), silent=TRUE)
    
    
    if(class(GarchDataBlock)[1]=="try-error") {
      Expc_LAWS_DYN <- c(Expc_LAWS_DYN, NaN)
      Expc_QB_DYN <- c(Expc_QB_DYN, NaN)
      Expc_LAWS_tau_hat_DYN <- c(Expc_LAWS_tau_hat_DYN, NaN)
      Expc_QB_tau_hat_DYN <- c(Expc_QB_tau_hat_DYN, NaN)
      Quan_DYN <- c(Quan_DYN, NaN)
      Quan_BIAS_DYN <- c(Quan_BIAS_DYN, NaN)
    }
    else{
      ################################################# 
      # EXTRACTS (standardized) RESIDUALS
      ResidualsST <- as.numeric(residuals(GarchDataBlock, standardize=TRUE)) 
      ResidualsNST <- as.numeric(residuals(GarchDataBlock, standardize=FALSE))
      
      # EXTRACTS CONDITIONAL-S/D
      Sigma <- as.numeric(sigma(GarchDataBlock))
      
      # PREDITCT CONDITIONAL MEAN at time t+1 (mu_{t+1})
      mu_np1 <- as.numeric(coef(GarchDataBlock)[names(coef(GarchDataBlock))=="mu"])
      
      # PREDITCT CONDITIONALS/D time n+1 -> 
      sigma_np1 <- sqrt(as.numeric(coef(GarchDataBlock)[names(coef(GarchDataBlock))=="omega"] +
                                     coef(GarchDataBlock)[names(coef(GarchDataBlock))=="alpha1"] * (ResidualsNST[block])^2 +
                                     coef(GarchDataBlock)[names(coef(GarchDataBlock))=="beta1"] * Sigma[block]^2))
      
      #################################################
      # EXPECTILE-LAWS with TAU'=alpha_n 
      temp <- predExpectiles(ResidualsST,
                             NULL, 
                             alpha_n, 
                             var=TRUE,
                             varType="asym-Ind-Log", k=k_fix)$ExpctHat
      Expc_LAWS_DYN <- c(Expc_LAWS_DYN, mu_np1 + sigma_np1 * temp)
      
      
      
      # EXPECTILE-QB  with TAU'=alpha_n 
      temp <- predExpectiles(ResidualsST,
                             NULL, 
                             alpha_n,
                             method=estMethod,
                             var=TRUE, varType="asym-Ind-Log", k=k_fix)$ExpctHat
      Expc_QB_DYN <- c(Expc_QB_DYN, mu_np1 + sigma_np1 * temp)
      
      
      #################################################
      
      # EXPECTILE-LAWS with TAU'=f (alpha_n) 
      temp <- predExpectiles(ResidualsST,
                             NULL, 
                             NULL,
                             var=TRUE, varType="asym-Ind-Log", k=k_fix, alpha_n=alpha_n)$ExpctHat
      Expc_LAWS_tau_hat_DYN <- c(Expc_LAWS_tau_hat_DYN, mu_np1 + sigma_np1 * temp)
      
      # EXPECTILE-QB with TAU'=f (alpha_n) 
      temp <- predExpectiles(ResidualsST,
                             NULL,
                             NULL,
                             method=estMethod, 
                             var=TRUE, varType="asym-Ind-Log", k=k_fix, alpha_n=alpha_n)$ExpctHat
      Expc_QB_tau_hat_DYN <- c(Expc_QB_tau_hat_DYN, mu_np1 + sigma_np1 * temp)
      #Exp_QB_Dyn_2
      
      #################################################
      
      # QUANTILE of LEVEL 1-tau_n' for RESIDUALS
      temp <- extQuantile(ResidualsST,
                          NULL, 
                          alpha_n,
                          var=FALSE, bias=TRUE, k=k_fix)$ExtQHat
      Quan_DYN <- c(Quan_DYN, mu_np1 + sigma_np1 * temp)
      #Qua_B_Dyn
      
      temp <- extQuantile(ResidualsST,
                          NULL,
                          alpha_n,
                          var=FALSE, bias=FALSE, k=k_fix)$ExtQHat
      Quan_BIAS_DYN <- c(Quan_BIAS_DYN, mu_np1 + sigma_np1 * temp)
      #Qua_NB_Dyn
      
      if(i%%125==0){
        cat(paste("iter:", round(i/(nsdj-block)*100, digits = 2), round(as.numeric(Sys.time()-start_time, units = "secs")),"seconds \n"))
        start_time <- Sys.time()}
      
    }
  }
  ###### EXCEEDANCES OVER RISK_MEASURES
  
  #INITIALIZE
  N_Expc_LAWS_DYN <- 0
  Date_Expc_LAWS_DYN <- as.Date(NULL, origin = "1970-01-01")
  NLR_Expc_LAWS_DYN <- NULL
  
  N_Expc_QB_DYN <- 0
  Date_Expc_QB_DYN <- as.Date(NULL, origin = "1970-01-01")
  NLR_Expc_QB_DYN <- NULL
  
  N_Expc_LAWS_tau_hat_DYN <- 0
  Date_Expc_LAWS_tau_hat_DYN <- as.Date(NULL, origin = "1970-01-01")
  NLR_Expc_LAWS_tau_hat_DYN <- NULL
  
  N_Expc_QB_tau_hat_DYN <- 0
  Date_Expc_QB_tau_hat_DYN <- as.Date(NULL, origin = "1970-01-01")
  NLR_Expc_QB_tau_hat_DYN <- NULL
  
  N_Quan_DYN <- 0
  Date_Quan_DYN <-  as.Date(NULL, origin = "1970-01-01")
  NLR_Quan_DYN <- NULL
  
  N_Quan_BIAS_DYN <- 0
  Date_Quan_BIAS_DYN <- as.Date(NULL, origin = "1970-01-01")
  NLR_Quan_BIAS_DYN <- NULL
  
  
  
  for(i in (block+1):(nsdj)){
    if(is.na(Expc_LAWS_DYN[i-block])==TRUE){
      cat("")
    }
    else{
      if(y$NlogReturn[i]>Expc_LAWS_DYN[i-block]){            
        Date_Expc_LAWS_DYN <- c(Date_Expc_LAWS_DYN, y$Date[i])
        NLR_Expc_LAWS_DYN <- c(NLR_Expc_LAWS_DYN, y$NlogReturn[i])
        N_Expc_LAWS_DYN <- N_Expc_LAWS_DYN + 1
      }
      
      if(y$NlogReturn[i]>Expc_QB_DYN[i-block]){
        Date_Expc_QB_DYN <- c(Date_Expc_QB_DYN, y$Date[i])
        NLR_Expc_QB_DYN <- c(NLR_Expc_QB_DYN, y$NlogReturn[i])
        N_Expc_QB_DYN <- N_Expc_QB_DYN + 1
      }
      
      if(y$NlogReturn[i]>Expc_LAWS_tau_hat_DYN[i-block]){            
        Date_Expc_LAWS_tau_hat_DYN <- c(Date_Expc_LAWS_tau_hat_DYN, y$Date[i])
        NLR_Expc_LAWS_tau_hat_DYN <- c(NLR_Expc_LAWS_tau_hat_DYN, y$NlogReturn[i])
        N_Expc_LAWS_tau_hat_DYN <- N_Expc_LAWS_tau_hat_DYN + 1
      }
      
      if(y$NlogReturn[i]>Expc_QB_tau_hat_DYN[i-block]){
        Date_Expc_QB_tau_hat_DYN <- c(Date_Expc_QB_tau_hat_DYN, y$Date[i])
        NLR_Expc_QB_tau_hat_DYN <- c(NLR_Expc_QB_tau_hat_DYN, y$NlogReturn[i])
        N_Expc_QB_tau_hat_DYN <- N_Expc_QB_tau_hat_DYN + 1
      }
      
      if(y$NlogReturn[i]>Quan_DYN[i-block]){            
        Date_Quan_DYN <- c(Date_Quan_DYN, y$Date[i])
        NLR_Quan_DYN <- c(NLR_Quan_DYN, y$NlogReturn[i])
        N_Quan_DYN <- N_Quan_DYN + 1
      }
      
      if(y$NlogReturn[i]>Quan_BIAS_DYN[i-block]){            
        Date_Quan_BIAS_DYN <- c(Date_Quan_BIAS_DYN, y$Date[i])
        NLR_Quan_BIAS_DYN <- c(NLR_Quan_BIAS_DYN, y$NlogReturn[i])
        N_Quan_BIAS_DYN <- N_Quan_BIAS_DYN + 1
      }
    }
  }
  
  
  #initialize
  setwd("/Users/andreateruzzi/Desktop/tesi/plots/Dynamic_risk_measures")
  pdf(paste0(coin,"_Dynamic_Errors_k",k_fix,".pdf" ), width=6, height=6)
  par(mai=c(.5,.5,.2,.1), mgp=c(1.6,.6,0))
  
  
  ################ ERRORS ################ 
  # 1.1] Errors 
  plot(Date[(nsdj-block+1):nsdj], 
       Errors, 
       type="l", xlab="Time", ylab="",cex.lab=1.3, cex.axis=1.3, main=paste(coin, " - RESIDUALS"),
       ylim=c(-4,7))
  legend("top", col=1, lwd=1, bty="n", cex=1.3, lty=1,
         legend="Residuals of Negative Log-Returns")
  
  # 1.2] ACF errors
  acf(Errors, 
      ylab="", main=paste("Autocorrelation Function Residuals -",coin), cex.lab=1.3, cex.axis=1.3)
  text(15,1, "Autocorrelation Function Residuals", cex=1.5)
  
  # 1.3] ACF errors^"
  acf(Errors^2,
      ylab="", main="Autocorrelation Function Residuals", cex.lab=1.3, cex.axis=1.3)
  text(15,1, "Autocorrelation Function Squared Residuals", cex=1.5)
  
  # 1.4] QQ plot 
  qqnorm(Errors,
         main="QQ-PLOT", cex.lab=1.3, cex.axis=1.3, pch=20, xlab="Normal Quantile")
  
  # 1.5] QQ line 
  qqline(Errors,
         lwd=2, col="blue")
  
  
  ################  TAIL ESTIMATORS PLOTS ################ 
  
  #2.1] TAIL INDEX MB, ML and HILL
  
  plot(k, 
       Tail_index_Mom, 
       type="l", ylim=c(-.1,1.5), lwd=2, ylab="",
       xlab="k", cex.lab=1.3, cex.axis=1.3, main=paste(coin, "RESIDUALS - TAIL ESTIMATORS"))
  lines(k,
        Tail_index_ML[,2],
        lwd=3, col=3, lty=3)
  lines(k,
        Tail_index_Hill[,2], 
        lwd=2, col=4, lty=2)
  abline(h=0, col=2, lty=4, lwd=2)
  legend("top", col=c(1,3,4), lwd=c(2,3,2), bty="n", cex=1.3, lty=c(1,3,2),
         legend=c(expression(hat(gamma)[n]~"MB"), expression(hat(gamma)[n]~"ML"),
                  expression(hat(gamma)[n]~"HILL")))
  
  #2.2] TAIL INDEX HILL and 95% CI
  plot(k, Tail_index_Hill[,2], type="l", ylim=c(-.1,1.5), lwd=2, ylab="",
       xlab="k", cex.lab=1.3, cex.axis=1.3, main=paste(coin, "RESIDUALS - HILL TAIL ESTIMATOR"))
  lines(k, Tail_index_Hill[,1], lty=2, lwd=2, col="blue")
  lines(k, Tail_index_Hill[,3], lty=2, lwd=2, col="blue")
  abline(h=0, col=2, lty=4, lwd=2)
  legend("top", lty=c(1,2), col=c("black", "blue"), lwd=c(2,2),
         bty="n", cex=1.3, legend=c(expression(hat(gamma)[n]~"HILL"),
                                    expression(paste(bold("95%-CI")))))
  
  
  
  ################ TAU'=alpha_n ################ 
  
  #3.1]  EXPECTILE-LAWS with 95%-CI-D and 95%-CI-IID 
  plot(k, 
       Exp_LAWS_1[,2], 
       type="l", ylim=c(0, 30), lwd=2, ylab="",
       xlab="k", cex.lab=1.3, cex.axis=1.3, col="limegreen", main=paste(coin,"RESIDUALS - TAU'=ALPHA"))
  
  lines(k,
        Exp_LAWS_1[,1], 
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_LAWS_1[,3], 
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_LAWS_1[,5],
        lwd=4, lty=3, col="gray60")
  lines(k,
        Exp_LAWS_1[,6], 
        lwd=4, lty=3, col="gray60")
  legend("topleft", lty=c(1,2,4,3), col=c("limegreen", "blue", "gray60"), lwd=c(2,2,4),
         legend=c(eval(bquote(expression(widetilde(xi^"*")[.(alpha_n)]))),
                  expression(paste("LAWS-D-", bold("95%-CI"))),
                  expression(paste("LAWS-IID-", bold("95%-CI")))),
         bty="n", cex=1.3)
  
  
  # 3.2]  EXPECTILE-QB with 95%-CI-D and 95%-CI-IID 
  plot(k, 
       Exp_QB_1[,2],
       type="l", ylim=c(0, 30), lwd=2, ylab="",
       xlab="k", cex.lab=1.3, cex.axis=1.3, col="limegreen", main=paste(coin,"RESIDUALS - TAU'=ALPHA"))
  lines(k,
        Exp_QB_1[,1],
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_QB_1[,3], lwd=2, lty=4, col="blue")
  lines(k,
        Exp_QB_1[,5], lwd=4, lty=3, col="gray60")
  lines(k, 
        Exp_QB_1[,6], lwd=4, lty=3, col="gray60")
  legend("topleft", lty=c(1,2,4,3), col=c("limegreen", "blue", "gray60"), lwd=c(2,2,4),
         legend=c(eval(bquote(expression(widehat(xi^"*")[.(alpha_n)]))),
                  expression(paste("QB-D-", bold("95%-CI"))),
                  expression(paste("QB-IID-", bold("95%-CI")))),
         bty="n", cex=1.3)
  
  
  # 3.3]  QUANTILE with 95%-CI-D, 95%-CI-IID and 95%-CI-D-ADJ
  plot(k, 
       Qua_WEIS[,2],
       type="l", ylim=c(0, 30), lwd=2, ylab="",
       xlab="k", cex.lab=1.3, cex.axis=1.3, col="limegreen", main=paste(coin,"RESIDUALS - QUANTILES TAU'=ALPHA"))
  lines(k,
        Qua_WEIS[,1],
        lwd=2, lty=4, col="blue")
  lines(k,
        Qua_WEIS[,3], 
        lwd=2, lty=4, col="blue")
  lines(k,
        Qua_WEIS[,5],
        lwd=4, lty=3, col="gray60")
  lines(k,
        Qua_WEIS[,6],
        lwd=4, lty=3, col="gray60")
  lines(k,
        Qua_WEIS[,7],
        lwd=2, lty=5, col="deepskyblue")
  lines(k,
        Qua_WEIS[,8], 
        lwd=2, lty=5, col="deepskyblue")
  legend("topleft", lty=c(1,4,3,5), col=c("limegreen", "blue", "gray60", "deepskyblue"), lwd=c(2,2,4,2),
         legend=c(eval(bquote(expression(widehat(q^"*")[.(alpha_n)]))),
                  eval(bquote(expression(widehat(q^"*")[.(alpha_n)]~bold("-D-95%-CI")))),
                  eval(bquote(expression(widehat(q^"*")[.(alpha_n)]~bold("-IID-95%-CI")))),
                  expression(paste("WEISS-D-ADJ-", bold("95%-CI")))),
         bty="n", cex=1.3)
  
  
  # 3.4]  DYNAMIC: TAU'= alpha_n
  # dynamicrisk measures
  axmax=min(max(Expc_LAWS_DYN[-1])+0.2,1.5)
  axmin=min(y$NlogReturn[(block+1):(nsdj+1)], na.rm = TRUE)
  plot(y$Date[(block+1):(nsdj+1)],  
       Expc_LAWS_DYN,                      # EXP-LAWS 
       type="l", ylim=c(axmin, axmax), lwd=1, ylab="", xlab="Time", cex.lab=1.3, cex.axis=1.3, col="limegreen", 
       main=paste(coin,"RESIDUALS - DYNAMIC PREDICTION TAU'=ALPHA"))
  lines(y$Date[(block+1):(nsdj+1)], 
        Expc_QB_DYN,                       # EXP-QB
        col="darkgray", lwd=1, lty=2)
  lines(y$Date[(block+1):(nsdj+1)],
        Quan_BIAS_DYN,                     # QUANTILE_BIAS
        col="blueviolet", lwd=1, lty=3)
  
  
  #static risk measures
  temp <- predExpectiles(y$NlogReturn[-1],
                         NULL, 
                         alpha_n, 
                         var=TRUE,
                         varType="asym-Ind-Log", k=k_fix)$ExpctHat
  
  abline(h=temp, col="red", lty=4, lwd=2)
  
  temp <- predExpectiles(y$NlogReturn[-1],
                         NULL, 
                         alpha_n, 
                         var=TRUE,
                         method=estMethod,
                         varType="asym-Ind-Log", k=k_fix)$ExpctHat
  
  abline(h=temp, col="blue", lty=4, lwd=2)
  
  temp <- extQuantile(y$NlogReturn[-1],
                      NULL, 
                      alpha_n,
                      var=FALSE, bias=TRUE, k=k_fix)$ExtQHat
  
  abline(h=temp, col="yellow", lty=4, lwd=2)
  
  
  #nlogReturns
  lines(y$Date[(block+1):(nsdj+1)], 
        y$NlogReturn[(block+1):(nsdj+1)],  # NLR
        col="black", lwd=1)
  
  #exceedances
  points(Date_Expc_LAWS_DYN, 
         NLR_Expc_LAWS_DYN,                # NLR 
         pch=8, col="deepskyblue", lwd=2)
  points(Date_Expc_QB_DYN,
         NLR_Expc_QB_DYN,                  # NLR 
         pch=4, col="darkorange", lwd=2)
  points(Date_Quan_BIAS_DYN,
         NLR_Quan_BIAS_DYN,                # NLR 
         pch=3, col="darkgreen", lwd=2)
  
  legend("topleft", lty=c(1,2,3,4,4,4), lwd=c(1,1,1,2,2,2),
         col=c("limegreen", "darkgray", "blueviolet","red", "blue", "yellow"),
         legend=c(eval(bquote(expression(widetilde(xi^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(widehat(xi^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(widehat(q^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(widetilde(xi^"*")[.(alpha_n)]))),
                  eval(bquote(expression(widehat(xi^"*")[.(alpha_n)]))),
                  eval(bquote(expression(widehat(q^"*")[.(alpha_n)])))),
         bty="n", cex=1.3)
  legend("topright", pch=c(8,4,3), lty=c(NA,NA,NA), lwd=c(2,2,2),
         col=c("deepskyblue", "darkorange", "darkgreen"),
         legend=c(eval(bquote(expression(italic(Y)[i]~">"~widetilde(xi^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(italic(Y)[i]~">"~widehat(xi^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(italic(Y)[i]~">"~widehat(q^"*")[.(alpha_n)](italic(Y)[n+1]~"|"~italic(F)[n]))))),
         bty="n", cex=1.3)
  
  
  ########################  TAU'= f(alpha_n) ######################## 
  
  # 4.1]  STATIC: LAWS with TAU'= f(alpha_n) with 95%-CI
  plot(k, 
       Exp_LAWS_2[,2],
       type="l", ylim=c(0, 30), lwd=2, ylab="", xlab="k", cex.lab=1.3, cex.axis=1.3, col="limegreen",
       main=paste(coin, " RESIDUALS - TAU'=TAU'(ALPHA)"))
  lines(k,
        Exp_LAWS_2[,1],
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_LAWS_2[,3],
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_LAWS_2[,5],
        lwd=4, lty=3, col="gray60")
  lines(k,
        Exp_LAWS_2[,6],
        lwd=4, lty=3, col="gray60")
  legend("topleft", lty=c(1,2,4,3), 
         col=c("limegreen", "blue", "gray60"), lwd=c(2,2,4),
         legend=c(eval(bquote(expression(widetilde(xi^"*")[widehat(tau~"'")[n](.(alpha_n))]))),
                  expression(paste("LAWS-D-", bold("95%-CI"))),
                  expression(paste("LAWS-IID-", bold("95%-CI")))),
         bty="n", cex=1.3)
  
  # 4.2]  STATIC: QB with TAU'= f(alpha_n) with 95%-CI
  plot(k,
       Exp_QB_2[,2],
       type="l", ylim=c(0, 30), lwd=2, ylab="", xlab="k", cex.lab=1.3, cex.axis=1.3, col="limegreen",
       main=paste(coin, "RESIDUALS - TAU'=TAU'(ALPHA)")
  )
  lines(k,
        Exp_QB_2[,1],
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_QB_2[,3],
        lwd=2, lty=4, col="blue")
  lines(k,
        Exp_QB_2[,5],
        lwd=4, lty=3, col="gray60")
  lines(k,
        Exp_QB_2[,6],
        lwd=4, lty=3, col="gray60")
  legend("topleft", lty=c(1,2,4,3), col=c("limegreen", "blue", "gray60"), lwd=c(2,2,4),
         legend=c(eval(bquote(expression(widehat(xi^"*")[widehat(tau~"'")[n](.(alpha_n))]))),
                  expression(paste("QB-D-", bold("95%-CI"))),
                  expression(paste("QB-IID-", bold("95%-CI")))),
         bty="n", cex=1.3)
  
  # 4.3] DYNAMIC: TAU'= f(alpha_n)
  
  # dyanmic risk measures
  axmax=min(max(Expc_LAWS_tau_hat_DYN[-1])+0.2,1.5)
  axmin=min(y$NlogReturn[(block+1):(nsdj+1)], na.rm = TRUE)
  plot(y$Date[(block+1):(nsdj+1)],
       Expc_LAWS_tau_hat_DYN,                #EXP-LAWS
       type="l", ylim=c(axmin,axmax), lwd=1, ylab="", xlab="Time", cex.lab=1.3, cex.axis=1.3, col="limegreen",
       main=paste(coin,"RESIDUALS - DYNAMIC PREDICTION TAU'= f(ALPHA)"))
  lines(y$Date[(block+1):(nsdj+1)],
        Expc_QB_tau_hat_DYN,                 #EXP-QB
        col="blueviolet", lwd=1)
  
  #static risk measures
  temp <- predExpectiles(y$NlogReturn[-1],
                         NULL, 
                         NULL,
                         var=TRUE, varType="asym-Ind-Log", 
                         k=k_fix, alpha_n=alpha_n)$ExpctHat
  
  abline(h=temp, col="red", lty=4, lwd=2)
  
  temp <- predExpectiles(y$NlogReturn[-1],
                         NULL, 
                         NULL,
                         method=estMethod, 
                         var=TRUE, varType="asym-Ind-Log", 
                         k=k_fix, alpha_n=alpha_n)$ExpctHat
  
  abline(h=temp, col="blue", lty=4, lwd=2)
  
  #nlogreturn
  lines(y$Date[(block+1):(nsdj+1)],
        y$NlogReturn[(block+1):(nsdj+1)],
        col="black", lwd=1)
  
  #exceedances
  points(Date_Expc_LAWS_tau_hat_DYN,
         NLR_Expc_LAWS_tau_hat_DYN,
         pch=3, col="deepskyblue", lwd=2)
  points(Date_Expc_QB_tau_hat_DYN,
         NLR_Expc_QB_tau_hat_DYN,
         pch=4, col="darkorange", lwd=2)
  
  
  legend("topleft", lty=c(1,1,4,4), lwd=c(1,1,2,2), col=c("limegreen", "blueviolet","red", "blue"),
         legend=c(eval(bquote(expression(widetilde(xi^"*")[widehat(tau~"'")[n](.(alpha_n))](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(widehat(xi^"*")[widehat(tau~"'")[n](.(alpha_n))](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(widetilde(xi^"*")[widehat(tau~"'")[n](.(alpha_n))]))),
                  eval(bquote(expression(widehat(xi^"*")[widehat(tau~"'")[n](.(alpha_n))])))),
         bty="n", cex=1.3)
  legend("topright",, pch=c(3,4), lty=c(NA,NA), lwd=c(2,2),
         col=c("deepskyblue", "darkorange"),
         legend=c(eval(bquote(expression(italic(Y)[i]~">"~widetilde(xi^"*")[widehat(tau~"'")[n](.(alpha_n))](italic(Y)[n+1]~"|"~italic(F)[n])))),
                  eval(bquote(expression(italic(Y)[i]~">"~widehat(xi^"*")[widehat(tau~"'")[n](.(alpha_n))](italic(Y)[n+1]~"|"~italic(F)[n]))))),
         bty="n", cex=1)
  
  
  ################
  dev.off()
  ################ END DYNAMIC PREDICTION DOW-JONES