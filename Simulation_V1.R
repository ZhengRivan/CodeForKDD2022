library("ggplot2")
library("rugarch")
library("forecast")
##======================================================================
## Simulation
##======================================================================
##----------------------------------------------------------------
## Set up
##----------------------------------------------------------------
##----------------------------------------------------------------
## True Parameters in ARMA(1,1)
##----------------------------------------------------------------
a = 0.7
b = 0.3

##----------------------------------------------------------------
## True Parameters in GARCH
##----------------------------------------------------------------
omega=0.9
alpha=0.4
beta=0.3

##----------------------------------------------------------------
## Start point, Total number of observation,Total number of anomalies
##----------------------------------------------------------------
X0 = 0
Pre_M = 24
N = 1000
Num_Anomaly = 20

##----------------------------------------------------------------
## Simulation Function for time series with anomalies
##----------------------------------------------------------------
TSO.sim = function(Pre_M, N, a, b, alpha, beta, omega, X0, Num_Anomaly)
{
  h = rep(NA,N+1)
  epsilon = rep(NA,N+1)
  
  ##----------------------------------------------------------------
  ## Observation without anomalies
  ##----------------------------------------------------------------
  X = rep(NA,N+1)
  
  ##----------------------------------------------------------------
  ## Observation with anomalies
  ##----------------------------------------------------------------
  Series = rep(NA,N+1)
  Add_Anomaly = rep(1,N+1)
  Anomaly_flag = rep(0,N+1)
  
  ##----------------------------------------------------------------
  ## Location of anomalies
  ##----------------------------------------------------------------
  set.seed(12345)
  Anomaly_loc = sort(sample((Pre_M+1):N, Num_Anomaly))
  Anomaly_flag[Anomaly_loc] = 1
  
  ##----------------------------------------------------------------
  ## Add omega to anomalies
  ##----------------------------------------------------------------
  Add_Anomaly[Anomaly_loc] = runif(Num_Anomaly, min = 4, max = 5)
  yita = rnorm(N+1, mean = 0, sd = 1)
  
  ##----------------------------------------------------------------
  ## Initialization
  ##----------------------------------------------------------------
  h[1]=omega/(alpha+beta)
  X[1] = X0
  epsilon[1] = rnorm(1,mean = 0, sd = sqrt(h[1]))
  
  ##----------------------------------------------------------------
  ## Generate a centered ARMA(1,1)-GARCH(1,1) process with additive anomalies
  ##----------------------------------------------------------------
  for(i in 2:(N+1)){
    
    h[i]=omega+alpha*(epsilon[i-1]^2)+beta*h[i-1]
    epsilon[i]=yita[i]*sqrt(h[i])
    X[i] = a*X[i-1] + epsilon[i] + b*epsilon[i-1]
    if(i %in% Anomaly_loc)
    {
      Series[i] = a*X[i-1] + sign(epsilon[i])*Add_Anomaly[i] + epsilon[i] - b*epsilon[i-1]
    }else
    {
      Series[i] = X[i]
    }
  }
  
  X = X[-1]
  Series = Series[-1]
  epsilon = epsilon[-1]
  Anomaly_flag = Anomaly_flag[-1]
  
  Sim = as.data.frame(cbind(X,Series, epsilon, Anomaly_flag))
  
  return(Sim)
}
##----------------------------------------------------------------
## Generate simulated a time series
##----------------------------------------------------------------
Data = TSO.sim(Pre_M,N,a,b,alpha,beta,omega,X0,Num_Anomaly)
##----------------------------------------------------------------
## Examine the Plausibility of simulation
##----------------------------------------------------------------
cbbPalette = c("red3", "blue4")
ggplot(data = Data, aes(x = 1:nrow(Data), y = Series, group = "1", color = factor(Anomaly_flag))) + 
    geom_line()+
    geom_point()+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'),
          legend.position=c(0,1))+
    theme_bw()+
   scale_colour_manual(label = c("Anomaly","Normal"),
                         breaks = c("1","0"),
                         name = "Label",
                         values = cbbPalette)

arma_spec_Series = auto.arima(Data$Series[1:800], max.p = 2, max.q = 2,max.D = 0, max.d = 0,
                                  ic = "aic", seasonal = FALSE)
arma_spec_X = auto.arima(Data$X[1:800], max.p = 2, max.q = 2,max.D = 0, max.d = 0,
                              ic = "aic", seasonal = FALSE)


ar1_Series = arma_spec_Series$coef["ar1"]
ma1_Series = arma_spec_Series$coef["ma1"]

ar1_X = arma_spec_X$coef["ar1"]
ma1_X = arma_spec_X$coef["ma1"]

myspec_r_anom=ugarchspec(variance.model = list(model = "sGARCH",
                                          garchOrder = c(1, 1),
                                          submodel = NULL,
                                          external.regressors = NULL,
                                          variance.targeting = TRUE),
                    mean.model = list(armaOrder = c(1,1),
                                      include.mean = FALSE,
                                      archm = FALSE,
                                      archpow = 1,
                                      arfima = FALSE,
                                      external.regressors = NULL,
                                      archex = FALSE),
                    fixed.pars = list(ar1 = ar1_Series,
                                      ma1 = ma1_Series),
                    distribution.model = "norm"
)
myspec_r_norm=ugarchspec(variance.model = list(model = "sGARCH",
                                               garchOrder = c(1, 1),
                                               submodel = NULL,
                                               external.regressors = NULL,
                                               variance.targeting = TRUE),
                         mean.model = list(armaOrder = c(1,1),
                                           include.mean = FALSE,
                                           archm = FALSE,
                                           archpow = 1,
                                           arfima = FALSE,
                                           external.regressors = NULL,
                                           archex = FALSE),
                         fixed.pars = list(ar1 = ar1_X,
                                           ma1 = ma1_X),
                         distribution.model = "norm"
)
modelfit_r_anom=ugarchfit(myspec_r_anom,data=Data$Series[1:800],solver="hybrid")
modelfit_r_norm=ugarchfit(myspec_r_norm,data=Data$X[1:800],solver="hybrid")
omega/(1 - alpha - beta)
as.numeric(coef(modelfit_r_anom)["omega"])/as.numeric(1 - coef(modelfit_r_anom)["alpha1"] - coef(modelfit_r_anom)["beta1"])
as.numeric(coef(modelfit_r_norm)["omega"])/as.numeric(1 - coef(modelfit_r_norm)["alpha1"] - coef(modelfit_r_norm)["beta1"])

##======================================================================
## Anomaly Detection
##======================================================================

##----------------------------------------------------------------
## Function of Computing the Z statistic of each normal observation in train set
##----------------------------------------------------------------
Z_stat = function(M, N, Train_ratio, data = Data)
{
  Var1 = c()
  Var2 = c()
  
  Z_stat = c()
  
  Anomaly_flag = c()
  
  for (i in (M+1):round(Train_ratio*N))
  {
    if(sum(data$Anomaly_flag[(i-M):(i-1)]) == 0)
    {
      arma_garch_spec = ugarchspec(variance.model = list(model = "sGARCH",
                                                         garchOrder = c(1, 1),
                                                         submodel = NULL,
                                                         external.regressors = NULL,
                                                         variance.targeting = TRUE),
                                   mean.model = list(armaOrder = c(1,1),
                                                     include.mean = FALSE,
                                                     archm = FALSE,
                                                     archpow = 1,
                                                     arfima = FALSE,
                                                     external.regressors = NULL,
                                                     archex = FALSE),
                                   fixed.pars = list(ar1 = ar1_Series,
                                                     ma1 = ma1_Series),
                                   distribution.model = "norm")
      
      arma_garch_model1 = ugarchfit(arma_garch_spec, data = data$Series[(i-M):(i-1)], solver="hybrid")
      arma_garch_model2 = ugarchfit(arma_garch_spec, data = data$Series[(i-M):i], solver="hybrid")
      
      var1 = as.numeric(coef(arma_garch_model1)["omega"])/(1 - as.numeric(coef(arma_garch_model1)["alpha1"])
                                                           - as.numeric(coef(arma_garch_model1)["beta1"]))
      var2 = as.numeric(coef(arma_garch_model2)["omega"])/(1 - as.numeric(coef(arma_garch_model2)["alpha1"])
                                                           - as.numeric(coef(arma_garch_model2)["beta1"]))
      a = length(residuals(arma_garch_model1)) - 1
      b = length(residuals(arma_garch_model2)) - 1
      
      z_stat = var1/var2
      p_value = 2*min(pf(z_stat, a, b), pf(z_stat, a, b, lower.tail = FALSE))
      
      Var1[i] = var1
      Var2[i] = var2
      
      Z_stat[i] = z_stat
      
      Anomaly_flag[i] = data$Anomaly_flag[i]
      
    }else{
      Var1[i] = 0
      Var2[i] = 0
      
      Z_stat[i] = 0
      Anomaly_flag[i] = data$Anomaly_flag[i]
    }
  }
  Train_result = as.data.frame(cbind(Z_stat, Anomaly_flag, Var1, Var2))
  
  return(Train_result)
}

##----------------------------------------------------------------
## Function of Anomaly detection in test set
##----------------------------------------------------------------
Anomaly_detect = function(est, Train_ratio, M, N, data = Data)
{
  Var1 = c()
  Var2 = c()
  
  Z_stat = c()
  P_value = c()
  
  Anomaly_detect = c()
  Raw_value = c()
  New_value = c()
  
  test_data = data$Series
  
  for (i in (round(Train_ratio*N) + 1):N)
  {
    arma_garch_spec = ugarchspec(variance.model = list(model = "sGARCH",
                                                       garchOrder = c(1, 1),
                                                       submodel = NULL,
                                                       external.regressors = NULL,
                                                       variance.targeting = TRUE),
                                 mean.model = list(armaOrder = c(1,1),
                                                   include.mean = FALSE,
                                                   archm = FALSE,
                                                   archpow = 1,
                                                   arfima = FALSE,
                                                   external.regressors = NULL,
                                                   archex = FALSE),
                                 fixed.pars = list(ar1 = ar1_Series,
                                                   ma1 = ma1_Series),
                                 distribution.model = "norm")
    arma_garch_model1 = ugarchfit(arma_garch_spec, data = test_data[(i-M):(i-1)], solver="hybrid")
    arma_garch_model2 = ugarchfit(arma_garch_spec, data = test_data[(i-M):i], solver="hybrid")
    
    var1 = as.numeric(coef(arma_garch_model1)["omega"])/(1 - as.numeric(coef(arma_garch_model1)["alpha1"])
                                                         - as.numeric(coef(arma_garch_model1)["beta1"]))
    var2 = as.numeric(coef(arma_garch_model2)["omega"])/(1 - as.numeric(coef(arma_garch_model2)["alpha1"])
                                                         - as.numeric(coef(arma_garch_model2)["beta1"]))
    a = length(residuals(arma_garch_model1)) - 1
    b = length(residuals(arma_garch_model2)) - 1
    
    z_stat = var1/var2
    
    p_value = 2*min(pf(z_stat, a, b), pf(z_stat, a, b, lower.tail = FALSE))
    
    if((z_stat < est$z_stat_min) | (z_stat > est$z_stat_max))
    {
      Anomaly_detect[i] = 1                                         
      Raw_value[i] = test_data[i]
      
      ##----------------------------------------------------------------
      ## Correction
      ##---------------------------------------------------------------
      
      Forecast_temp = ugarchforecast(arma_garch_model1,
                                     n.ahead = 1,
                                     data = test_data[(i-M):(i-1)])
      test_data[i] = as.numeric(Forecast_temp@forecast$seriesFor) 
      New_value[i] = test_data[i]
    }else{
      Anomaly_detect[i] = 0
      Raw_value[i] = test_data[i]
      New_value[i] = 0
    }
    Var1[i] = var1
    Var2[i] = var2
    
    Z_stat[i] = z_stat
    
    P_value[i] = p_value
  }
  Test_result = as.data.frame(cbind(Z_stat, P_value, Var1, Var2, Anomaly_detect, Raw_value, New_value))
  return(Test_result)
}
#----------------------------------------------------------------
## Function of Using Bagging Strategy
#----------------------------------------------------------------
Bagging = function(M_list, N, Train_ratio, data = Data, Trim_num)
{
  M_num = length(M_list)
  Result_List = list()
  for (i in 1:M_num)
  {
    Train_Result = Z_stat(M = M_list[i], N, Train_ratio, data = Data)
    ##----------------------------------------------------------------
    ## Extension for Reducing False Positive
    ##---------------------------------------------------------------
    Hist_Data = na.omit(Train_Result[(Train_Result$Z_stat != 0),])
    
    p = ggplot(data = Hist_Data,
           aes(x = Z_stat, fill = factor(Anomaly_flag), color = factor(Anomaly_flag))) + 
           geom_histogram(bins = 10)+
           theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
           theme_bw()+
           scale_colour_manual(label = c("Anomaly","Normal"),
                          breaks = c("1","0"),
                          name = "Label",
                          values = c("black","black"))+
          scale_fill_discrete(name = "Label",
                              label = c("Anomaly","Normal"),
                              breaks = c("1","0"))+
          guides(color=FALSE)
    print(p)
    
    est = list(z_stat_max = sort(Train_Result$Z_stat[(Train_Result$Z_stat != 0) & (Train_Result$Anomaly_flag == 0)], decreasing = TRUE)[1]*(1+Extend_Ratio),
               z_stat_min = sort(Train_Result$Z_stat[(Train_Result$Z_stat != 0) & (Train_Result$Anomaly_flag == 0)], decreasing = FALSE)[1]*(1-Extend_Ratio))
    ##----------------------------------------------------------------
    ## Store Training Results
    ##----------------------------------------------------------------
    temp = Anomaly_detect(est, Train_ratio, M =  M_list[i], N, data = Data)
    Test_result = temp[(is.na(temp$Z_stat) !=TRUE), ]
    
    hist(Test_result$Z_stat, main = paste("Histogram of m = ", M_list[[i]], "in Testing Set", sep = " "), breaks = 10)
    Result_List[[i]] = cbind(No. = which(Test_result$Anomaly_detect == 1) + N*Train_ratio,
                             Test_result[which(Test_result$Anomaly_detect == 1), ])
  }
  return(Result_List)
}
##======================================================================
##Implement Anomaly detection
##======================================================================
##----------------------------------------------------------------
## Ratio of Training set
##----------------------------------------------------------------
Train_ratio = 0.8
##----------------------------------------------------------------
## Bagging Strategy for Anomaly Detection
##----------------------------------------------------------------
##----------------------------------------------------------------
## Different M for Bagging, make sure that the first M observation are not anomalies
##----------------------------------------------------------------
Extend_Ratio = 0.02

M_list = c(110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125)
Result_List = Bagging(M_list, N, Train_ratio, data = Data, Trim_num)

No_record = c()
for (k in 1:length(M_list))
{
  no_record = Result_List[[k]]$No.
  No_record = c(No_record, no_record)
  
}
Freq_table = table(No_record)/length(M_list)
Freq_table[which((Freq_table>=0.5))]

Data_Corr = Result_List[[2]]$New_value




