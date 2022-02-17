library("ggplot2")
library("rugarch")
library("forecast")

##======================================================================
## Set up
##======================================================================
data = read.csv("real_1.csv", header = TRUE)
data$timestamp = as.POSIXct(data$timestamp, origin = "1970-01-01")
names(data) = c("timestamp", "Series", "anomaly_flag")

ggplot(data, aes(x = timestamp, y = Series, colour = anomaly_flag))+
  geom_line()+
  geom_point()
##----------------------------------------------------------------
## Ratio of Training set
##----------------------------------------------------------------
Train_ratio = 0.8

##----------------------------------------------------------------
## Identify Main Body Function
##----------------------------------------------------------------
arma_spec = autoarfima(data$Series[1:(round(nrow(data)*Train_ratio))], 
                       criterion = "BIC", 
                       method = "full", 
                       include.mean = TRUE)
arma_spec

ar1_Series = as.numeric(arma_spec$fit@fit$coef['ar1'])
ma1_Series = as.numeric(arma_spec$fit@fit$coef['ma1'])
mu_Series = as.numeric(arma_spec$fit@fit$coef['mu'])

N = nrow(data)


arma_garch_spec = ugarchspec(variance.model = list(model = "sGARCH",
                                                   garchOrder = c(1, 1),
                                                   submodel = NULL,
                                                   external.regressors = NULL,
                                                   variance.targeting = TRUE),
                             mean.model = list(armaOrder = c(1, 1),
                                               include.mean = TRUE,
                                               archm = FALSE,
                                               archpow = 1,
                                               arfima = FALSE,
                                               external.regressors = NULL,
                                               archex = FALSE),
                             fixed.pars = list(ar1 = ar1_Series,
                                               ma1 = ma1_Series,
                                               mu  = mu_Series),
                             
                             distribution.model = "norm")

arma_garch_model = ugarchfit(arma_garch_spec, data = data$Series[1:round(Train_ratio*N)], solver="hybrid")
as.numeric(arma_garch_model@fit$coef["omega"]/(1 - arma_garch_model@fit$coef["alpha1"] - arma_garch_model@fit$coef["beta1"]))

##======================================================================
## Anomaly Detection
##======================================================================

##----------------------------------------------------------------
## Function of Computing the Z Statistic of Each Normal Observation in Train Set
##----------------------------------------------------------------
Z_stat = function(M, N, Train_ratio, data)
{
  Var1 = c()
  Var2 = c()
  
  Z_stat = c()
  
  P_value = c()
  
  Stationary = c()
  
  for (i in (M+1):round(Train_ratio*N))
  {
    if(sum(data$anomaly_flag[(i-M):i]) == 0)
    {
      arma_garch_spec = ugarchspec(variance.model = list(model = "sGARCH",
                                                         garchOrder = c(1, 1),
                                                         submodel = NULL,
                                                         external.regressors = NULL,
                                                         variance.targeting = TRUE),
                                   mean.model = list(armaOrder = c(1, 1),
                                                     include.mean = TRUE,
                                                     archm = FALSE,
                                                     archpow = 1,
                                                     arfima = FALSE,
                                                     external.regressors = NULL,
                                                     archex = FALSE),
                                   fixed.pars = list(ar1 = ar1_Series,
                                                     ma1 = ma1_Series,
                                                     mu  = mu_Series),
                                  
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
      
      P_value[i] = p_value
      
      Stationary[i] = as.numeric(coef(arma_garch_model1)["alpha1"]) + as.numeric(coef(arma_garch_model1)["beta1"])
      
    }else{
      Var1[i] = 0
      Var2[i] = 0
      
      Z_stat[i] = 0
      P_value[i] = 0
    }
  }
  Train_result = as.data.frame(cbind(Z_stat, P_value, Var1, Var2, Stationary))
  
  return(Train_result)
}

##----------------------------------------------------------------
## Function of Anomaly Detection in Test Set
##----------------------------------------------------------------
anomaly_detect = function(est, Train_ratio, M, N, data)
{
  Var1 = c()
  Var2 = c()
  
  Z_stat = c()
  P_value = c()
  
  anomaly_detect = c()
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
                                 mean.model = list(armaOrder = c(1, 1),
                                                   include.mean = TRUE,
                                                   archm = FALSE,
                                                   archpow = 1,
                                                   arfima = FALSE,
                                                   external.regressors = NULL,
                                                   archex = FALSE),
                                 fixed.pars = list(ar1 = ar1_Series,
                                                   ma1 = ma1_Series,
                                                   mu  = mu_Series),
                                 
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
      anomaly_detect[i] = 1                                         
      Raw_value[i] = test_data[i]
      
      Forecast_temp = ugarchforecast(arma_garch_model1,
                                     n.ahead = 1,
                                     data = test_data[(i-M):(i-1)])
      test_data[i] = as.numeric(Forecast_temp@forecast$seriesFor) 
      New_value[i] = test_data[i]
    }else{
      anomaly_detect[i] = 0
      Raw_value[i] = test_data[i]
      New_value[i] = 0
    }
    Var1[i] = var1
    Var2[i] = var2
    
    Z_stat[i] = z_stat
    
    P_value[i] = p_value
  }
  Test_result = as.data.frame(cbind(Z_stat, P_value, Var1, Var2, anomaly_detect, Raw_value, New_value))
  return(Test_result)
}
#----------------------------------------------------------------
## Function of Using Bagging Strategy
#----------------------------------------------------------------
Bagging = function(M_list, N, Train_ratio, data, Trim_num)
{
  M_num = length(M_list)
  Result_List = list()
  for (i in 1:M_num)
  {
    Train_Result = Z_stat(M = M_list[i], N, Train_ratio, data)
    
    
    est = list(z_stat_max = sort(Train_Result$Z_stat[(Train_Result$Z_stat != 0)], decreasing = TRUE)[1]*(1+Extend_Ratio),
               z_stat_min = sort(Train_Result$Z_stat[(Train_Result$Z_stat != 0)], decreasing = FALSE)[1]*(1-Extend_Ratio))
    ##----------------------------------------------------------------
    ## Store Training Results
    ##----------------------------------------------------------------
    temp = anomaly_detect(est, Train_ratio, M =  M_list[i], N, data)
    Test_result = temp[(is.na(temp$Z_stat) !=TRUE), ]
    
    hist(Test_result$Z_stat, main = paste("Histogram of m = ", M_list[[i]], sep = " "), breaks = 10)
    Result_List[[i]] = cbind(No. = which(Test_result$anomaly_detect == 1) + round(N*Train_ratio),
                             Test_result[which(Test_result$anomaly_detect == 1), ])
  }
  return(Result_List)
}

#----------------------------------------------------------------
## Function of Computing AUC
#----------------------------------------------------------------
AUC = function(Detect_Result, data, Train_ratio, N)
{
  Score = data.frame(Score = rep(0,times = N))
  Score[as.numeric(names(Detect_Result)),] = as.numeric(Detect_Result)
  
  data[,4] = Score
  
  Test_Result = data[((round(N*Train_ratio) + 1):N),]
  
  Pos = Test_Result[which(Test_Result$anomaly_flag == 1),]
  Neg = Test_Result[which(Test_Result$anomaly_flag == 0),]
  
  n_pos = nrow(Pos)
  n_neg = nrow(Neg)
  
  count = 0
  
  for(i in 1:n_pos)
  {
    for(j in 1:n_neg)
    {
      if(Pos$Score[i]>Neg$Score[j])
      {
        count = count + 1
      }else if(Pos$Score[i] == Neg$Score[j])
      {
        count = count + 0.5
      }else{
        count = count
      }
    }
  }
  
  AUC = count/(n_pos*n_neg)
  
  return(AUC)
}

##======================================================================
##Implement Anomaly Detection
##======================================================================
##----------------------------------------------------------------
## Bagging Strategy for anomaly Detection
##----------------------------------------------------------------
##----------------------------------------------------------------
## Different M for Bagging. Make Sure That the First M Observation Are Not Anomalies
##----------------------------------------------------------------
Extend_Ratio = 0.05

M_list = c(110, 120, 130, 140, 150, 160, 170, 180, 190)
Result_List = Bagging(M_list, N, Train_ratio, data, Trim_num)

No_record = c()
for (k in 1:length(M_list))
{
  no_record = Result_List[[k]]$No.
  No_record = c(No_record, no_record)
  
}
Freq_table = table(No_record)/length(M_list)
Freq_table[which((Freq_table>=0.5))]

Detect_Result = Freq_table[which((Freq_table>=0.5))]
AUC(Detect_Result, data, Train_ratio, N)

