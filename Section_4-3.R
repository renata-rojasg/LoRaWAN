library(tidyverse)
library(readr)
library(xts)
library(lmtest)
library(forecast)
library("reticulate")
library(tensorflow)
library(keras)
library(TSLSTM)
######################
## Data preparation ##
######################
data <- read_delim("combined_hourly_data.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                              origin="1970-01-01 00:00:00"))
# available at https://github.com/emanueleg/lora-rssi/blob/master/vineyard-2021_data/combined_hourly_data.csv
data<-data[1:1870,]
n<-round(dim(data)[1]*.8)
dim(data)[1]-n

#########################
## Train and test sets ##
#########################
data$RSSI_03[is.na(data$RSSI_03)]<-
  mean(na.omit(data$RSSI_03))
data$RSSI_04[is.na(data$RSSI_04)]<-
  mean(na.omit(data$RSSI_04))
data$RSSI_05[is.na(data$RSSI_05)]<-
  mean(na.omit(data$RSSI_05))
data$RSSI_06[is.na(data$RSSI_06)]<-
  mean(na.omit(data$RSSI_06))
data$RSSI_08[335]<-
  mean(na.omit(data$RSSI_08))
datatrain<-cbind(data[1:n,])  
colnames(datatrain)<-paste0(colnames(datatrain),"_train")

datatest<-cbind(data[(n+1):(dim(data)[1]),])  
colnames(datatest)<-paste0(colnames(datatest),"_test")

attach(datatrain)
attach(datatest)
attach(data)
#############################
## Fitting ML-based Models ##
#############################
MAE<-MAPE<-RMSE<-COR<-matrix(NA,8,5)
colnames(MAE)<-colnames(MAPE)<-colnames(RMSE)<-
  colnames(COR)<-c("ANN","RF","SVM","LSTM","ARIMA-Temp")
rownames(MAPE)<-rownames(RMSE)<-
  rownames(COR)<-names(data)[6:13]
X<-cbind(temp_train,hum_train,bar_train,rain_train)
Xtest<-cbind(temp_test,hum_test,bar_test,rain_test)
Xchoosed<-X[,1]
Xchoosedt<-Xtest[,1]
m=8
normalize <- function(x) {
  return((x-min(x)) / (max(x)-min(x)))
}
denormalize <- function(x, original_min, original_max) {
  return (x * (original_max - original_min) + original_min)
}

for(i in 1:8){
  if(i==8){n=dim(datatrain)[1]-149}
  RSSI <- get(paste0("RSSI_0",i))#
  RSSIdata<-data.frame(RSSI,temp,hum,bar,rain,
                       lag1=lag.xts(RSSI,1),
                       lag2=lag.xts(RSSI,2),
                       lag3=lag.xts(RSSI,3),
                       lag4=lag.xts(RSSI,4),
                       lag5=lag.xts(RSSI,5),
                       lag6=lag.xts(RSSI,6),
                       lag7=lag.xts(RSSI,7),
                       lag8=lag.xts(RSSI,8)
                       
  )[-c(1:(sum(is.na(RSSI))+m)),]
  
  a01<-assign(paste0("arimax0",i), 
              auto.arima(RSSI[1:n],
                         xreg = X[c(sum(is.na(RSSI))+1:n),],
                         allowdrift=FALSE))
  a04<-Arima(RSSI[1:n],order=arimaorder(a01),
             xreg=Xchoosed[c(sum(is.na(RSSI))+1:n)])
  
  RSSI_norm <- as.data.frame(lapply(RSSIdata,normalize))
  RSSI_train <- RSSI_norm[1:(n-m), ]
  RSSI_test <- RSSI_norm[(n-m+1):(dim(RSSI_norm)[1]),]
  
  set.seed(10)
  set_random_seed(10)
  annt04<-neuralnet::neuralnet(RSSI ~ ., data = RSSI_train)
  rf04<-randomForest::randomForest(RSSI~.,data=RSSI_train)
  svmt04<-e1071::svm(RSSI~.,data=RSSI_train)
  
  if(i==8){RSSI<-RSSI[-c(1:149)]
  TSLSTM<-ts.lstm(ts=(RSSI),xreg = cbind(temp,hum,bar,rain)[-c(1:149),],
                  tsLag=4,xregLag = 4,
                  SplitRatio = 0.7825,
                  LSTMUnits=5, Epochs=20)
  }else{
    TSLSTM<-ts.lstm(ts=RSSI,xreg = cbind(temp,hum,bar,rain),
                    # SplitRatio = 0.7995,
                    tsLag=max(arimaorder(a01)),xregLag = 0,
                    LSTMUnits=5, Epochs=20)
  }
  
  # forecasting
  RSSItest<-get(paste0("RSSI_0",i,"_test"))
  new4<-assign(paste0("arima_cov2star0",i),
               Arima(RSSItest,xreg=Xchoosedt,
                     model=a04)) #one-step-ahead
  f_ann04p<-predict(annt04,RSSI_test)
  f_ann04<-denormalize(f_ann04p,min(RSSIdata$RSSI),max(RSSIdata$RSSI))
  f_rf04p<-predict(rf04,RSSI_test)
  f_rf04<-denormalize(f_rf04p,min(RSSIdata$RSSI),max(RSSIdata$RSSI))
  f_svm04p<-predict(svmt04,RSSI_test)
  f_svm04<-denormalize(f_svm04p,min(RSSIdata$RSSI),max(RSSIdata$RSSI))
  f_lstm<-(TSLSTM$TestPredictedValue) # the function normalizes internally and does not predict the last RSSItest value
  
  MAPE[i,]<-c(forecast::accuracy(RSSItest,f_ann04)[5],
              forecast::accuracy(RSSItest,f_rf04)[5],
              forecast::accuracy(RSSItest,f_svm04)[5],
              forecast::accuracy(RSSItest,f_lstm)[5],
              forecast::accuracy(RSSItest,new4$fitted)[5]
  )
  RMSE[i,]<-c(forecast::accuracy(RSSItest,f_ann04)[2],
              forecast::accuracy(RSSItest,f_rf04)[2],
              forecast::accuracy(RSSItest,f_svm04)[2],
              forecast::accuracy(RSSItest,f_lstm)[2],
              forecast::accuracy(RSSItest,new4$fitted)[2]
  )
  COR[i,]<-c(cor(RSSI_test$RSSI,f_ann04),
             cor(RSSI_test$RSSI,f_rf04),
             cor(RSSI_test$RSSI,f_svm04),
             cor(RSSI_test$RSSI[-374],f_lstm),
             cor(RSSI_test$RSSI,new4$fitted)
  )
  MAE[i,]<-c(forecast::accuracy(RSSItest,f_ann04)[3],
             forecast::accuracy(RSSItest,f_rf04)[3],
             forecast::accuracy(RSSItest,f_svm04)[3],
             forecast::accuracy(RSSItest,f_lstm)[3],
             forecast::accuracy(RSSItest,new4$fitted)[3]
  )
  assign(paste0("result0",i),
         t(data.frame(MAE=MAE[i,],MAPE=MAPE[i,],RMSE=RMSE[i,],COR=COR[i,]))
  )
}
print(MAPE)

# Calculating the percentage difference 
# of  with  ARIMAX-ANN respect to ARIMAX_Temp
MAE_AUM<-(MAE[,5]-MAE[,1:4])/MAE[,5]*100
M_AUM<-(MAPE[,5]-MAPE[,1:4])/MAPE[,5]*100
RMSE_AUM<-(RMSE[,5]-RMSE[,1:4])/RMSE[,5]*100
COR_AUM<-(COR[,1:4]-COR[,5])/COR[,5]*100

# organizing the table
result<- cbind(result01,rbind(
  MAE_AUM[1,],M_AUM[1,],RMSE_AUM[1,],COR_AUM[1,]
))
for(i in 2:8){
  r<-cbind(get(paste0("result0",i)),rbind(
    MAE_AUM[i,],M_AUM[i,],RMSE_AUM[i,],COR_AUM[i,]
  )
  )
  result<-abind::abind(result,r,along = 1)
}
print(round(result,3))

# Counting the times the models were the best option
count<-apply(cbind(apply(result01[1:3,], 1, rank)==1,
                   COR=rank(result01[4,])==5),1,sum)
for(i in 2:8){
  r<-get(paste0("result0",i))
  r<-apply(cbind(apply(r[1:3,], 1, rank)==1,
                 COR=rank(r[4,])==5),1,sum)
  count<-abind::abind(count,r,along = 2)
}
count<-abind::abind(count,apply(count,1,sum),along = 2)
colnames(count)<-c(rownames(MAPE),"Overall")

pdif <- as.data.frame(result[,6:9]) %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = -"Metric",  
               names_to = "Model", values_to = "Value")

result_plot<- ggplot(pdif, aes(x=factor(Model,levels=colnames(result[,6:9])),y = Value)) +
  geom_boxplot() +
  labs(title="",x="AI algorithm",
       y = "Percentage difference") +
  geom_hline(yintercept=0, linetype=2, 
             color = "grey0", size=.3)+
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold",size=15),
        plot.title = element_text(face="bold",size=15),
        legend.text = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold", color="black",
                                    size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    size=15),
        axis.text.x = element_text(face="bold", color="black",#angle = 50,
                                   size=13),
        axis.text.y = element_text(face="bold", color="black",
                                   size=15),
        panel.background = element_rect(fill = "white", colour = "black"))


print(t(count))
print(result)
result_plot
