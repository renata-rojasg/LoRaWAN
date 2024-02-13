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

############################
## Fitting the algorithms ##
############################
MAE<-MAPE<-RMSE<-COR<-matrix(NA,8,5)
colnames(MAE)<-colnames(MAPE)<-colnames(RMSE)<-
  colnames(COR)<-c("ARIMA-ANN","ARIMA-RF","ARIMA-SVM","ARIMA-LSTM","ARIMA-COV**")
rownames(MAPE)<-rownames(RMSE)<-
  rownames(COR)<-names(data)[6:13]

X<-cbind(temp_train,hum_train,bar_train,rain_train)
Xtest<-cbind(temp_test,hum_test,bar_test,rain_test)
Xchoosed<-X[,1]
Xchoosedt<-Xtest[,1]
m<-4 

for(i in 1:8){
  if(i==8){n=dim(datatrain)[1]-149}
  RSSIcomplete <- get(paste0("RSSI_0",i))#
  RSSI <- get(paste0("RSSI_0",i,"_train"))#
  # fitting the algorithms
  a01<-assign(paste0("arimax0",i), 
              auto.arima(RSSI,xreg = X,allowdrift=FALSE))
  a04<-Arima(RSSI,order=arimaorder(a01),xreg=Xchoosed)
  res04<-a04$residuals
  res04_data<- na.omit(data.frame(res=res04,
                                  lag1=lag.xts(res04,1),
                                  lag2=lag.xts(res04,2),
                                  lag3=lag.xts(res04,3),
                                  lag4=lag.xts(res04,4)
  ))[-c(1:m),]
  set.seed(10)
  set_random_seed(10)
  annt04<-neuralnet::neuralnet(res~.,data=res04_data)
  rf04<-randomForest::randomForest(res~.,data=res04_data)
  svmt04<-e1071::svm(res~.,data=res04_data)
  # forecasting
  RSSI_test <- get(paste0("RSSI_0",i,"_test"))
  new4<-assign(paste0("arima_cov2star0",i),
               Arima(RSSI_test,xreg=Xchoosedt,
                     model=a04)) #one-step-ahead
  res_test04<-c(a04$residuals[(length(RSSI)-m+1):length(RSSI)],
                new4$residuals)
  res_test04_data<-data.frame(res=res_test04,
                              lag1=lag.xts(res_test04,1),
                              lag2=lag.xts(res_test04,2),
                              lag3=lag.xts(res_test04,3),
                              lag4=lag.xts(res_test04,4)
  )[-c(1:m),]
  f_ann04<-predict(annt04,res_test04_data)+new4$fitted
  f_rf04<-predict(rf04,res_test04_data)+new4$fitted
  f_svm04<-predict(svmt04,res_test04_data)+new4$fitted
  
  f_complete<-Arima(RSSIcomplete,xreg=temp,
                    model=a04)
  res_test<-f_complete$residuals
  if(i==8){
    res_test<-res_test[-c(1:149)]
    TSLSTM<-ts.lstm(ts=res_test,
                    tsLag=4,xregLag = 0,
                    SplitRatio = 0.7825,
                    LSTMUnits=5, Epochs=20)
  }else{
    TSLSTM<-ts.lstm(ts=res_test,
                    tsLag=4,xregLag = 0,
                    LSTMUnits=5, Epochs=20)
  }
  f_lstm04<-(TSLSTM$TestPredictedValue)+new4$fitted[-374]
  
  MAPE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[5],
              forecast::accuracy(RSSI_test,f_rf04)[5],
              forecast::accuracy(RSSI_test,f_svm04)[5],
              forecast::accuracy(RSSI_test,f_lstm04)[5],
              forecast::accuracy(RSSI_test,new4$fitted)[5]
  )
  RMSE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[2],
              forecast::accuracy(RSSI_test,f_rf04)[2],
              forecast::accuracy(RSSI_test,f_svm04)[2],
              forecast::accuracy(RSSI_test,f_lstm04)[2],
              forecast::accuracy(RSSI_test,new4$fitted)[2]
  )
  COR[i,]<-c(cor(RSSI_test,f_ann04),
             cor(RSSI_test,f_rf04),
             cor(RSSI_test,f_svm04),
             cor(RSSI_test[-374],f_lstm04),
             cor(RSSI_test,new4$fitted)
  )
  MAE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[3],
             forecast::accuracy(RSSI_test,f_rf04)[3],
             forecast::accuracy(RSSI_test,f_svm04)[3],
             forecast::accuracy(RSSI_test,f_lstm04)[3],
             forecast::accuracy(RSSI_test,new4$fitted)[3]
  )
  
  assign(paste0("result0",i),
         t(data.frame(MAE=MAE[i,],MAPE=MAPE[i,],RMSE=RMSE[i,],COR=COR[i,]))
  )
}

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

# Calculating the percentage difference 
# of  with  ARIMAX_all respect to ARIMAX_Temp
MAE_all<-data.frame(values=c(MAE_AUM[,1],M_AUM[,1],
                             RMSE_AUM[,1],COR_AUM[,1]),
                    measure=c(rep("MAE",8),rep("MAPE",8),
                              rep("RMSE",8),rep("COR",8)),
                    model=rep(colnames(data[,6:13]),4)
)

MAE_all <- MAE_all %>% 
  mutate(plot_text=
           case_when(values<0 ~ .2,
                     values>0 ~ values+.2),
         
  )

plot_result1 <- ggplot(MAE_all,
                       aes(y=values,x=model,
                           fill=factor(measure,
                                       levels = c("MAE","MAPE","RMSE","COR")))
)+
  geom_bar(stat="identity",position = "dodge")+
  labs(fill="",y="Percentage differences",x="") +
  scale_fill_manual(values=c("grey85","grey70","grey45","grey30") )+
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold",size=15),
        plot.title = element_text(face="bold",size=15),
        legend.text = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold", color="black",
                                    size=15),
        axis.title.x = element_text(face="bold", color="black",
                                    size=15),
        axis.text.x = element_text(face="bold", color="black",
                                   size=15),
        axis.text.y = element_text(face="bold", color="black",
                                   size=15),
        panel.background = element_rect(fill = "white", colour = "white"))

pdif <- as.data.frame(result[,6:9]) %>%
  rownames_to_column(var = "Metric") %>%
  pivot_longer(cols = starts_with("ARIMA"), 
               names_to = "Model", values_to = "Value")

plot_result2 <- ggplot(pdif, aes(x=factor(Model,levels=colnames(result[,6:9])),y = Value)) +
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
plot_result1
plot_result2