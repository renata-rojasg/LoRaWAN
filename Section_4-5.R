rm(list=ls())
gc()
library(tidyverse)
library(readr)
library(xts)
library(lmtest)
library(forecast)
library(xtable)
######################
## Data preparation ##
######################
data <- read_delim("combined_hourly_data.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                              origin="1970-01-01 00:00:00"))
# available at https://github.com/emanueleg/lora-rssi/blob/master/vineyard-2021_data/combined_hourly_data.csv
#
data<-data[1:1870,]
n<-round(dim(data)[1]*.8) 
dim(data)[1]-n
data$RSSI_03[is.na(data$RSSI_03)]<-
  mean(na.omit(data$RSSI_03))
data$RSSI_04[is.na(data$RSSI_04)]<-
  mean(na.omit(data$RSSI_04))
data$RSSI_05[is.na(data$RSSI_05)]<-
  mean(na.omit(data$RSSI_05))
data$RSSI_06[is.na(data$RSSI_06)]<-
  mean(na.omit(data$RSSI_06))
#########################
## Train and test sets ##
#########################
datatrain<-cbind(data[1:n,])  
colnames(datatrain)<-paste0(colnames(datatrain),"_train")

datatest<-cbind(data[(n+1):(dim(data)[1]),])  
colnames(datatest)<-paste0(colnames(datatest),"_test")

attach(datatrain)
attach(datatest)
attach(data)

#############################
## Missing data experiment ##
#############################
X<-cbind(temp_train,hum_train, bar_train,rain_train)
Xtest<-cbind(temp_test,hum_test,bar_test,rain_test)
Xchoosed<-temp_train
Xchoosedt<-temp_test
m<-4

for(i in 1:8){
  # training the algorithms
  RSSI_train <- get(paste0("RSSI_0",i,"_train"))#
  RSSI <- get(paste0("RSSI_0",i))#
  # fitting the algorithms
  a01<-assign(paste0("arimax0",i), 
              auto.arima(RSSI_train,xreg = X,allowdrift=FALSE))
  a04<-Arima(RSSI_train,order=arimaorder(a01),xreg=Xchoosed)
  RSSI_test <- get(paste0("RSSI_0",i,"_test"))#
  res04<-a04$residuals
  res04_data<- na.omit(data.frame(res=res04,
                                  lag1=lag.xts(res04,1),
                                  lag2=lag.xts(res04,2),
                                  lag3=lag.xts(res04,3),
                                  lag4=lag.xts(res04,4)
  ))[-c(1:m),]
  
  set.seed(10)
  annt04<-neuralnet::neuralnet(res~.,data=res04_data)
  
  hat1<-hat2<-hat3<-matrix(NA,length(RSSI_test),3)
  colnames(hat1)<-colnames(hat2)<-colnames(hat3)<-c("RMSE","MAE","MAPE")
  
  steps_values<-round(length(RSSI_test)*c(.01,.05,.1))
  steps<-length(steps_values)
  meanRMSE<-meanMAE<-meanMAPE<-
    rRMSE<-rMAE<-rMAPE<-matrix(NA,steps,3)
  
  for(j in 1:steps){
    n_ahead<-steps_values[j]
    hat1<-hat2<-hat3<-matrix(NA,length(RSSI_test)-n_ahead+1,3)
    colnames(hat1)<-colnames(hat2)<-colnames(hat3)<-c("RMSE","MAE","MAPE")
    for(k in 0:(length(RSSI_test)-n_ahead)){
      # updating the data
      RSSI_test1<-RSSI[1:(n+k)]
      n_new<-length(RSSI_test1)
      X1<-as.matrix(data[1:(n+k),2])
      Xtest1<-as.matrix(data[(n+k+1):(n+k+n_ahead),2])
      # fitted values
      new4_12<-Arima(RSSI_test1,model=a04,xreg = X1)
      # last n_ahead observations
      RSSI_prev<-RSSI_test1[(n_new-n_ahead+1):(n_new)] 
      # original values 
      RSSI_true<-RSSI[(n_new+1):(n_new+n_ahead)]
      # moving average
      RSSI_hat1<-rep(mean(RSSI_prev),n_ahead)
      # ARIMA-COV**
      RSSI_hat2<-predict(new4_12, 
                         n.ahead = n_ahead,
                         newxreg=Xtest1[,1])$pred 
      # ARIMA-ANN
      res_hat2<-RSSI_hat2-RSSI_true
      res_test04<-c(new4_12$residuals[(length(new4_12$residuals)-m+1):
                                        length(new4_12$residuals)],
                    res_hat2)
      res_test04_data<-data.frame(res=res_test04,
                                  lag1=lag.xts(res_test04,1),
                                  lag2=lag.xts(res_test04,2),
                                  lag3=lag.xts(res_test04,3),
                                  lag4=lag.xts(res_test04,4)
      )[-c(1:m),]
      RSSI_hat3<-predict(annt04,res_test04_data)+RSSI_hat2
      # accuracy measures
      hat1[k+1,]<-forecast::accuracy(RSSI_hat1, RSSI_true)[c(2,3,5)]
      hat2[k+1,]<-forecast::accuracy(RSSI_hat2, RSSI_true)[c(2,3,5)]
      hat3[k+1,]<-forecast::accuracy(RSSI_hat3, RSSI_true)[c(2,3,5)]
    }
    
    RMSE<-cbind(hat1[,1],hat2[,1],hat3[,1])
    MAE <-cbind(hat1[,2],hat2[,2],hat3[,2])
    MAPE<-cbind(hat1[,3],hat2[,3],hat3[,3])
    colnames(meanRMSE)<-colnames(RMSE)<-
      colnames(MAE)<-colnames(MAPE)<-
      c("MA-imputation", "ARIMA-COV** imputation", "ARIMA-ANN imputation")
    rownames(meanRMSE)<-c("1%","5%","10%")
    meanRMSE[j,]<-apply(RMSE,2,mean)
    meanMAE[j,]<-apply(MAE,2,mean)
    meanMAPE[j,]<-apply(MAPE,2,mean)
    rRMSE[j,]<-apply(apply(RMSE, 1, rank)==1,1,sum)
    rMAE[j,]<-apply(apply(MAE, 1, rank)==1,1,sum)
    rMAPE[j,]<-apply(apply(MAPE, 1, rank)==1,1,sum)
  }
  assign(paste0("RMSE_0",i),data.frame(
    model=c(rep("MA",3), 
            rep("ARIMA-COV**",3), 
            rep("ARIMA-ANN",3)),
    gap=rep(c("1%","5%","10%"),3),
    value=as.vector(meanRMSE)
  )
  )
  assign(paste0("p",i), ggplot(get(paste0("RMSE_0",i)), aes(x=gap,y=value,fill=model))+
           geom_bar(stat="identity",position = "dodge")+
           labs(fill="",x="Gap sizes", y="RMSE")  +
           # ylim(0,4)+
           geom_text(aes(x=gap,y=value/2,
                         label=round(value,3)), 
                     fontface="bold", 
                     angle=90,
                     position = position_dodge(width = 1),  
                     # color="black", 
                     size=6)+ 
           scale_x_discrete(limits=c("1%","5%","10%"))+
           scale_fill_manual(labels=c("ARIMA-ANN",
                                      expression(bold(ARIMAX['Temp'])), 
                                      "MA"),
                             values=c("grey90","grey65","grey45") )+
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
                 panel.background = element_rect(fill = "white", colour = "white")))
  
}

# Figure 7
p1 + guides(fill=guide_legend(ncol=2))

p2 + guides(fill=guide_legend(ncol=2))

p3 + guides(fill=guide_legend(ncol=2))

p4 + guides(fill=guide_legend(ncol=2))

p5 + guides(fill=guide_legend(ncol=2))

p6 + guides(fill=guide_legend(ncol=2))

p7 + guides(fill=guide_legend(ncol=2))

p8 + guides(fill=guide_legend(ncol=2))

