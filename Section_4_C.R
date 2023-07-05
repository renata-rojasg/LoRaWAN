rm(list=ls())
gc()
library(tidyverse)
library(readr)
library(xts)
library(lmtest)
library(forecast)
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
##########################
## Descreptive analysis ##
##########################
summary1<-bind_rows(apply(data[,2:21],2,summary))[7]
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

####################
## Fitting ARIMAX ##
####################
MAE<-MAPE<-RMSE<-COR<-matrix(NA,8,4)
colnames(MAE)<-colnames(MAPE)<-colnames(RMSE)<-
  colnames(COR)<-c("ARIMA-ANN","ARIMA-RF*","ARIMA-SVM","ARIMA-COV**")
rownames(MAPE)<-rownames(RMSE)<-
  rownames(COR)<-names(data)[6:13]

X<-cbind(temp_train,hum_train,bar_train,rain_train)
Xtest<-cbind(temp_test,hum_test,bar_test,rain_test)
Xchoosed<-X[,1]
Xchoosedt<-Xtest[,1]
m<-4 

for(i in 1:8){
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
  
  MAPE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[5],
              forecast::accuracy(RSSI_test,f_rf04)[5],
              forecast::accuracy(RSSI_test,f_svm04)[5],
              forecast::accuracy(RSSI_test,new4$fitted)[5]
  )
  RMSE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[2],
              forecast::accuracy(RSSI_test,f_rf04)[2],
              forecast::accuracy(RSSI_test,f_svm04)[2],
              forecast::accuracy(RSSI_test,new4$fitted)[2]
  )
  COR[i,]<-c(cor(RSSI_test,f_ann04),
             cor(RSSI_test,f_rf04),
             cor(RSSI_test,f_svm04),
             cor(RSSI_test,new4$fitted)
  )
  MAE[i,]<-c(forecast::accuracy(RSSI_test,f_ann04)[3],
             forecast::accuracy(RSSI_test,f_rf04)[3],
             forecast::accuracy(RSSI_test,f_svm04)[3],
             forecast::accuracy(RSSI_test,new4$fitted)[3]
  )
  
  RSSI_test<-xts(RSSI_test, order.by=data$timestamp[(n+1):(dim(data)[1])])
  new1fit<-xts(as.numeric(f_ann04), order.by=data$timestamp[(n+1):(dim(data)[1])])
  new4fit<-xts(new4$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])
  
  limits<-c(min(RSSI_test[275:374],new1fit[275:374],new4fit[275:374])-1,
            max(RSSI_test[275:374],new1fit[275:374],new4fit[275:374])+1.5)
  if(i==2){limits=c(-94,-72)}
  if(i==5){limits=c(-83,-63)}
  
  assign(paste0("p",i), {
    plot(RSSI_test[275:374],main="",
         yaxis.right=FALSE, grid.col = "white",
         format.labels="%b-%Y", main.timespan = FALSE,
         lwd=0.5,cex.lab=1.3,cex.axis=1.3,
         ylim=limits)
    lines(new1fit[275:374],col=2,lty = 2,lwd = 1.5)
    lines(new4fit[275:374],col=4,lty = 5,lwd = 1.5)
    addLegend("topright",
              legend.names=c("Original data","ARIMA-ANN",expression(ARIMAX['Temp'])),
              col=c(1,2,4), cex=1.2,lty=c(1,2,4),
              lwd=c(1,1.5,1.5),
              ncol=2,
              bg="white")
  }
  )
  
  assign(paste0("result0",i),
         t(data.frame(MAE=MAE[i,],MAPE=MAPE[i,],RMSE=RMSE[i,],COR=COR[i,]))
  )
}

# Calculating the percentage difference 
# of  with  ARIMAX-ANN respect to ARIMAX_Temp
MAE_AUM<-(MAE[,4]-MAE[,1:3])/MAE[,4]*100
M_AUM<-(MAPE[,4]-MAPE[,1:3])/MAPE[,4]*100
RMSE_AUM<-(RMSE[,4]-RMSE[,1:3])/RMSE[,4]*100
COR_AUM<-(COR[,1:3]-COR[,4])/COR[,4]*100

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
                   COR=rank(result01[4,])==4),1,sum)
for(i in 2:8){
  r<-get(paste0("result0",i))
  r<-apply(cbind(apply(r[1:3,], 1, rank)==1,
                 COR=rank(r[4,])==4),1,sum)
  count<-abind::abind(count,r,along = 2)
}
count<-abind::abind(count,apply(count,1,sum),along = 2)
colnames(count)<-c(rownames(MAPE),"Overall")

print(t(count))

# fitted versus true values
p1; p2; p3
p4; p5; p6
p7; p8

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

ggplot(MAE_all %>% 
         filter(model %in% c("RSSI_01","RSSI_03","RSSI_08")), 
       aes(y=values,x=model,
           fill=factor(measure,
                       levels = c("MAE","MAPE","RMSE","COR")))
)+
  geom_bar(stat="identity",position = "dodge")+
  labs(fill="",y="Percentage differences",x="") +
  geom_text(aes(x=model,y=plot_text,
                label=round(values,3)),
            fontface="bold",
            angle=50,
            position = position_dodge(width = 1),
            color="black",
            size=5)+
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


ggplot(MAE_all,
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

