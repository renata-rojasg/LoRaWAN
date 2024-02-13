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
data<-data[1:1870,]
n<-round(dim(data)[1]*.8) 
dim(data)[1]-n
#####################################
## Descreptive analysis - TABLE II ##
#####################################
summary1<-bind_rows(apply(data[,2:13],2,summary))[7]
data$RSSI_03[is.na(data$RSSI_03)]<-
  mean(na.omit(data$RSSI_03))
data$RSSI_04[is.na(data$RSSI_04)]<-
  mean(na.omit(data$RSSI_04))
data$RSSI_05[is.na(data$RSSI_05)]<-
  mean(na.omit(data$RSSI_05))
data$RSSI_06[is.na(data$RSSI_06)]<-
  mean(na.omit(data$RSSI_06))
summary2<-bind_rows(apply(data[,2:13],2,summary))
sums<-as.matrix(cbind(bind_rows(summary2)[-c(2,5,7)],
                      SD=apply(na.omit(data[,2:13]),2,sd),
                      `NA`=summary1
                      ))
sums[is.na(sums)] = 0
rownames(sums)<-colnames(data[,2:13])
sums # RSSI_08 starts 150 hours after


#######################################
## Pearson's correlation - TABLE III ##
#######################################
correlation<-data[,c(2:12)]
mcor<- psych::corr.test(correlation)
lower<-cbind(paste0(round(mcor$r[,1],3)," (",round(mcor$p[,1],3),")"),
             paste0(round(mcor$r[,2],3)," (",round(mcor$p[,2],3),")"),
             paste0(round(mcor$r[,3],3)," (",round(mcor$p[,3],3),")"),
             paste0(round(mcor$r[,4],3)," (",round(mcor$p[,4],3),")"))
lower[upper.tri(lower)]<-""
lower<-as.data.frame(lower)
rownames(lower)<-colnames(data[,c(2:12)])

correlation8<-na.omit(data[,c(2:5,13)])
mcor8<- psych::corr.test(correlation8)
results08<-t(as.matrix(paste0(round(mcor8$r[5,1:4],3)," (",
       round(mcor8$p[5,1:4],3),")"),ncol=4))
rownames(results08)<-c("RSSI_08")

table3<-rbind(lower,results08)
table3


#########################
## Time series Figures ##
#########################
RSSI_01 <- xts(data$RSSI_01, order.by=data$timestamp)
RSSI_02 <- xts(data$RSSI_02, order.by=data$timestamp)
RSSI_03 <- xts(data$RSSI_03, order.by=data$timestamp)
RSSI_04 <- xts(data$RSSI_04, order.by=data$timestamp)
RSSI_05 <- xts(data$RSSI_05, order.by=data$timestamp)
RSSI_06 <- xts(data$RSSI_06, order.by=data$timestamp)
RSSI_07 <- xts(data$RSSI_07, order.by=data$timestamp)
RSSI_08 <- xts(data$RSSI_08, order.by=data$timestamp)
{plot(RSSI_01,main="", yaxis.right=FALSE, grid.col = "white",
               format.labels="%b-%Y", main.timespan = FALSE,
               cex.axis=1.2,
               lwd=0.5,ylim=c(-95,-42),ylab="",cex.lab=1.2)
  par(cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2) 
  lines(RSSI_02,main="RSSI 02",col=2)
  lines(RSSI_03,main="RSSI 02",col=3)
  lines(RSSI_04,main="RSSI 02",col=4)
  addLegend("topright",
            legend.names=c("RSSI 01","RSSI 02","RSSI 03","RSSI 04"),
            col=1:4, cex=1.2,
            lwd=rep(.5,4),
            ncol=2,
            bg="white")
}
{
  plot(RSSI_05,main="", yaxis.right=FALSE, grid.col = "white",
       format.labels="%b-%Y", main.timespan = FALSE,
       cex.axis=1.2,
       lwd=0.5,ylim=c(-95,-42),ylab="",cex.lab=1.2)
  par(cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2) 
  lines(RSSI_06,main="",col=2)
  lines(RSSI_07,main="",col=3)
  lines(RSSI_08,main="",col=4)
  addLegend("topright",
            legend.names=c("RSSI 05","RSSI 06","RSSI 07","RSSI 08"),
            col=1:4, cex=1.2,
            lwd=rep(.5,4),
            ncol=2,
            bg="white")
}


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
order_arima<-matrix(NA,8,3)
MAE<-MAPE<-RMSE<-COR<-matrix(NA,8,4)
colnames(MAE)<-colnames(MAPE)<-colnames(RMSE)<-
  colnames(COR)<-c("ARIMA-COV","ARIMA-COV*","ARIMA-COV**","ARIMA")
rownames(MAPE)<-rownames(RMSE)<-rownames(order_arima)<-
  rownames(COR)<-names(data)[6:13]
Xsig<-values<-sinal<-matrix(0,8,4)
rownames(Xsig)<-colnames(data[,6:13])


coeffs<-matrix(0,8,6)
X<-cbind(temp_train,hum_train,bar_train,rain_train)
Xtest<-cbind(temp_test,hum_test,bar_test,rain_test)
Xchoosed<-X[,1]
Xchoosedt<-Xtest[,1]

for(i in 1:8){
  RSSI <- get(paste0("RSSI_0",i,"_train"))#
  # fitting the algorithms
  a01<-assign(paste0("arimax0",i), auto.arima(RSSI,xreg = X,allowdrift=FALSE))
  tcoef<-(coeftest(a01)<0.05)[(length(a01$coef)-dim(X)[2]+1):length(a01$coef),4]
  Xnew<-X[,tcoef]
  order_arima[i,]<-arimaorder(a01)
  Xsig[i,]<-c(c("T","RH","Bar","Rain")[tcoef],rep("",4-sum(tcoef)))
  sinal[i,]<-(coef(a01)<0)[(length(a01$coef)-dim(X)[2]+1):length(a01$coef)]
  values[i,]<-(coef(a01))[(length(a01$coef)-dim(X)[2]+1):length(a01$coef)]
  Xnewt<-Xtest[,tcoef]
  a02<-Arima(RSSI,arimaorder(a01),xreg=Xnew)
  a03<-Arima(RSSI,arimaorder(a01))
  a04<-Arima(RSSI,order=arimaorder(a01),xreg=Xchoosed)
  # forecasting
  RSSI_test <- get(paste0("RSSI_0",i,"_test"))
  new1<-assign(paste0("arima_cov0",i),
               Arima(RSSI_test,xreg = Xtest,model=a01)) #one-step-ahead
  new2<-assign(paste0("arima_covstar",i),
               Arima(RSSI_test,xreg = Xnewt,model=a02)) #one-step-ahead
  new3<-assign(paste0("arima_pred0",i),
               Arima(RSSI_test,model=a03)) #one-step-ahead
  new4<-assign(paste0("arima_cov2star0",i),
               Arima(RSSI_test,xreg=Xchoosedt,model=a04)) #one-step-ahead
   MAPE[i,]<-c(forecast::accuracy(RSSI_test,new1$fitted)[5],
              forecast::accuracy(RSSI_test,new2$fitted)[5],
              forecast::accuracy(RSSI_test,new4$fitted)[5],
              forecast::accuracy(RSSI_test,new3$fitted)[5]
  )
  RMSE[i,]<-c(forecast::accuracy(RSSI_test,new1$fitted)[2],
              forecast::accuracy(RSSI_test,new2$fitted)[2],
              forecast::accuracy(RSSI_test,new4$fitted)[2],
              forecast::accuracy(RSSI_test,new3$fitted)[2]
  )
  COR[i,]<-c(cor(RSSI_test,new1$fitted),
             cor(RSSI_test,new2$fitted),
             cor(RSSI_test,new4$fitted),
             cor(RSSI_test,new3$fitted)
  )
  MAE[i,]<-c(forecast::accuracy(RSSI_test,new1$fitted)[3],
             forecast::accuracy(RSSI_test,new2$fitted)[3],
             forecast::accuracy(RSSI_test,new4$fitted)[3],
             forecast::accuracy(RSSI_test,new3$fitted)[3]
  )

  RSSI_test<-xts(RSSI_test, order.by=data$timestamp[(n+1):(dim(data)[1])])
  new1fit<-xts(new1$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])
  new2fit<-xts(new2$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])
  new3fit<-xts(new3$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])
  new4fit<-xts(new4$fitted, order.by=data$timestamp[(n+1):(dim(data)[1])])
  
  assign(paste0("result0",i),
         t(data.frame(MAE=MAE[i,],MAPE=MAPE[i,],RMSE=RMSE[i,],COR=COR[i,]))
  )
}

print(cbind(order_arima,Xsig)) # TABLE IV

# Boxplots of the estimated beta-coefficients in the ARIMAX
colnames(values)<-c("T","RH","Bar","Rain")
values<-as.data.frame(values)
ggplot(stack(values), aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(title="",x="Weather parameter", 
       y = expression(paste(beta,"-coefficient estimates"))) +
  geom_hline(yintercept=0, linetype=2, 
             color = "grey0", size=.3)+
  theme(axis.title.y = element_text(color=1,size=15),
        axis.title.x = element_text(color=1,size=15),
        axis.text.x = element_text(color=1,size=15),
        axis.text.y = element_text(color=1,size=15),
        panel.background = element_rect(fill = "white", 
                                        colour = "black"))


# Calculating the percentage difference with respect to ARIMA
MAE_AUM<-(MAE[,4]-MAE[,1:3])/MAE[,4]
M_AUM<-(MAPE[,4]-MAPE[,1:3])/MAPE[,4]
RMSE_AUM<-(RMSE[,4]-RMSE[,1:3])/RMSE[,4]
COR_AUM<-(COR[,1:3]-COR[,4])/COR[,4]

# organizing the table
result<- cbind(result01,rbind(
  MAE_AUM[1,],M_AUM[1,],RMSE_AUM[1,],COR_AUM[1,]
  )*100
  )
for(i in 2:8){
  r<-cbind(get(paste0("result0",i)),rbind(
           MAE_AUM[i,],M_AUM[i,],RMSE_AUM[i,],COR_AUM[i,]
  )*100
  )
  result<-abind::abind(result,r,along = 1)
}
print(result,digits=3) # TABLE V

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

print(t(count)) # TABLE VI


# Calculating the percentage difference 
# of  with  ARIMAX_all respect to ARIMAX_Temp
MAE_all<-data.frame(values=c((MAE[,3]-MAE[,1])/MAE[,3],
                             (MAPE[,3]-MAPE[,1])/MAPE[,3],
                             (RMSE[,3]-RMSE[,1])/RMSE[,3],
                             (COR[,1]-COR[,3])/COR[,3])*100,
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
  ylim(-1.1,1.1)+
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



