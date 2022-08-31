#Reproducable R-codes for "...." 2022
#Melanie Trouse & Vitaly Ganusov

data1 #dataframe with time and lengths
#required libraries
library(minpack.lm)
library(reshape2)
library(tidyr)
library(nlme)

#################################################
#######    Exponential Model Fits   #############
#################################################

ecg.fits<-NULL
etg.fits<-NULL

for(i in 2:length(data1)) {                   #cells start at column 2
  sub<-na.omit(data1[,c(1,i)])                #subset one cell 
  names(sub)<-c("hour","length")
  Td<-max(sub[,1])                                    #time of division
  
  #exponential constant growth rate model fit
  fitecg<-lm(log(length) ~ hour, data = sub)
  Lb<-exp(coef(fitecg)[1])                     #length at birth(e^intercept)
  lambda<-coef(fitecg)[2]                      #growth rate(slope)
  S<-shapiro.test(resid(fitecg))[2]            #shapiro p-value for residuals
  DeltaL<-Lb*(exp(lambda*Td)-1)               #change in length
  ecg.one<-data.frame(names(data1)[i],Td,       #save parameters from fit
                      Lb, lambda, DeltaL, S, AIC(fitecg))
  ecg.fits<-rbind(ecg.fits, ecg.one)
  
  #exponenital time-dependent growth rate model fit
  hour2<-(sub$hour)^2                 #square times for time-dependent growth
  fitetg<-lm(log(length) ~ hour+hour2, data = sub)
  Lb<-exp(coef(fitetg)[1])
  lambda<-coef(fitetg)[2]
  Dlambda<-coef(fitetg)[3]                  #Change in growth rate
  S<-shapiro.test(resid(fitetg))[2]
  DeltaL<-Lb*(exp(lambda*Td+Dlambda*Td^2)-1)
  
  
  #ANOVA for cell
  aa<-anova(fitecg,fitetg)                 #compare the two nested models
  etg.one<-data.frame(names(data1)[i], Td, Lb, lambda, 
                      Dlambda, DeltaL, S,AIC(fitetg), aa[2,6])
  etg.fits<-rbind(expquad,dataq)
  
}


############################################
########   Linear Model Fits ###############
############################################

lcg.fits<-NULL
ltg.fits<-NULL

for(i in 2:length(data1)) {
  sub<-na.omit(data1[,c(1,i)])
  names(sub)<-c("hour","length")
  Td<-max(sub[,1])
  
  #linear constant growth rate model fit
  fitlcg<-nlsLM(log(length)~log(Lb*(1+hour*lambda)), data=sub, 
               start=list(Lb=1.5, lambda=0.2))
  S<-shapiro.test(resid(fitlcg))[2]
  Lb<-coef(fitlcg)[1]
  lambda<-coef(fitlcg)[2]
  DeltaL<-Lb*lambda*Td
  ecg.one<-data.frame(names(data1)[i], Td, Lb,lambda, DeltaL, S, AIC(fitlcg))
  lcg.fits<-rbind(lcg.fits, datal)
  
  
  #linear quadratic model fit
  hour2<-(sub$hour)^2
  fitltg<-nlsLM(log(length)~log(Lb*(1+hour*lambda + hour2*Dlambda)), data=sub, 
              start=list(Lb=1.5, lambda=0.2, Dlambda=0.01))
  S<-shapiro.test(resid(fitltg))[2]
  Lb<-coef(fitltg)[1]
  lambda<-coef(fitltg)[2]
  Dlambda<-coef(fitltg)[3]
  DeltaL<-Lb*(lambda*Td+Dlambda*Td^2)
  
  #ANOVA for cell
  aa<-anova(fitlcg,fitltg)
  ltg.one<-data.frame(names(data1)[i], Lb, lambda, 
                    Dlambda, DeltaL,S, AIC(fitltg), aa[2,6])
  ltg.fits<-rbind(ltg.fits, ltg.one)
}

##############################################
####### Mixed Effects Models #################
##############################################
#dataframe to long form

data1.long<-na.omit(pivot_longer(data1, cols=2:length(1), 
                         names_to = "hour", values_to  = "length"))

#make an hour^2 column
data1.long$hour2<-data1.long$hour^2

#group data
data1.new<-groupedData(length~hour+hour2 | Cell.ID, data=data1.long)


#add control 
nlmeControl(maxIter = 100, msVerbose = TRUE, minScale = .0001)
#ECG fit
ecg.me<-nlme(log(length)~log(Lb) + r*hour,
             data = data1.new,
             fixed = Lb + r~ 1,
             start=c(Lb= 1.5, r = 0.2))

#ETG fit
etg.me<-nlme(log(length)~log(Lb) + r*hour + Deltar*hour2, 
             data=data1.new,
             fixed=Lb + r +Deltar~1,
             start=c(Lb=1.5, r=0.2, Deltar=0.05))

#LCG
lcg.me<-nlme(log(length)~log(Lb*(1+ r*hour)),
             data= data1.new, 
             fixed= Lb + r ~1,
             start=c(Lb=1.2, r=0.1))

#LTG
ltg.me<-nlme(log(length)~log(Lb*(1+r*hour+Deltar*hour2)),
             data=data1.new,
             fixed= Lb +r +Deltar ~1,
             start=c(Lb=1.5, r=0.2, Deltar=0.02))