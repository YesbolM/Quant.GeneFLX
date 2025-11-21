library(MCMCglmm)

#Import data
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

data <- read.csv("data.csv")
str(data)
data<-as.data.frame(data) # 6858 obs. of  13 variables
names(data)[1]<-"animal"
for(x in c(1:5,12:13))data[,x]<-as.factor(data[,x])
data$mLA.active<-ifelse(data$sex=="male",data$LA_active,NA)
meanActive.m<-mean(data$mLA.activ, na.rm=TRUE)
data$meanActive.m<-data$mLA.active/meanActive.m #mean standardization

data$mLA.passive<-ifelse(data$sex=="male",data$LA_passive,NA)
data$fLA.active<-ifelse(data$sex=="female",data$LA_active,NA)
meanActive.f<-mean(data$fLA.activ, na.rm=TRUE)
data$meanActive.f<-data$fLA.active/meanActive.f

data$fLA.passive<-ifelse(data$sex=="female",data$LA_passive,NA)
data$mWS<-ifelse(data$sex=="male",data$WS,NA)
meanWs.m<-mean(data$mWS,na.rm=TRUE)
data$meanWs.m<-data$mWS/meanWs.m

data$fWS<-ifelse(data$sex=="female",data$WS,NA)
meanWS.f<-mean(data$fWS,na.rm=TRUE)
data$meanWs.f<-data$fWS/meanWS.f

data$mDR<-ifelse(data$sex=="male",data$DR ,NA)
meanDR.m<-mean(data$mDR,na.rm=TRUE)
data$meanDr.m<-data$mDR/meanDR.m

data$fDR<-ifelse(data$sex=="female",data$DR,NA)
meanDR.f<-mean(data$fDR,na.rm=TRUE)
data$meanDr.f<-data$fDR/meanDR.f

colnames(data)[14] = "new.sec"

rep1<-data[data$replicat=="rep1",]
FLX1<-rep1[rep1$treatment=="FLX",]
CFM1<-rep1[rep1$treatment=="CFM",]
CWT1<-rep1[rep1$treatment=="CWT",]

rep2<-data[data$replicat=="rep2",]
FLX2<-rep2[rep2$treatment=="FLX",]
CFM2<-rep2[rep2$treatment=="CFM",]
CWT2<-rep2[rep2$treatment=="CWT",]

rep1FLX_p <- read.csv("rep1FLX.p.csv")
ped1<-rep1FLX_p
head(ped1)
str(ped1)
ped1<-as.data.frame(ped1)

rep1CFM_p <- read.csv("rep1CFM.p.csv")
ped2<-rep1CFM_p
head(ped2)
str(ped2)
ped2<-as.data.frame(ped2)

rep1CWT_p <- read.csv("rep1CWT.p.csv")
ped3<-rep1CWT_p
head(ped3)
str(ped3)
ped3<-as.data.frame(ped3)

rep2FLX_p <- read.csv("rep2FLX.p.csv")
ped4<-rep2FLX_p
head(ped4)
str(ped4)
ped1<-as.data.frame(ped4)

rep2CFM_p <- read.csv("rep2CFM.p.csv")
ped5<-rep2CFM_p
head(ped5)
str(ped5)
ped5<-as.data.frame(ped5)

rep2CWT_p <- read.csv("rep2CWT.p.csv")
ped6<-rep2CWT_p
head(ped6)
str(ped6)
ped6<-as.data.frame(ped6)


#Rep1 FLX

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1flx<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                   random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped1,data=FLX1,prior=p2,verbose=FALSE,
                   family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep1flx, file = "rep1flx")
rep1flx<-readRDS("rep1flx")

#Male Va wing size
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep1flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep1flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep1flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep1flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep1flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep1flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep1flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)

#Rep1 CFM

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1cfm<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                  random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped2,data=CFM1,prior=p2,verbose=FALSE,
                  family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep1cfm, file = "rep1cfm")
rep1cfm<-readRDS("rep1cfm")

#Male Va wing size
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep1cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep1cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep1cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep1cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep1cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep1cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep1cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)

#Rep1 CWT

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1cwt<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                  random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped3,data=CWT1,prior=p2,verbose=FALSE,
                  family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep1cwt, file = "rep1cwt")
rep1cwt<-readRDS("rep1cwt")

#Male Va wing size
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep1cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep1cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep1cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep1cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep1cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep1cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep1cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep1cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep1cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)

#rep2 FLX

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2flx<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                  random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped4,data=FLX2,prior=p2,verbose=FALSE,
                  family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep2flx, file = "rep2flx")
rep2flx<-readRDS("rep2flx")

#Male Va wing size
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep2flx$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep2flx$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep2flx$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2flx$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep2flx$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2flx$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep2flx$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep2flx$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep2flx$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2flx$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)

#Rep2 CFM

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2cfm<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                  random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped5,data=CFM2,prior=p2,verbose=FALSE,
                  family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep2cfm, file = "rep2cfm")
rep2cfm<-readRDS("rep2cfm")

#Male Va wing size
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep2cfm$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep2cfm$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep2cfm$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cfm$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep2cfm$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cfm$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep2cfm$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep2cfm$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep2cfm$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2cfm$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)

#Rep2 CWT

p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2cwt<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                  random=~us(trait):animal+us(trait):new.sec,rcov=~us(trait):units,pedigree=ped6,data=CWT2,prior=p2,verbose=FALSE,
                  family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
saveRDS(rep2cwt, file = "rep2cwt")
rep2cwt<-readRDS("rep2cwt")

#Male Va wing size
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"],0.95)
#Male Vr wing size
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"],0.95)
#Male heritability WS
h2.mWS <- rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]/
  (rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]+
     rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.new.sec"]+
     rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"],0.95)
#Female Vr WS
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"],0.95)
#Female heritability WS
h2.fWS <- rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]/
  (rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]+
     rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.new.sec"]+
     rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Male Va DR
posterior.mode(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"],0.95)
#Male Vr DR
posterior.mode(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"],0.95)
#Male heritability DR
h2.mDR <- rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]/
  (rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]+
     rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.new.sec"]+
     rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"],0.95)
#Female Vr DR
posterior.mode(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"],0.95)
#Female heritability DR
h2.fDR <- rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]/
  (rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]+
     rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.new.sec"]+
     rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Male Va LA
posterior.mode(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"],0.95)
#Male Vr LA
posterior.mode(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"],0.95)
#Male heritability LA
h2.mLA <- rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]/
  (rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]+
     rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.new.sec"]+
     rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"],0.95)
#Female Vr LA
posterior.mode(rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
HPDinterval(rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"],0.95)
#Female heritability LA
h2.fLA <- rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]/
  (rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]+
     rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.new.sec"]+
     rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Intersexual genetic covariance WS
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"],0.95)
#Intersexual genetic correlation WS
rmf.WS<-rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"],0.95)
#Intersexual genetic correlation DR
rmf.DR<-rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"],0.95)
#Intersexual genetic correlation LA
rmf.LA<-rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]*
           rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)

#Intra-sex cross-trait covariances
#mWS-mDR
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"],0.95)
#mWS-mLA
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"],0.95)
#mDR-mLA
posterior.mode(rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"],0.95)
#fWS-fDR
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"],0.95)
#fWS-fLA
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"],0.95)
#fDR-fLA
posterior.mode(rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"],0.95)
#Intra-sex cross-trait correlations
#mWS-mDR
rmWS.mDR<-rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rmWS.mDR)
HPDinterval(rmWS.mDR,0.95)
#mWS-mLA
rmWS.mLA<-rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmWS.mLA)
HPDinterval(rmWS.mLA,0.95)
#mDR-mLA
rmDR.mLA<-rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rmDR.mLA)
HPDinterval(rmDR.mLA,0.95)
#fWS-fDR
rfWS.fDR<-rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rfWS.fDR)
HPDinterval(rfWS.fDR,0.95)
#fWS-fLA
rfWS.fLA<-rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfWS.fLA)
HPDinterval(rfWS.fLA,0.95)
#fDR-fLA
rfDR.fLA<-rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rfDR.fLA)
HPDinterval(rfDR.fLA,0.95)

#Cross-sex cross-trait covariances
#mWS-fDR
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"],0.95)
#mWS-fLA
posterior.mode(rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"],0.95)
#mDR-fLA
posterior.mode(rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"],0.95)
#fWS-mDR
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"],0.95)
#fWS-mLA
posterior.mode(rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"],0.95)
#fDR-mLA
posterior.mode(rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"])
HPDinterval(rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"],0.95)
#Cross-sex cross-trait correlations
#mWS-fDR
rmWS.fDR<-rep2cwt$VCV[,"traitmeanWs.m:traitmeanDr.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"])
posterior.mode(rmWS.fDR)
HPDinterval(rmWS.fDR,0.95)
#mWS-fLA
rmWS.fLA<-rep2cwt$VCV[,"traitmeanWs.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]*
           rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmWS.fLA)
HPDinterval(rmWS.fLA,0.95)
#mDR-fLA
rmDR.fLA<-rep2cwt$VCV[,"traitmeanDr.m:traitmeanActive.f.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]*
           rep2cwt$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"])
posterior.mode(rmDR.fLA)
HPDinterval(rmDR.fLA,0.95)
#fWS-mDR
rfWS.mDR<-rep2cwt$VCV[,"traitmeanWs.f:traitmeanDr.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cwt$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"])
posterior.mode(rfWS.mDR)
HPDinterval(rfWS.mDR,0.95)
#fWS-mLA
rfWS.mLA<-rep2cwt$VCV[,"traitmeanWs.f:traitmeanActive.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]*
           rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfWS.mLA)
HPDinterval(rfWS.mLA,0.95)
#fDR-mLA
rfDR.mLA<-rep2cwt$VCV[,"traitmeanDr.f:traitmeanActive.m.animal"]/
  + sqrt(rep2cwt$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]*
           rep2cwt$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"])
posterior.mode(rfDR.mLA)
HPDinterval(rfDR.mLA,0.95)
