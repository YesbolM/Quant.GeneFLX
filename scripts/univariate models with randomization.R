library(MCMCglmm)

#Import dataset and subset by population
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

data <- read.csv("data.csv")
str(data)
data<-as.data.frame(data)
names(data)[1]<-"animal"
for(x in c(1:5,12:13))data[,x]<-as.factor(data[,x])
data$mLA.active<-ifelse(data$sex=="male",data$LA_active,NA)
data$mLA.passive<-ifelse(data$sex=="male",data$LA_passive,NA)
data$fLA.active<-ifelse(data$sex=="female",data$LA_active,NA)
data$fLA.passive<-ifelse(data$sex=="female",data$LA_passive,NA)
data$mWS<-ifelse(data$sex=="male",data$WS,NA)
data$fWS<-ifelse(data$sex=="female",data$WS,NA)
data$mDR<-ifelse(data$sex=="male",data$DR ,NA)
data$fDR<-ifelse(data$sex=="female",data$DR,NA)

rep1<-data[data$replicat=="rep1",]
FLX1<-rep1[rep1$treatment=="FLX",]
CFM1<-rep1[rep1$treatment=="CFM",]
CWT1<-rep1[rep1$treatment=="CWT",]

rep2<-data[data$replicat=="rep2",]
FLX2<-rep2[rep2$treatment=="FLX",]
CFM2<-rep2[rep2$treatment=="CFM",]
CWT2<-rep2[rep2$treatment=="CWT",]

#Import pedigree files
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

rep1CWT_p <- read.csv("data 2023 01 23/rep1CWT.p.csv")
ped3<-rep1CWT_p
head(ped3)
str(ped3)
ped3<-as.data.frame(ped3)

rep2FLX_p <- read.csv("rep2FLX.p.csv")
ped4<-rep2FLX_p
head(ped4)
str(ped4)
ped4<-as.data.frame(ped4)

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


#Rep1FLX
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(FLX1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                        family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
                        #For shorter run time:
                        #family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.mWS,file = "rep1FLX.mWS")

rep1FLX.mWS <- readRDS("rep1FLX.mWS")

#Male Va wing size
posterior.mode(rep1FLX.mWS$VCV[,"animal"])
HPDinterval(rep1FLX.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep1FLX.mWS$VCV[,"units"])
HPDinterval(rep1FLX.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep1FLX.mWS$VCV[,"animal"]/
  (rowSums(rep1FLX.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.fWS,file = "rep1FLX.fWS")

rep1FLX.fWS <- readRDS("rep1FLX.fWS")

#Female Va wing size
posterior.mode(rep1FLX.fWS$VCV[,"animal"])
HPDinterval(rep1FLX.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep1FLX.fWS$VCV[,"units"])
HPDinterval(rep1FLX.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep1FLX.fWS$VCV[,"animal"]/
  (rowSums(rep1FLX.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")


p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.mDR,file = "rep1FLX.mDR")

rep1FLX.mDR <- readRDS("rep1FLX.mDR")

#Male Va desiccation resistance
posterior.mode(rep1FLX.mDR$VCV[,"animal"])
HPDinterval(rep1FLX.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep1FLX.mDR$VCV[,"units"])
HPDinterval(rep1FLX.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep1FLX.mDR$VCV[,"animal"]/
  (rowSums(rep1FLX.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.fDR,file = "rep1FLX.fDR")

rep1FLX.fDR <- readRDS("rep1FLX.fDR")

#Female Va desiccation resistance
posterior.mode(rep1FLX.fDR$VCV[,"animal"])
HPDinterval(rep1FLX.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep1FLX.fDR$VCV[,"units"])
HPDinterval(rep1FLX.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep1FLX.fDR$VCV[,"animal"]/
  (rowSums(rep1FLX.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(FLX1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.mLA,file = "rep1FLX.mLA")

rep1FLX.mLA <- readRDS("rep1FLX.mLA")

#Male Va Locomotor activity
posterior.mode(rep1FLX.mLA$VCV[,"animal"])
HPDinterval(rep1FLX.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep1FLX.mLA$VCV[,"units"])
HPDinterval(rep1FLX.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep1FLX.mLA$VCV[,"animal"]/
  (rowSums(rep1FLX.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1FLX.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped1,data=FLX1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1FLX.fLA,file = "rep1FLX.fLA")

rep1FLX.fLA <- readRDS("rep1FLX.fLA")

#Female Va Locomotor activity
posterior.mode(rep1FLX.fLA$VCV[,"animal"])
HPDinterval(rep1FLX.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep1FLX.fLA$VCV[,"units"])
HPDinterval(rep1FLX.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep1FLX.fLA$VCV[,"animal"]/
  (rowSums(rep1FLX.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)


#Rep1CFM
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CFM1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.mWS,file = "rep1CFM.mWS")

rep1CFM.mWS <- readRDS("rep1CFM.mWS")

#Male Va wing size
posterior.mode(rep1CFM.mWS$VCV[,"animal"])
HPDinterval(rep1CFM.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep1CFM.mWS$VCV[,"units"])
HPDinterval(rep1CFM.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep1CFM.mWS$VCV[,"animal"]/
  (rowSums(rep1CFM.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.fWS,file = "rep1CFM.fWS")

rep1CFM.fWS <- readRDS("rep1CFM.fWS")

#Female Va wing size
posterior.mode(rep1CFM.fWS$VCV[,"animal"])
HPDinterval(rep1CFM.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep1CFM.fWS$VCV[,"units"])
HPDinterval(rep1CFM.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep1CFM.fWS$VCV[,"animal"]/
  (rowSums(rep1CFM.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.mDR,file = "rep1CFM.mDR")

rep1CFM.mDR <- readRDS("rep1CFM.mDR")

#Male Va desiccation resistance
posterior.mode(rep1CFM.mDR$VCV[,"animal"])
HPDinterval(rep1CFM.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep1CFM.mDR$VCV[,"units"])
HPDinterval(rep1CFM.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep1CFM.mDR$VCV[,"animal"]/
  (rowSums(rep1CFM.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.fDR,file = "rep1CFM.fDR")

rep1CFM.fDR <- readRDS("rep1CFM.fDR")

#Female Va desiccation resistance
posterior.mode(rep1CFM.fDR$VCV[,"animal"])
HPDinterval(rep1CFM.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep1CFM.fDR$VCV[,"units"])
HPDinterval(rep1CFM.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep1CFM.fDR$VCV[,"animal"]/
  (rowSums(rep1CFM.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CFM1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.mLA,file = "rep1CFM.mLA")

rep1CFM.mLA <- readRDS("rep1CFM.mLA")

#Male Va Locomotor activity
posterior.mode(rep1CFM.mLA$VCV[,"animal"])
HPDinterval(rep1CFM.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep1CFM.mLA$VCV[,"units"])
HPDinterval(rep1CFM.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep1CFM.mLA$VCV[,"animal"]/
  (rowSums(rep1CFM.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CFM.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped2,data=CFM1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CFM.fLA,file = "rep1CFM.fLA")

rep1CFM.fLA <- readRDS("rep1CFM.fLA")

#Female Va Locomotor activity
posterior.mode(rep1CFM.fLA$VCV[,"animal"])
HPDinterval(rep1CFM.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep1CFM.fLA$VCV[,"units"])
HPDinterval(rep1CFM.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep1CFM.fLA$VCV[,"animal"]/
  (rowSums(rep1CFM.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Rep1CWT
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CWT1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.mWS,file = "rep1CWT.mWS")

rep1CWT.mWS <- readRDS("rep1CWT.mWS")

#Male Va wing size
posterior.mode(rep1CWT.mWS$VCV[,"animal"])
HPDinterval(rep1CWT.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep1CWT.mWS$VCV[,"units"])
HPDinterval(rep1CWT.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep1CWT.mWS$VCV[,"animal"]/
  (rowSums(rep1CWT.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.fWS,file = "rep1CWT.fWS")

rep1CWT.fWS <- readRDS("rep1CWT.fWS")

#Female Va wing size
posterior.mode(rep1CWT.fWS$VCV[,"animal"])
HPDinterval(rep1CWT.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep1CWT.fWS$VCV[,"units"])
HPDinterval(rep1CWT.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep1CWT.fWS$VCV[,"animal"]/
  (rowSums(rep1CWT.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")


p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.mDR,file = "rep1CWT.mDR")

rep1CWT.mDR <- readRDS("rep1CWT.mDR")

#Male Va desiccation resistance
posterior.mode(rep1CWT.mDR$VCV[,"animal"])
HPDinterval(rep1CWT.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep1CWT.mDR$VCV[,"units"])
HPDinterval(rep1CWT.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep1CWT.mDR$VCV[,"animal"]/
  (rowSums(rep1CWT.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.fDR,file = "rep1CWT.fDR")

rep1CWT.fDR <- readRDS("rep1CWT.fDR")

#Female Va desiccation resistance
posterior.mode(rep1CWT.fDR$VCV[,"animal"])
HPDinterval(rep1CWT.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep1CWT.fDR$VCV[,"units"])
HPDinterval(rep1CWT.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep1CWT.fDR$VCV[,"animal"]/
  (rowSums(rep1CWT.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CWT1)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.mLA,file = "rep1CWT.mLA")

rep1CWT.mLA <- readRDS("rep1CWT.mLA")

#Male Va Locomotor activity
posterior.mode(rep1CWT.mLA$VCV[,"animal"])
HPDinterval(rep1CWT.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep1CWT.mLA$VCV[,"units"])
HPDinterval(rep1CWT.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep1CWT.mLA$VCV[,"animal"]/
  (rowSums(rep1CWT.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep1CWT.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped3,data=CWT1,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1CWT.fLA,file = "rep1CWT.fLA")

rep1CWT.fLA <- readRDS("rep1CWT.fLA")

#Female Va Locomotor activity
posterior.mode(rep1CWT.fLA$VCV[,"animal"])
HPDinterval(rep1CWT.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep1CWT.fLA$VCV[,"units"])
HPDinterval(rep1CWT.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep1CWT.fLA$VCV[,"animal"]/
  (rowSums(rep1CWT.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#rep2FLX
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(FLX2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.mWS,file = "rep2FLX.mWS")

rep2FLX.mWS <- readRDS("rep2FLX.mWS")

#Male Va wing size
posterior.mode(rep2FLX.mWS$VCV[,"animal"])
HPDinterval(rep2FLX.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep2FLX.mWS$VCV[,"units"])
HPDinterval(rep2FLX.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep2FLX.mWS$VCV[,"animal"]/
  (rowSums(rep2FLX.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.fWS,file = "rep2FLX.fWS")

rep2FLX.fWS <- readRDS("rep2FLX.fWS")

#Female Va wing size
posterior.mode(rep2FLX.fWS$VCV[,"animal"])
HPDinterval(rep2FLX.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep2FLX.fWS$VCV[,"units"])
HPDinterval(rep2FLX.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep2FLX.fWS$VCV[,"animal"]/
  (rowSums(rep2FLX.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")


p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.mDR,file = "rep2FLX.mDR")

rep2FLX.mDR <- readRDS("rep2FLX.mDR")

#Male Va desiccation resistance
posterior.mode(rep2FLX.mDR$VCV[,"animal"])
HPDinterval(rep2FLX.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep2FLX.mDR$VCV[,"units"])
HPDinterval(rep2FLX.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep2FLX.mDR$VCV[,"animal"]/
  (rowSums(rep2FLX.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.fDR,file = "rep2FLX.fDR")

rep2FLX.fDR <- readRDS("rep2FLX.fDR")

#Female Va desiccation resistance
posterior.mode(rep2FLX.fDR$VCV[,"animal"])
HPDinterval(rep2FLX.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep2FLX.fDR$VCV[,"units"])
HPDinterval(rep2FLX.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep2FLX.fDR$VCV[,"animal"]/
  (rowSums(rep2FLX.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(FLX2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.mLA,file = "rep2FLX.mLA")

rep2FLX.mLA <- readRDS("rep2FLX.mLA")

#Male Va Locomotor activity
posterior.mode(rep2FLX.mLA$VCV[,"animal"])
HPDinterval(rep2FLX.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep2FLX.mLA$VCV[,"units"])
HPDinterval(rep2FLX.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep2FLX.mLA$VCV[,"animal"]/
  (rowSums(rep2FLX.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2FLX.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped4,data=FLX2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2FLX.fLA,file = "rep2FLX.fLA")

rep2FLX.fLA <- readRDS("rep2FLX.fLA")

#Female Va Locomotor activity
posterior.mode(rep2FLX.fLA$VCV[,"animal"])
HPDinterval(rep2FLX.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep2FLX.fLA$VCV[,"units"])
HPDinterval(rep2FLX.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep2FLX.fLA$VCV[,"animal"]/
  (rowSums(rep2FLX.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)


#rep2CFM
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CFM2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.mWS,file = "rep2CFM.mWS")

rep2CFM.mWS <- readRDS("rep2CFM.mWS")

#Male Va wing size
posterior.mode(rep2CFM.mWS$VCV[,"animal"])
HPDinterval(rep2CFM.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep2CFM.mWS$VCV[,"units"])
HPDinterval(rep2CFM.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep2CFM.mWS$VCV[,"animal"]/
  (rowSums(rep2CFM.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.fWS,file = "rep2CFM.fWS")

rep2CFM.fWS <- readRDS("rep2CFM.fWS")

#Female Va wing size
posterior.mode(rep2CFM.fWS$VCV[,"animal"])
HPDinterval(rep2CFM.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep2CFM.fWS$VCV[,"units"])
HPDinterval(rep2CFM.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep2CFM.fWS$VCV[,"animal"]/
  (rowSums(rep2CFM.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")


p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.mDR,file = "rep2CFM.mDR")

rep2CFM.mDR <- readRDS("rep2CFM.mDR")

#Male Va desiccation resistance
posterior.mode(rep2CFM.mDR$VCV[,"animal"])
HPDinterval(rep2CFM.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep2CFM.mDR$VCV[,"units"])
HPDinterval(rep2CFM.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep2CFM.mDR$VCV[,"animal"]/
  (rowSums(rep2CFM.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.fDR,file = "rep2CFM.fDR")

rep2CFM.fDR <- readRDS("rep2CFM.fDR")

#Female Va desiccation resistance
posterior.mode(rep2CFM.fDR$VCV[,"animal"])
HPDinterval(rep2CFM.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep2CFM.fDR$VCV[,"units"])
HPDinterval(rep2CFM.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep2CFM.fDR$VCV[,"animal"]/
  (rowSums(rep2CFM.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CFM2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.mLA,file = "rep2CFM.mLA")

rep2CFM.mLA <- readRDS("rep2CFM.mLA")

#Male Va Locomotor activity
posterior.mode(rep2CFM.mLA$VCV[,"animal"])
HPDinterval(rep2CFM.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep2CFM.mLA$VCV[,"units"])
HPDinterval(rep2CFM.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep2CFM.mLA$VCV[,"animal"]/
  (rowSums(rep2CFM.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CFM.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped5,data=CFM2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CFM.fLA,file = "rep2CFM.fLA")

rep2CFM.fLA <- readRDS("rep2CFM.fLA")

#Female Va Locomotor activity
posterior.mode(rep2CFM.fLA$VCV[,"animal"])
HPDinterval(rep2CFM.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep2CFM.fLA$VCV[,"units"])
HPDinterval(rep2CFM.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep2CFM.fLA$VCV[,"animal"]/
  (rowSums(rep2CFM.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#rep2CWT
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CWT2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.mWS<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.mWS,file = "rep2CWT.mWS")

rep2CWT.mWS <- readRDS("rep2CWT.mWS")

#Male Va wing size
posterior.mode(rep2CWT.mWS$VCV[,"animal"])
HPDinterval(rep2CWT.mWS$VCV[,"animal"],0.95)
#Male Vr wing size
posterior.mode(rep2CWT.mWS$VCV[,"units"])
HPDinterval(rep2CWT.mWS$VCV[,"units"],0.95)
#Male heritability WS
h2.mWS <- rep2CWT.mWS$VCV[,"animal"]/
  (rowSums(rep2CWT.mWS$VCV))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.fWS<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.fWS,file = "rep2CWT.fWS")

rep2CWT.fWS <- readRDS("rep2CWT.fWS")

#Female Va wing size
posterior.mode(rep2CWT.fWS$VCV[,"animal"])
HPDinterval(rep2CWT.fWS$VCV[,"animal"],0.95)
#Female Vr wing size
posterior.mode(rep2CWT.fWS$VCV[,"units"])
HPDinterval(rep2CWT.fWS$VCV[,"units"],0.95)
#Female heritability WS
h2.fWS <- rep2CWT.fWS$VCV[,"animal"]/
  (rowSums(rep2CWT.fWS$VCV))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)


#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")


p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.mDR<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.mDR,file = "rep2CWT.mDR")

rep2CWT.mDR <- readRDS("rep2CWT.mDR")

#Male Va desiccation resistance
posterior.mode(rep2CWT.mDR$VCV[,"animal"])
HPDinterval(rep2CWT.mDR$VCV[,"animal"],0.95)
#Male Vr desiccation resistance
posterior.mode(rep2CWT.mDR$VCV[,"units"])
HPDinterval(rep2CWT.mDR$VCV[,"units"],0.95)
#Male heritability DR
h2.mDR <- rep2CWT.mDR$VCV[,"animal"]/
  (rowSums(rep2CWT.mDR$VCV))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.fDR<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="gaussian",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.fDR,file = "rep2CWT.fDR")

rep2CWT.fDR <- readRDS("rep2CWT.fDR")

#Female Va desiccation resistance
posterior.mode(rep2CWT.fDR$VCV[,"animal"])
HPDinterval(rep2CWT.fDR$VCV[,"animal"],0.95)
#Female Vr desiccation resistance
posterior.mode(rep2CWT.fDR$VCV[,"units"])
HPDinterval(rep2CWT.fDR$VCV[,"units"],0.95)
#Female heritability DR
h2.fDR <- rep2CWT.fDR$VCV[,"animal"]/
  (rowSums(rep2CWT.fDR$VCV))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)

#Locomotion
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

str(CWT2)
p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.mLA<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.mLA,file = "rep2CWT.mLA")

rep2CWT.mLA <- readRDS("rep2CWT.mLA")

#Male Va Locomotor activity
posterior.mode(rep2CWT.mLA$VCV[,"animal"])
HPDinterval(rep2CWT.mLA$VCV[,"animal"],0.95)
#Male Vr Locomotor activity
posterior.mode(rep2CWT.mLA$VCV[,"units"])
HPDinterval(rep2CWT.mLA$VCV[,"units"],0.95)
#Male heritability LA
h2.mLA <- rep2CWT.mLA$VCV[,"animal"]/
  (rowSums(rep2CWT.mLA$VCV))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)

p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
rep2CWT.fLA<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped6,data=CWT2,prior=p2,verbose=TRUE,
                      family="poisson",nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family="gaussian",nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2CWT.fLA,file = "rep2CWT.fLA")

rep2CWT.fLA <- readRDS("rep2CWT.fLA")

#Female Va Locomotor activity
posterior.mode(rep2CWT.fLA$VCV[,"animal"])
HPDinterval(rep2CWT.fLA$VCV[,"animal"],0.95)
#Female Vr Locomotor activity
posterior.mode(rep2CWT.fLA$VCV[,"units"])
HPDinterval(rep2CWT.fLA$VCV[,"units"],0.95)
#Female heritability LA
h2.fLA <- rep2CWT.fLA$VCV[,"animal"]/
  (rowSums(rep2CWT.fLA$VCV))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)

#Randomizations
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")

#FLX1 male wing size
rep1FLX.mWS <- readRDS("rep1FLX.mWS")
rep1FLX.mWS$VCV
rep1FLX.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.mWS.rand.VCV) <- dimnames(rep1FLX.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
    #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped1,data=FLX1m.rand,prior=p2,verbose=FALSE,
                        family="gaussian",nitt=550000,thin=2000,burnin=150000)
  #For shorter run time:
  #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1FLX.mWS.rand.VCV[i,] <- posterior.mode(rep1FLX.mWS.rand$VCV)
  write.csv(rep1FLX.mWS.rand.VCV, "rep1FLX.mWS.rand.VCV.csv")
}

rep1FLX.mWS.rand.VCV <- read.csv("rep1FLX.mWS.rand.VCV.csv")
head(rep1FLX.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep1FLX.mWS$VCV[,"animal"])
HPDinterval(rep1FLX.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.mWS.rand.VCV[,"animal"],0.95))
sum(rep1FLX.mWS.rand.VCV$"animal">=posterior.mode(rep1FLX.mWS$VCV[,"animal"]))
hist(rep1FLX.mWS.rand.VCV$"animal")


#FLX1 female wing size
rep1FLX.fWS <- readRDS("rep1FLX.fWS")
rep1FLX.fWS$VCV
rep1FLX.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.fWS.rand.VCV) <- dimnames(rep1FLX.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped1,data=FLX1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1FLX.fWS.rand.VCV[i,] <- posterior.mode(rep1FLX.fWS.rand$VCV)
  write.csv(rep1FLX.fWS.rand.VCV, "rep1FLX.fWS.rand.VCV.csv")
}


rep1FLX.fWS.rand.VCV <- read.csv("rep1FLX.fWS.rand.VCV.csv")
head(rep1FLX.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep1FLX.fWS$VCV[,"animal"])
HPDinterval(rep1FLX.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.fWS.rand.VCV[,"animal"],0.95))
sum(rep1FLX.fWS.rand.VCV$"animal">=posterior.mode(rep1FLX.fWS$VCV[,"animal"]))
hist(rep1FLX.fWS.rand.VCV$"animal")


#FLX1 male DR
rep1FLX.mDR <- readRDS("rep1FLX.mDR")
rep1FLX.mDR$VCV
rep1FLX.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.mDR.rand.VCV) <- dimnames(rep1FLX.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped1,data=FLX1m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1FLX.mDR.rand.VCV[i,] <- posterior.mode(rep1FLX.mDR.rand$VCV)
  write.csv(rep1FLX.mDR.rand.VCV, "rep1FLX.mDR.rand.VCV.csv")
}


rep1FLX.mDR.rand.VCV <- read.csv("rep1FLX.mDR.rand.VCV.csv")
head(rep1FLX.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep1FLX.mDR$VCV[,"animal"])
HPDinterval(rep1FLX.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.mDR.rand.VCV[,"animal"],0.95))
sum(rep1FLX.mDR.rand.VCV$"animal">=posterior.mode(rep1FLX.mDR$VCV[,"animal"]))
hist(rep1FLX.mDR.rand.VCV$"animal")


#FLX1 female DR
rep1FLX.fDR <- readRDS("rep1FLX.fDR")
rep1FLX.fDR$VCV
rep1FLX.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.fDR.rand.VCV) <- dimnames(rep1FLX.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped1,data=FLX1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1FLX.fDR.rand.VCV[i,] <- posterior.mode(rep1FLX.fDR.rand$VCV)
  write.csv(rep1FLX.fDR.rand.VCV, "rep1FLX.fDR.rand.VCV.csv")
}


rep1FLX.fDR.rand.VCV <- read.csv("rep1FLX.fDR.rand.VCV.csv")
head(rep1FLX.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep1FLX.fDR$VCV[,"animal"])
HPDinterval(rep1FLX.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.fDR.rand.VCV[,"animal"],0.95))
sum(rep1FLX.fDR.rand.VCV$"animal">=posterior.mode(rep1FLX.fDR$VCV[,"animal"]))
hist(rep1FLX.fDR.rand.VCV$"animal")


#FLX1 male LA
rep1FLX.mLA <- readRDS("rep1FLX.mLA")
rep1FLX.mLA$VCV
rep1FLX.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.mLA.rand.VCV) <- dimnames(rep1FLX.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped1,data=FLX1m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1FLX.mLA.rand.VCV[i,] <- posterior.mode(rep1FLX.mLA.rand$VCV)
  write.csv(rep1FLX.mLA.rand.VCV, "rep1FLX.mLA.rand.VCV.csv")
}


rep1FLX.mLA.rand.VCV <- read.csv("rep1FLX.mLA.rand.VCV.csv")
head(rep1FLX.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep1FLX.mLA$VCV[,"animal"])
HPDinterval(rep1FLX.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.mLA.rand.VCV[,"animal"],0.95))
sum(rep1FLX.mLA.rand.VCV$"animal">=posterior.mode(rep1FLX.mLA$VCV[,"animal"]))
hist(rep1FLX.mLA.rand.VCV$"animal")


#FLX1 female LA
rep1FLX.fLA <- readRDS("rep1FLX.fLA")
rep1FLX.fLA$VCV
rep1FLX.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1FLX.fLA.rand.VCV) <- dimnames(rep1FLX.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1FLX.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped1,data=FLX1f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1FLX.fLA.rand.VCV[i,] <- posterior.mode(rep1FLX.fLA.rand$VCV)
  write.csv(rep1FLX.fLA.rand.VCV, "rep1FLX.fLA.rand.VCV.csv")
}


rep1FLX.fLA.rand.VCV <- read.csv("rep1FLX.fLA.rand.VCV.csv")
head(rep1FLX.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep1FLX.fLA$VCV[,"animal"])
HPDinterval(rep1FLX.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1FLX.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1FLX.fLA.rand.VCV[,"animal"],0.95))
sum(rep1FLX.fLA.rand.VCV$"animal">=posterior.mode(rep1FLX.fLA$VCV[,"animal"]))
hist(rep1FLX.fLA.rand.VCV$"animal")


#CFM1 male wing size
rep1CFM.mWS <- readRDS("rep1CFM.mWS")
rep1CFM.mWS$VCV
rep1CFM.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.mWS.rand.VCV) <- dimnames(rep1CFM.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped2,data=CFM1m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CFM.mWS.rand.VCV[i,] <- posterior.mode(rep1CFM.mWS.rand$VCV)
  write.csv(rep1CFM.mWS.rand.VCV, "rep1CFM.mWS.rand.VCV.csv")
}


rep1CFM.mWS.rand.VCV <- read.csv("rep1CFM.mWS.rand.VCV.csv")
head(rep1CFM.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep1CFM.mWS$VCV[,"animal"])
HPDinterval(rep1CFM.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.mWS.rand.VCV[,"animal"],0.95))
sum(rep1CFM.mWS.rand.VCV$"animal">=posterior.mode(rep1CFM.mWS$VCV[,"animal"]))
hist(rep1CFM.mWS.rand.VCV$"animal")


#CFM1 female wing size
rep1CFM.fWS <- readRDS("rep1CFM.fWS")
rep1CFM.fWS$VCV
rep1CFM.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.fWS.rand.VCV) <- dimnames(rep1CFM.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped2,data=CFM1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CFM.fWS.rand.VCV[i,] <- posterior.mode(rep1CFM.fWS.rand$VCV)
  write.csv(rep1CFM.fWS.rand.VCV, "rep1CFM.fWS.rand.VCV.csv")
}


rep1CFM.fWS.rand.VCV <- read.csv("rep1CFM.fWS.rand.VCV.csv")
head(rep1CFM.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep1CFM.fWS$VCV[,"animal"])
HPDinterval(rep1CFM.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.fWS.rand.VCV[,"animal"],0.95))
sum(rep1CFM.fWS.rand.VCV$"animal">=posterior.mode(rep1CFM.fWS$VCV[,"animal"]))
hist(rep1CFM.fWS.rand.VCV$"animal")


#CFM1 male DR
rep1CFM.mDR <- readRDS("rep1CFM.mDR")
rep1CFM.mDR$VCV
rep1CFM.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.mDR.rand.VCV) <- dimnames(rep1CFM.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped2,data=CFM1m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CFM.mDR.rand.VCV[i,] <- posterior.mode(rep1CFM.mDR.rand$VCV)
  write.csv(rep1CFM.mDR.rand.VCV, "rep1CFM.mDR.rand.VCV.csv")
}


rep1CFM.mDR.rand.VCV <- read.csv("rep1CFM.mDR.rand.VCV.csv")
head(rep1CFM.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep1CFM.mDR$VCV[,"animal"])
HPDinterval(rep1CFM.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.mDR.rand.VCV[,"animal"],0.95))
sum(rep1CFM.mDR.rand.VCV$"animal">=posterior.mode(rep1CFM.mDR$VCV[,"animal"]))
hist(rep1CFM.mDR.rand.VCV$"animal")


#CFM1 female DR
rep1CFM.fDR <- readRDS("rep1CFM.fDR")
rep1CFM.fDR$VCV
rep1CFM.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.fDR.rand.VCV) <- dimnames(rep1CFM.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped2,data=CFM1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CFM.fDR.rand.VCV[i,] <- posterior.mode(rep1CFM.fDR.rand$VCV)
  write.csv(rep1CFM.fDR.rand.VCV, "rep1CFM.fDR.rand.VCV.csv")
}


rep1CFM.fDR.rand.VCV <- read.csv("rep1CFM.fDR.rand.VCV.csv")
head(rep1CFM.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep1CFM.fDR$VCV[,"animal"])
HPDinterval(rep1CFM.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.fDR.rand.VCV[,"animal"],0.95))
sum(rep1CFM.fDR.rand.VCV$"animal">=posterior.mode(rep1CFM.fDR$VCV[,"animal"]))
hist(rep1CFM.fDR.rand.VCV$"animal")


#CFM1 male LA
rep1CFM.mLA <- readRDS("rep1CFM.mLA")
rep1CFM.mLA$VCV
rep1CFM.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.mLA.rand.VCV) <- dimnames(rep1CFM.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped2,data=CFM1m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1CFM.mLA.rand.VCV[i,] <- posterior.mode(rep1CFM.mLA.rand$VCV)
  write.csv(rep1CFM.mLA.rand.VCV, "rep1CFM.mLA.rand.VCV.csv")
}


rep1CFM.mLA.rand.VCV <- read.csv("rep1CFM.mLA.rand.VCV.csv")
head(rep1CFM.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep1CFM.mLA$VCV[,"animal"])
HPDinterval(rep1CFM.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.mLA.rand.VCV[,"animal"],0.95))
sum(rep1CFM.mLA.rand.VCV$"animal">=posterior.mode(rep1CFM.mLA$VCV[,"animal"]))
hist(rep1CFM.mLA.rand.VCV$"animal")
#Male Vr LA
posterior.mode(rep1CFM.mLA$VCV[,"units"])
HPDinterval(rep1CFM.mLA$VCV[,"units"],0.95)
sum(rep1CFM.mLA.rand.VCV$"units"<=posterior.mode(rep1CFM.mLA$VCV[,"units"]))
hist(rep1CFM.mLA.rand.VCV$"units")
#Male heritability LA
h2.mLA <- rep1CFM.mLA$VCV[,"animal"]/(sum(posterior.mode(rep1CFM.mLA$VCV)))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
h2.mLA.rand <- rep1CFM.mLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CFM.mLA.rand.VCV[,2:4]))))
hist(h2.mLA.rand)
sum(h2.mLA.rand>=posterior.mode(h2.mLA))
#All parameters p<0.01.

#CFM1 female LA
rep1CFM.fLA <- readRDS("rep1CFM.fLA")
rep1CFM.fLA$VCV
rep1CFM.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1CFM.fLA.rand.VCV) <- dimnames(rep1CFM.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CFM.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped2,data=CFM1f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1CFM.fLA.rand.VCV[i,] <- posterior.mode(rep1CFM.fLA.rand$VCV)
  write.csv(rep1CFM.fLA.rand.VCV, "rep1CFM.fLA.rand.VCV.csv")
}


rep1CFM.fLA.rand.VCV <- read.csv("rep1CFM.fLA.rand.VCV.csv")
head(rep1CFM.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep1CFM.fLA$VCV[,"animal"])
HPDinterval(rep1CFM.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CFM.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CFM.fLA.rand.VCV[,"animal"],0.95))
sum(rep1CFM.fLA.rand.VCV$"animal">=posterior.mode(rep1CFM.fLA$VCV[,"animal"]))
hist(rep1CFM.fLA.rand.VCV$"animal")
#Female Vr LA
posterior.mode(rep1CFM.fLA$VCV[,"units"])
HPDinterval(rep1CFM.fLA$VCV[,"units"],0.95)
sum(rep1CFM.fLA.rand.VCV$"units"<=posterior.mode(rep1CFM.fLA$VCV[,"units"]))
hist(rep1CFM.fLA.rand.VCV$"units")
#Female heritability LA
h2.fLA <- rep1CFM.fLA$VCV[,"animal"]/(sum(posterior.mode(rep1CFM.fLA$VCV)))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
h2.fLA.rand <- rep1CFM.fLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CFM.fLA.rand.VCV[,2:4]))))
hist(h2.fLA.rand)
sum(h2.fLA.rand>=posterior.mode(h2.fLA))
#All parameters p<0.01.

#CWT1 male wing size
rep1CWT.mWS <- readRDS("rep1CWT.mWS")
rep1CWT.mWS$VCV
rep1CWT.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.mWS.rand.VCV) <- dimnames(rep1CWT.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped3,data=CWT1m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CWT.mWS.rand.VCV[i,] <- posterior.mode(rep1CWT.mWS.rand$VCV)
  write.csv(rep1CWT.mWS.rand.VCV, "rep1CWT.mWS.rand.VCV.csv")
}


rep1CWT.mWS.rand.VCV <- read.csv("rep1CWT.mWS.rand.VCV.csv")
head(rep1CWT.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep1CWT.mWS$VCV[,"animal"])
HPDinterval(rep1CWT.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.mWS.rand.VCV[,"animal"],0.95))
sum(rep1CWT.mWS.rand.VCV$"animal">=posterior.mode(rep1CWT.mWS$VCV[,"animal"]))
hist(rep1CWT.mWS.rand.VCV$"animal")
#Male Vr wing size
posterior.mode(rep1CWT.mWS$VCV[,"units"])
HPDinterval(rep1CWT.mWS$VCV[,"units"],0.95)
sum(rep1CWT.mWS.rand.VCV$"units"<=posterior.mode(rep1CWT.mWS$VCV[,"units"]))
hist(rep1CWT.mWS.rand.VCV$"units")
#Male heritability WS
h2.mWS <- rep1CWT.mWS$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.mWS$VCV)))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
h2.mWS.rand <- rep1CWT.mWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.mWS.rand.VCV[,2:4]))))
hist(h2.mWS.rand)
sum(h2.mWS.rand>=posterior.mode(h2.mWS))
#All parameters p<0.01.

#CWT1 female wing size
rep1CWT.fWS <- readRDS("rep1CWT.fWS")
rep1CWT.fWS$VCV
rep1CWT.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.fWS.rand.VCV) <- dimnames(rep1CWT.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped3,data=CWT1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CWT.fWS.rand.VCV[i,] <- posterior.mode(rep1CWT.fWS.rand$VCV)
  write.csv(rep1CWT.fWS.rand.VCV, "rep1CWT.fWS.rand.VCV.csv")
}


rep1CWT.fWS.rand.VCV <- read.csv("rep1CWT.fWS.rand.VCV.csv")
head(rep1CWT.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep1CWT.fWS$VCV[,"animal"])
HPDinterval(rep1CWT.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.fWS.rand.VCV[,"animal"],0.95))
sum(rep1CWT.fWS.rand.VCV$"animal">=posterior.mode(rep1CWT.fWS$VCV[,"animal"]))
hist(rep1CWT.fWS.rand.VCV$"animal")
#Female Vr wing size
posterior.mode(rep1CWT.fWS$VCV[,"units"])
HPDinterval(rep1CWT.fWS$VCV[,"units"],0.95)
sum(rep1CWT.fWS.rand.VCV$"units"<=posterior.mode(rep1CWT.fWS$VCV[,"units"]))
hist(rep1CWT.fWS.rand.VCV$"units")
#Female heritability WS
h2.fWS <- rep1CWT.fWS$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.fWS$VCV)))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
h2.fWS.rand <- rep1CWT.fWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.fWS.rand.VCV[,2:4]))))
hist(h2.fWS.rand)
sum(h2.fWS.rand>=posterior.mode(h2.fWS))
#All parameters p<0.01.

#CWT1 male DR
rep1CWT.mDR <- readRDS("rep1CWT.mDR")
rep1CWT.mDR$VCV
rep1CWT.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.mDR.rand.VCV) <- dimnames(rep1CWT.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped3,data=CWT1m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CWT.mDR.rand.VCV[i,] <- posterior.mode(rep1CWT.mDR.rand$VCV)
  write.csv(rep1CWT.mDR.rand.VCV, "rep1CWT.mDR.rand.VCV.csv")
}


rep1CWT.mDR.rand.VCV <- read.csv("rep1CWT.mDR.rand.VCV.csv")
head(rep1CWT.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep1CWT.mDR$VCV[,"animal"])
HPDinterval(rep1CWT.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.mDR.rand.VCV[,"animal"],0.95))
sum(rep1CWT.mDR.rand.VCV$"animal">=posterior.mode(rep1CWT.mDR$VCV[,"animal"]))
hist(rep1CWT.mDR.rand.VCV$"animal")
#Male Vr DR
posterior.mode(rep1CWT.mDR$VCV[,"units"])
HPDinterval(rep1CWT.mDR$VCV[,"units"],0.95)
sum(rep1CWT.mDR.rand.VCV$"units"<=posterior.mode(rep1CWT.mDR$VCV[,"units"]))
hist(rep1CWT.mDR.rand.VCV$"units")
#Male heritability DR
h2.mDR <- rep1CWT.mDR$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.mDR$VCV)))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
h2.mDR.rand <- rep1CWT.mDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.mDR.rand.VCV[,2:4]))))
hist(h2.mDR.rand)
sum(h2.mDR.rand>=posterior.mode(h2.mDR))
#All parameters p<0.01.

#CWT1 female DR
rep1CWT.fDR <- readRDS("rep1CWT.fDR")
rep1CWT.fDR$VCV
rep1CWT.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.fDR.rand.VCV) <- dimnames(rep1CWT.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped3,data=CWT1f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep1CWT.fDR.rand.VCV[i,] <- posterior.mode(rep1CWT.fDR.rand$VCV)
  write.csv(rep1CWT.fDR.rand.VCV, "rep1CWT.fDR.rand.VCV.csv")
}


rep1CWT.fDR.rand.VCV <- read.csv("rep1CWT.fDR.rand.VCV.csv")
head(rep1CWT.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep1CWT.fDR$VCV[,"animal"])
HPDinterval(rep1CWT.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.fDR.rand.VCV[,"animal"],0.95))
sum(rep1CWT.fDR.rand.VCV$"animal">=posterior.mode(rep1CWT.fDR$VCV[,"animal"]))
hist(rep1CWT.fDR.rand.VCV$"animal")
#Female Vr DR
posterior.mode(rep1CWT.fDR$VCV[,"units"])
HPDinterval(rep1CWT.fDR$VCV[,"units"],0.95)
sum(rep1CWT.fDR.rand.VCV$"units"<=posterior.mode(rep1CWT.fDR$VCV[,"units"]))
hist(rep1CWT.fDR.rand.VCV$"units")
#Female heritability DR
h2.fDR <- rep1CWT.fDR$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.fDR$VCV)))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
h2.fDR.rand <- rep1CWT.fDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.fDR.rand.VCV[,2:4]))))
hist(h2.fDR.rand)
sum(h2.fDR.rand>=posterior.mode(h2.fDR))
#All parameters p<0.01.

#CWT1 male LA
rep1CWT.mLA <- readRDS("rep1CWT.mLA")
rep1CWT.mLA$VCV
rep1CWT.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.mLA.rand.VCV) <- dimnames(rep1CWT.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT1m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped3,data=CWT1m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1CWT.mLA.rand.VCV[i,] <- posterior.mode(rep1CWT.mLA.rand$VCV)
  write.csv(rep1CWT.mLA.rand.VCV, "rep1CWT.mLA.rand.VCV.csv")
}


rep1CWT.mLA.rand.VCV <- read.csv("rep1CWT.mLA.rand.VCV.csv")
head(rep1CWT.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep1CWT.mLA$VCV[,"animal"])
HPDinterval(rep1CWT.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.mLA.rand.VCV[,"animal"],0.95))
sum(rep1CWT.mLA.rand.VCV$"animal">=posterior.mode(rep1CWT.mLA$VCV[,"animal"]))
hist(rep1CWT.mLA.rand.VCV$"animal")
#Male Vr LA
posterior.mode(rep1CWT.mLA$VCV[,"units"])
HPDinterval(rep1CWT.mLA$VCV[,"units"],0.95)
sum(rep1CWT.mLA.rand.VCV$"units"<=posterior.mode(rep1CWT.mLA$VCV[,"units"]))
hist(rep1CWT.mLA.rand.VCV$"units")
#Male heritability LA
h2.mLA <- rep1CWT.mLA$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.mLA$VCV)))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
h2.mLA.rand <- rep1CWT.mLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.mLA.rand.VCV[,2:4]))))
hist(h2.mLA.rand)
sum(h2.mLA.rand>=posterior.mode(h2.mLA))
#All parameters p<0.01.

#CWT1 female LA
rep1CWT.fLA <- readRDS("rep1CWT.fLA")
rep1CWT.fLA$VCV
rep1CWT.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep1CWT.fLA.rand.VCV) <- dimnames(rep1CWT.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT1f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT1f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep1CWT.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped3,data=CWT1f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep1CWT.fLA.rand.VCV[i,] <- posterior.mode(rep1CWT.fLA.rand$VCV)
  write.csv(rep1CWT.fLA.rand.VCV, "rep1CWT.fLA.rand.VCV.csv")
}


rep1CWT.fLA.rand.VCV <- read.csv("rep1CWT.fLA.rand.VCV.csv")
head(rep1CWT.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep1CWT.fLA$VCV[,"animal"])
HPDinterval(rep1CWT.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep1CWT.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep1CWT.fLA.rand.VCV[,"animal"],0.95))
sum(rep1CWT.fLA.rand.VCV$"animal">=posterior.mode(rep1CWT.fLA$VCV[,"animal"]))
hist(rep1CWT.fLA.rand.VCV$"animal")
#Female Vr LA
posterior.mode(rep1CWT.fLA$VCV[,"units"])
HPDinterval(rep1CWT.fLA$VCV[,"units"],0.95)
sum(rep1CWT.fLA.rand.VCV$"units"<=posterior.mode(rep1CWT.fLA$VCV[,"units"]))
hist(rep1CWT.fLA.rand.VCV$"units")
#Female heritability LA
h2.fLA <- rep1CWT.fLA$VCV[,"animal"]/(sum(posterior.mode(rep1CWT.fLA$VCV)))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
h2.fLA.rand <- rep1CWT.fLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep1CWT.fLA.rand.VCV[,2:4]))))
hist(h2.fLA.rand)
sum(h2.fLA.rand>=posterior.mode(h2.fLA))
#All parameters p<0.01.

#FLX2 male wing size
rep2FLX.mWS <- readRDS("rep2FLX.mWS")
rep2FLX.mWS$VCV
rep2FLX.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.mWS.rand.VCV) <- dimnames(rep2FLX.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped4,data=FLX2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2FLX.mWS.rand.VCV[i,] <- posterior.mode(rep2FLX.mWS.rand$VCV)
  write.csv(rep2FLX.mWS.rand.VCV, "rep2FLX.mWS.rand.VCV.csv")
}


rep2FLX.mWS.rand.VCV <- read.csv("rep2FLX.mWS.rand.VCV.csv")
head(rep2FLX.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep2FLX.mWS$VCV[,"animal"])
HPDinterval(rep2FLX.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.mWS.rand.VCV[,"animal"],0.95))
sum(rep2FLX.mWS.rand.VCV$"animal">=posterior.mode(rep2FLX.mWS$VCV[,"animal"]))
hist(rep2FLX.mWS.rand.VCV$"animal")
#Male Vr wing size
posterior.mode(rep2FLX.mWS$VCV[,"units"])
HPDinterval(rep2FLX.mWS$VCV[,"units"],0.95)
sum(rep2FLX.mWS.rand.VCV$"units"<=posterior.mode(rep2FLX.mWS$VCV[,"units"]))
hist(rep2FLX.mWS.rand.VCV$"units")
#Male heritability WS
h2.mWS <- rep2FLX.mWS$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.mWS$VCV)))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
h2.mWS.rand <- rep2FLX.mWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.mWS.rand.VCV[,2:4]))))
hist(h2.mWS.rand)
sum(h2.mWS.rand>=posterior.mode(h2.mWS))
#All parameters p<0.01.

#FLX2 female wing size
rep2FLX.fWS <- readRDS("rep2FLX.fWS")
rep2FLX.fWS$VCV
rep2FLX.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.fWS.rand.VCV) <- dimnames(rep2FLX.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped4,data=FLX2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2FLX.fWS.rand.VCV[i,] <- posterior.mode(rep2FLX.fWS.rand$VCV)
  write.csv(rep2FLX.fWS.rand.VCV, "rep2FLX.fWS.rand.VCV.csv")
}


rep2FLX.fWS.rand.VCV <- read.csv("rep2FLX.fWS.rand.VCV.csv")
head(rep2FLX.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep2FLX.fWS$VCV[,"animal"])
HPDinterval(rep2FLX.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.fWS.rand.VCV[,"animal"],0.95))
sum(rep2FLX.fWS.rand.VCV$"animal">=posterior.mode(rep2FLX.fWS$VCV[,"animal"]))
hist(rep2FLX.fWS.rand.VCV$"animal")
#Female Vr wing size
posterior.mode(rep2FLX.fWS$VCV[,"units"])
HPDinterval(rep2FLX.fWS$VCV[,"units"],0.95)
sum(rep2FLX.fWS.rand.VCV$"units"<=posterior.mode(rep2FLX.fWS$VCV[,"units"]))
hist(rep2FLX.fWS.rand.VCV$"units")
#Female heritability WS
h2.fWS <- rep2FLX.fWS$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.fWS$VCV)))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
h2.fWS.rand <- rep2FLX.fWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.fWS.rand.VCV[,2:4]))))
hist(h2.fWS.rand)
sum(h2.fWS.rand>=posterior.mode(h2.fWS))
#All parameters p<0.01.

#FLX2 male DR
rep2FLX.mDR <- readRDS("rep2FLX.mDR")
rep2FLX.mDR$VCV
rep2FLX.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.mDR.rand.VCV) <- dimnames(rep2FLX.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped4,data=FLX2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2FLX.mDR.rand.VCV[i,] <- posterior.mode(rep2FLX.mDR.rand$VCV)
  write.csv(rep2FLX.mDR.rand.VCV, "rep2FLX.mDR.rand.VCV.csv")
}


rep2FLX.mDR.rand.VCV <- read.csv("rep2FLX.mDR.rand.VCV.csv")
head(rep2FLX.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep2FLX.mDR$VCV[,"animal"])
HPDinterval(rep2FLX.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.mDR.rand.VCV[,"animal"],0.95))
sum(rep2FLX.mDR.rand.VCV$"animal">=posterior.mode(rep2FLX.mDR$VCV[,"animal"]))
hist(rep2FLX.mDR.rand.VCV$"animal")
#Male Vr DR
posterior.mode(rep2FLX.mDR$VCV[,"units"])
HPDinterval(rep2FLX.mDR$VCV[,"units"],0.95)
sum(rep2FLX.mDR.rand.VCV$"units"<=posterior.mode(rep2FLX.mDR$VCV[,"units"]))
hist(rep2FLX.mDR.rand.VCV$"units")
#Male heritability DR
h2.mDR <- rep2FLX.mDR$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.mDR$VCV)))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
h2.mDR.rand <- rep2FLX.mDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.mDR.rand.VCV[,2:4]))))
hist(h2.mDR.rand)
sum(h2.mDR.rand>=posterior.mode(h2.mDR))
#All parameters p<0.01.

#FLX2 female DR
rep2FLX.fDR <- readRDS("rep2FLX.fDR")
rep2FLX.fDR$VCV
rep2FLX.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.fDR.rand.VCV) <- dimnames(rep2FLX.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped4,data=FLX2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2FLX.fDR.rand.VCV[i,] <- posterior.mode(rep2FLX.fDR.rand$VCV)
  write.csv(rep2FLX.fDR.rand.VCV, "rep2FLX.fDR.rand.VCV.csv")
}


rep2FLX.fDR.rand.VCV <- read.csv("rep2FLX.fDR.rand.VCV.csv")
head(rep2FLX.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep2FLX.fDR$VCV[,"animal"])
HPDinterval(rep2FLX.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.fDR.rand.VCV[,"animal"],0.95))
sum(rep2FLX.fDR.rand.VCV$"animal">=posterior.mode(rep2FLX.fDR$VCV[,"animal"]))
hist(rep2FLX.fDR.rand.VCV$"animal")
#Female Vr DR
posterior.mode(rep2FLX.fDR$VCV[,"units"])
HPDinterval(rep2FLX.fDR$VCV[,"units"],0.95)
sum(rep2FLX.fDR.rand.VCV$"units"<=posterior.mode(rep2FLX.fDR$VCV[,"units"]))
hist(rep2FLX.fDR.rand.VCV$"units")
#Female heritability DR
h2.fDR <- rep2FLX.fDR$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.fDR$VCV)))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
h2.fDR.rand <- rep2FLX.fDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.fDR.rand.VCV[,2:4]))))
hist(h2.fDR.rand)
sum(h2.fDR.rand>=posterior.mode(h2.fDR))
#All parameters p<0.01.

#FLX2 male LA
rep2FLX.mLA <- readRDS("rep2FLX.mLA")
rep2FLX.mLA$VCV
rep2FLX.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.mLA.rand.VCV) <- dimnames(rep2FLX.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  FLX2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped4,data=FLX2m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2FLX.mLA.rand.VCV[i,] <- posterior.mode(rep2FLX.mLA.rand$VCV)
  write.csv(rep2FLX.mLA.rand.VCV, "rep2FLX.mLA.rand.VCV.csv")
}


rep2FLX.mLA.rand.VCV <- read.csv("rep2FLX.mLA.rand.VCV.csv")
head(rep2FLX.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep2FLX.mLA$VCV[,"animal"])
HPDinterval(rep2FLX.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.mLA.rand.VCV[,"animal"],0.95))
sum(rep2FLX.mLA.rand.VCV$"animal">=posterior.mode(rep2FLX.mLA$VCV[,"animal"]))
hist(rep2FLX.mLA.rand.VCV$"animal")
#Male Vr LA
posterior.mode(rep2FLX.mLA$VCV[,"units"])
HPDinterval(rep2FLX.mLA$VCV[,"units"],0.95)
sum(rep2FLX.mLA.rand.VCV$"units"<=posterior.mode(rep2FLX.mLA$VCV[,"units"]))
hist(rep2FLX.mLA.rand.VCV$"units")
#Male heritability LA
h2.mLA <- rep2FLX.mLA$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.mLA$VCV)))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
h2.mLA.rand <- rep2FLX.mLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.mLA.rand.VCV[,2:4]))))
hist(h2.mLA.rand)
sum(h2.mLA.rand>=posterior.mode(h2.mLA))
#All parameters p<0.01.

#FLX2 female LA
rep2FLX.fLA <- readRDS("rep2FLX.fLA")
rep2FLX.fLA$VCV
rep2FLX.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2FLX.fLA.rand.VCV) <- dimnames(rep2FLX.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  FLX2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(FLX2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2FLX.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped4,data=FLX2f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2FLX.fLA.rand.VCV[i,] <- posterior.mode(rep2FLX.fLA.rand$VCV)
  write.csv(rep2FLX.fLA.rand.VCV, "rep2FLX.fLA.rand.VCV.csv")
}


rep2FLX.fLA.rand.VCV <- read.csv("rep2FLX.fLA.rand.VCV.csv")
head(rep2FLX.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep2FLX.fLA$VCV[,"animal"])
HPDinterval(rep2FLX.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2FLX.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2FLX.fLA.rand.VCV[,"animal"],0.95))
sum(rep2FLX.fLA.rand.VCV$"animal">=posterior.mode(rep2FLX.fLA$VCV[,"animal"]))
hist(rep2FLX.fLA.rand.VCV$"animal")
#Female Vr LA
posterior.mode(rep2FLX.fLA$VCV[,"units"])
HPDinterval(rep2FLX.fLA$VCV[,"units"],0.95)
sum(rep2FLX.fLA.rand.VCV$"units"<=posterior.mode(rep2FLX.fLA$VCV[,"units"]))
hist(rep2FLX.fLA.rand.VCV$"units")
#Female heritability LA
h2.fLA <- rep2FLX.fLA$VCV[,"animal"]/(sum(posterior.mode(rep2FLX.fLA$VCV)))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
h2.fLA.rand <- rep2FLX.fLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2FLX.fLA.rand.VCV[,2:4]))))
hist(h2.fLA.rand)
sum(h2.fLA.rand>=posterior.mode(h2.fLA))
#All parameters p<0.01.

#CFM2 male wing size
rep2CFM.mWS <- readRDS("rep2CFM.mWS")
rep2CFM.mWS$VCV
rep2CFM.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.mWS.rand.VCV) <- dimnames(rep2CFM.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped5,data=CFM2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CFM.mWS.rand.VCV[i,] <- posterior.mode(rep2CFM.mWS.rand$VCV)
  write.csv(rep2CFM.mWS.rand.VCV, "rep2CFM.mWS.rand.VCV.csv")
}


rep2CFM.mWS.rand.VCV <- read.csv("rep2CFM.mWS.rand.VCV.csv")
head(rep2CFM.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep2CFM.mWS$VCV[,"animal"])
HPDinterval(rep2CFM.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.mWS.rand.VCV[,"animal"],0.95))
sum(rep2CFM.mWS.rand.VCV$"animal">=posterior.mode(rep2CFM.mWS$VCV[,"animal"]))
hist(rep2CFM.mWS.rand.VCV$"animal")
#Male Vr wing size
posterior.mode(rep2CFM.mWS$VCV[,"units"])
HPDinterval(rep2CFM.mWS$VCV[,"units"],0.95)
sum(rep2CFM.mWS.rand.VCV$"units"<=posterior.mode(rep2CFM.mWS$VCV[,"units"]))
hist(rep2CFM.mWS.rand.VCV$"units")
#Male heritability WS
h2.mWS <- rep2CFM.mWS$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.mWS$VCV)))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
h2.mWS.rand <- rep2CFM.mWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.mWS.rand.VCV[,2:4]))))
hist(h2.mWS.rand)
sum(h2.mWS.rand>=posterior.mode(h2.mWS))
#All parameters p<0.01.

#CFM2 female wing size
rep2CFM.fWS <- readRDS("rep2CFM.fWS")
rep2CFM.fWS$VCV
rep2CFM.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.fWS.rand.VCV) <- dimnames(rep2CFM.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped5,data=CFM2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CFM.fWS.rand.VCV[i,] <- posterior.mode(rep2CFM.fWS.rand$VCV)
  write.csv(rep2CFM.fWS.rand.VCV, "rep2CFM.fWS.rand.VCV.csv")
}


rep2CFM.fWS.rand.VCV <- read.csv("rep2CFM.fWS.rand.VCV.csv")
head(rep2CFM.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep2CFM.fWS$VCV[,"animal"])
HPDinterval(rep2CFM.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.fWS.rand.VCV[,"animal"],0.95))
sum(rep2CFM.fWS.rand.VCV$"animal">=posterior.mode(rep2CFM.fWS$VCV[,"animal"]))
hist(rep2CFM.fWS.rand.VCV$"animal")
#Female Vr wing size
posterior.mode(rep2CFM.fWS$VCV[,"units"])
HPDinterval(rep2CFM.fWS$VCV[,"units"],0.95)
sum(rep2CFM.fWS.rand.VCV$"units"<=posterior.mode(rep2CFM.fWS$VCV[,"units"]))
hist(rep2CFM.fWS.rand.VCV$"units")
#Female heritability WS
h2.fWS <- rep2CFM.fWS$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.fWS$VCV)))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
h2.fWS.rand <- rep2CFM.fWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.fWS.rand.VCV[,2:4]))))
hist(h2.fWS.rand)
sum(h2.fWS.rand>=posterior.mode(h2.fWS))
#All parameters p<0.01.

#CFM2 male DR
rep2CFM.mDR <- readRDS("rep2CFM.mDR")
rep2CFM.mDR$VCV
rep2CFM.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.mDR.rand.VCV) <- dimnames(rep2CFM.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped5,data=CFM2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CFM.mDR.rand.VCV[i,] <- posterior.mode(rep2CFM.mDR.rand$VCV)
  write.csv(rep2CFM.mDR.rand.VCV, "rep2CFM.mDR.rand.VCV.csv")
}


rep2CFM.mDR.rand.VCV <- read.csv("rep2CFM.mDR.rand.VCV.csv")
head(rep2CFM.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep2CFM.mDR$VCV[,"animal"])
HPDinterval(rep2CFM.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.mDR.rand.VCV[,"animal"],0.95))
sum(rep2CFM.mDR.rand.VCV$"animal">=posterior.mode(rep2CFM.mDR$VCV[,"animal"]))
hist(rep2CFM.mDR.rand.VCV$"animal")
#Male Vr DR
posterior.mode(rep2CFM.mDR$VCV[,"units"])
HPDinterval(rep2CFM.mDR$VCV[,"units"],0.95)
sum(rep2CFM.mDR.rand.VCV$"units"<=posterior.mode(rep2CFM.mDR$VCV[,"units"]))
hist(rep2CFM.mDR.rand.VCV$"units")
#Male heritability DR
h2.mDR <- rep2CFM.mDR$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.mDR$VCV)))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
h2.mDR.rand <- rep2CFM.mDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.mDR.rand.VCV[,2:4]))))
hist(h2.mDR.rand)
sum(h2.mDR.rand>=posterior.mode(h2.mDR))
#All parameters p<0.01.

#CFM2 female DR
rep2CFM.fDR <- readRDS("rep2CFM.fDR")
rep2CFM.fDR$VCV
rep2CFM.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.fDR.rand.VCV) <- dimnames(rep2CFM.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped5,data=CFM2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CFM.fDR.rand.VCV[i,] <- posterior.mode(rep2CFM.fDR.rand$VCV)
  write.csv(rep2CFM.fDR.rand.VCV, "rep2CFM.fDR.rand.VCV.csv")
}


rep2CFM.fDR.rand.VCV <- read.csv("rep2CFM.fDR.rand.VCV.csv")
head(rep2CFM.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep2CFM.fDR$VCV[,"animal"])
HPDinterval(rep2CFM.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.fDR.rand.VCV[,"animal"],0.95))
sum(rep2CFM.fDR.rand.VCV$"animal">=posterior.mode(rep2CFM.fDR$VCV[,"animal"]))
hist(rep2CFM.fDR.rand.VCV$"animal")
#Female Vr DR
posterior.mode(rep2CFM.fDR$VCV[,"units"])
HPDinterval(rep2CFM.fDR$VCV[,"units"],0.95)
sum(rep2CFM.fDR.rand.VCV$"units"<=posterior.mode(rep2CFM.fDR$VCV[,"units"]))
hist(rep2CFM.fDR.rand.VCV$"units")
#Female heritability DR
h2.fDR <- rep2CFM.fDR$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.fDR$VCV)))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
h2.fDR.rand <- rep2CFM.fDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.fDR.rand.VCV[,2:4]))))
hist(h2.fDR.rand)
sum(h2.fDR.rand>=posterior.mode(h2.fDR))
#All parameters p<0.01.

#CFM2 male LA
rep2CFM.mLA <- readRDS("rep2CFM.mLA")
rep2CFM.mLA$VCV
rep2CFM.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.mLA.rand.VCV) <- dimnames(rep2CFM.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CFM2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped5,data=CFM2m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2CFM.mLA.rand.VCV[i,] <- posterior.mode(rep2CFM.mLA.rand$VCV)
  write.csv(rep2CFM.mLA.rand.VCV, "rep2CFM.mLA.rand.VCV.csv")
}


rep2CFM.mLA.rand.VCV <- read.csv("rep2CFM.mLA.rand.VCV.csv")
head(rep2CFM.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep2CFM.mLA$VCV[,"animal"])
HPDinterval(rep2CFM.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.mLA.rand.VCV[,"animal"],0.95))
sum(rep2CFM.mLA.rand.VCV$"animal">=posterior.mode(rep2CFM.mLA$VCV[,"animal"]))
hist(rep2CFM.mLA.rand.VCV$"animal")
#Male Vr LA
posterior.mode(rep2CFM.mLA$VCV[,"units"])
HPDinterval(rep2CFM.mLA$VCV[,"units"],0.95)
sum(rep2CFM.mLA.rand.VCV$"units"<=posterior.mode(rep2CFM.mLA$VCV[,"units"]))
hist(rep2CFM.mLA.rand.VCV$"units")
#Male heritability LA
h2.mLA <- rep2CFM.mLA$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.mLA$VCV)))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
h2.mLA.rand <- rep2CFM.mLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.mLA.rand.VCV[,2:4]))))
hist(h2.mLA.rand)
sum(h2.mLA.rand>=posterior.mode(h2.mLA))
#All parameters p<0.01.

#CFM2 female LA
rep2CFM.fLA <- readRDS("rep2CFM.fLA")
rep2CFM.fLA$VCV
rep2CFM.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2CFM.fLA.rand.VCV) <- dimnames(rep2CFM.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CFM2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CFM2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CFM.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped5,data=CFM2f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2CFM.fLA.rand.VCV[i,] <- posterior.mode(rep2CFM.fLA.rand$VCV)
  write.csv(rep2CFM.fLA.rand.VCV, "rep2CFM.fLA.rand.VCV.csv")
}


rep2CFM.fLA.rand.VCV <- read.csv("rep2CFM.fLA.rand.VCV.csv")
head(rep2CFM.fLA.rand.VCV)
str(rep2CFM.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep2CFM.fLA$VCV[,"animal"])
HPDinterval(rep2CFM.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CFM.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CFM.fLA.rand.VCV[,"animal"],0.95))
sum(rep2CFM.fLA.rand.VCV$"animal">=posterior.mode(rep2CFM.fLA$VCV[,"animal"]))
hist(rep2CFM.fLA.rand.VCV$"animal")
#Female Vr LA
posterior.mode(rep2CFM.fLA$VCV[,"units"])
HPDinterval(rep2CFM.fLA$VCV[,"units"],0.95)
sum(rep2CFM.fLA.rand.VCV$"units"<=posterior.mode(rep2CFM.fLA$VCV[,"units"]))
hist(rep2CFM.fLA.rand.VCV$"units")
#Female heritability LA
h2.fLA <- rep2CFM.fLA$VCV[,"animal"]/(sum(posterior.mode(rep2CFM.fLA$VCV)))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
h2.fLA.rand <- rep2CFM.fLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CFM.fLA.rand.VCV[,2:4]))))
hist(h2.fLA.rand)
sum(h2.fLA.rand>=posterior.mode(h2.fLA))
#All parameters p<0.01.

#CWT2 male wing size
rep2CWT.mWS <- readRDS("rep2CWT.mWS")
rep2CWT.mWS$VCV
rep2CWT.mWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.mWS.rand.VCV) <- dimnames(rep2CWT.mWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.mWS.rand<-MCMCglmm(mWS~round, random=~animal+round.sec,pedigree=ped6,data=CWT2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CWT.mWS.rand.VCV[i,] <- posterior.mode(rep2CWT.mWS.rand$VCV)
  write.csv(rep2CWT.mWS.rand.VCV, "rep2CWT.mWS.rand.VCV.csv")
}


rep2CWT.mWS.rand.VCV <- read.csv("rep2CWT.mWS.rand.VCV.csv")
head(rep2CWT.mWS.rand.VCV)

#Significance testing
#Male Va wing size
posterior.mode(rep2CWT.mWS$VCV[,"animal"])
HPDinterval(rep2CWT.mWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.mWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.mWS.rand.VCV[,"animal"],0.95))
sum(rep2CWT.mWS.rand.VCV$"animal">=posterior.mode(rep2CWT.mWS$VCV[,"animal"]))
hist(rep2CWT.mWS.rand.VCV$"animal")
#Male Vr wing size
posterior.mode(rep2CWT.mWS$VCV[,"units"])
HPDinterval(rep2CWT.mWS$VCV[,"units"],0.95)
sum(rep2CWT.mWS.rand.VCV$"units"<=posterior.mode(rep2CWT.mWS$VCV[,"units"]))
hist(rep2CWT.mWS.rand.VCV$"units")
#Male heritability WS
h2.mWS <- rep2CWT.mWS$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.mWS$VCV)))
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
h2.mWS.rand <- rep2CWT.mWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.mWS.rand.VCV[,2:4]))))
hist(h2.mWS.rand)
sum(h2.mWS.rand>=posterior.mode(h2.mWS))
#All parameters p<0.01.

#CWT2 female wing size
rep2CWT.fWS <- readRDS("rep2CWT.fWS")
rep2CWT.fWS$VCV
rep2CWT.fWS.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.fWS.rand.VCV) <- dimnames(rep2CWT.fWS$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.fWS.rand<-MCMCglmm(fWS~round, random=~animal+round.sec,pedigree=ped6,data=CWT2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CWT.fWS.rand.VCV[i,] <- posterior.mode(rep2CWT.fWS.rand$VCV)
  write.csv(rep2CWT.fWS.rand.VCV, "rep2CWT.fWS.rand.VCV.csv")
}


rep2CWT.fWS.rand.VCV <- read.csv("rep2CWT.fWS.rand.VCV.csv")
head(rep2CWT.fWS.rand.VCV)

#Significance testing
#Female Va wing size
posterior.mode(rep2CWT.fWS$VCV[,"animal"])
HPDinterval(rep2CWT.fWS$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.fWS.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.fWS.rand.VCV[,"animal"],0.95))
sum(rep2CWT.fWS.rand.VCV$"animal">=posterior.mode(rep2CWT.fWS$VCV[,"animal"]))
hist(rep2CWT.fWS.rand.VCV$"animal")
#Female Vr wing size
posterior.mode(rep2CWT.fWS$VCV[,"units"])
HPDinterval(rep2CWT.fWS$VCV[,"units"],0.95)
sum(rep2CWT.fWS.rand.VCV$"units"<=posterior.mode(rep2CWT.fWS$VCV[,"units"]))
hist(rep2CWT.fWS.rand.VCV$"units")
#Female heritability WS
h2.fWS <- rep2CWT.fWS$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.fWS$VCV)))
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
h2.fWS.rand <- rep2CWT.fWS.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.fWS.rand.VCV[,2:4]))))
hist(h2.fWS.rand)
sum(h2.fWS.rand>=posterior.mode(h2.fWS))
#All parameters p<0.01.

#CWT2 male DR
rep2CWT.mDR <- readRDS("rep2CWT.mDR")
rep2CWT.mDR$VCV
rep2CWT.mDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.mDR.rand.VCV) <- dimnames(rep2CWT.mDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.mDR.rand<-MCMCglmm(mDR~round, random=~animal+round.sec,pedigree=ped6,data=CWT2m.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CWT.mDR.rand.VCV[i,] <- posterior.mode(rep2CWT.mDR.rand$VCV)
  write.csv(rep2CWT.mDR.rand.VCV, "rep2CWT.mDR.rand.VCV.csv")
}


rep2CWT.mDR.rand.VCV <- read.csv("rep2CWT.mDR.rand.VCV.csv")
head(rep2CWT.mDR.rand.VCV)

#Significance testing
#Male Va DR
posterior.mode(rep2CWT.mDR$VCV[,"animal"])
HPDinterval(rep2CWT.mDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.mDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.mDR.rand.VCV[,"animal"],0.95))
sum(rep2CWT.mDR.rand.VCV$"animal">=posterior.mode(rep2CWT.mDR$VCV[,"animal"]))
hist(rep2CWT.mDR.rand.VCV$"animal")
#Male Vr DR
posterior.mode(rep2CWT.mDR$VCV[,"units"])
HPDinterval(rep2CWT.mDR$VCV[,"units"],0.95)
sum(rep2CWT.mDR.rand.VCV$"units"<=posterior.mode(rep2CWT.mDR$VCV[,"units"]))
hist(rep2CWT.mDR.rand.VCV$"units")
#Male heritability DR
h2.mDR <- rep2CWT.mDR$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.mDR$VCV)))
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
h2.mDR.rand <- rep2CWT.mDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.mDR.rand.VCV[,2:4]))))
hist(h2.mDR.rand)
sum(h2.mDR.rand>=posterior.mode(h2.mDR))
#All parameters p<0.01.

#CWT2 female DR
rep2CWT.fDR <- readRDS("rep2CWT.fDR")
rep2CWT.fDR$VCV
rep2CWT.fDR.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.fDR.rand.VCV) <- dimnames(rep2CWT.fDR$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.fDR.rand<-MCMCglmm(fDR~round, random=~animal+round.sec,pedigree=ped6,data=CWT2f.rand,prior=p2,verbose=FALSE,
                             family="gaussian",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="gaussian",nitt=1000,thin=10,burnin=100) 
  rep2CWT.fDR.rand.VCV[i,] <- posterior.mode(rep2CWT.fDR.rand$VCV)
  write.csv(rep2CWT.fDR.rand.VCV, "rep2CWT.fDR.rand.VCV.csv")
}


rep2CWT.fDR.rand.VCV <- read.csv("rep2CWT.fDR.rand.VCV.csv")
head(rep2CWT.fDR.rand.VCV)

#Significance testing
#Fmeale Va DR
posterior.mode(rep2CWT.fDR$VCV[,"animal"])
HPDinterval(rep2CWT.fDR$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.fDR.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.fDR.rand.VCV[,"animal"],0.95))
sum(rep2CWT.fDR.rand.VCV$"animal">=posterior.mode(rep2CWT.fDR$VCV[,"animal"]))
hist(rep2CWT.fDR.rand.VCV$"animal")
#Female Vr DR
posterior.mode(rep2CWT.fDR$VCV[,"units"])
HPDinterval(rep2CWT.fDR$VCV[,"units"],0.95)
sum(rep2CWT.fDR.rand.VCV$"units"<=posterior.mode(rep2CWT.fDR$VCV[,"units"]))
hist(rep2CWT.fDR.rand.VCV$"units")
#Female heritability DR
h2.fDR <- rep2CWT.fDR$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.fDR$VCV)))
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
h2.fDR.rand <- rep2CWT.fDR.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.fDR.rand.VCV[,2:4]))))
hist(h2.fDR.rand)
sum(h2.fDR.rand>=posterior.mode(h2.fDR))
#All parameters p<0.01.

#CWT2 male LA
rep2CWT.mLA <- readRDS("rep2CWT.mLA")
rep2CWT.mLA$VCV
rep2CWT.mLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.mLA.rand.VCV) <- dimnames(rep2CWT.mLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mDR <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mDR <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mDR <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mDR <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  #Combine randomized subsets
  CWT2m.rand <- rbind(r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2m.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.mLA.rand<-MCMCglmm(mLA.active~round, random=~animal+round.sec,pedigree=ped6,data=CWT2m.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2CWT.mLA.rand.VCV[i,] <- posterior.mode(rep2CWT.mLA.rand$VCV)
  write.csv(rep2CWT.mLA.rand.VCV, "rep2CWT.mLA.rand.VCV.csv")
}


rep2CWT.mLA.rand.VCV <- read.csv("rep2CWT.mLA.rand.VCV.csv")
head(rep2CWT.mLA.rand.VCV)

#Significance testing
#Male Va LA
posterior.mode(rep2CWT.mLA$VCV[,"animal"])
HPDinterval(rep2CWT.mLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.mLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.mLA.rand.VCV[,"animal"],0.95))
sum(rep2CWT.mLA.rand.VCV$"animal">=posterior.mode(rep2CWT.mLA$VCV[,"animal"]))
hist(rep2CWT.mLA.rand.VCV$"animal")
#Male Vr LA
posterior.mode(rep2CWT.mLA$VCV[,"units"])
HPDinterval(rep2CWT.mLA$VCV[,"units"],0.95)
sum(rep2CWT.mLA.rand.VCV$"units"<=posterior.mode(rep2CWT.mLA$VCV[,"units"]))
hist(rep2CWT.mLA.rand.VCV$"units")
#Male heritability LA
h2.mLA <- rep2CWT.mLA$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.mLA$VCV)))
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
h2.mLA.rand <- rep2CWT.mLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.mLA.rand.VCV[,2:4]))))
hist(h2.mLA.rand)
sum(h2.mLA.rand>=posterior.mode(h2.mLA))
#All parameters p<0.01.

#CWT2 female LA
rep2CWT.fLA <- readRDS("rep2CWT.fLA")
rep2CWT.fLA$VCV
rep2CWT.fLA.rand.VCV <- matrix(NA,200,3)
colnames(rep2CWT.fLA.rand.VCV) <- dimnames(rep2CWT.fLA$VCV)[[2]]
for (i in 1:200) {
  print(i)
  #4 rounds, 1 sex = 4 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$fDR <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active <- r1F[rand.order,17]
  r1F.rand$fLA.passive <- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$fDR <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active <- r2F[rand.order,17]
  r2F.rand$fLA.passive <- r2F[rand.order,18]
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$fDR <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active <- r3F[rand.order,17]
  r3F.rand$fLA.passive <- r3F[rand.order,18]
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$fDR <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active <- r4F[rand.order,17]
  r4F.rand$fLA.passive <- r4F[rand.order,18]
  #Combine randomized subsets
  CWT2f.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand)
  #str(CWT2f.rand)
  
  #Analysis
  p2<-list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1,nu = 0.002), G2 = list(V = 1, nu = 0.002)))
  rep2CWT.fLA.rand<-MCMCglmm(fLA.active~round, random=~animal+round.sec,pedigree=ped6,data=CWT2f.rand,prior=p2,verbose=FALSE,
                             family="poisson",nitt=550000,thin=2000,burnin=150000)
                             #For shorter run time:
                             #family="poisson",nitt=1000,thin=10,burnin=100) 
  rep2CWT.fLA.rand.VCV[i,] <- posterior.mode(rep2CWT.fLA.rand$VCV)
  write.csv(rep2CWT.fLA.rand.VCV, "rep2CWT.fLA.rand.VCV.csv")
}


rep2CWT.fLA.rand.VCV <- read.csv("rep2CWT.fLA.rand.VCV.csv")
head(rep2CWT.fLA.rand.VCV)

#Significance testing
#Fmeale Va LA
posterior.mode(rep2CWT.fLA$VCV[,"animal"])
HPDinterval(rep2CWT.fLA$VCV[,"animal"],0.95)
posterior.mode(as.mcmc(rep2CWT.fLA.rand.VCV[,"animal"]))
HPDinterval(as.mcmc(rep2CWT.fLA.rand.VCV[,"animal"],0.95))
sum(rep2CWT.fLA.rand.VCV$"animal">=posterior.mode(rep2CWT.fLA$VCV[,"animal"]))
hist(rep2CWT.fLA.rand.VCV$"animal")
#Female Vr LA
posterior.mode(rep2CWT.fLA$VCV[,"units"])
HPDinterval(rep2CWT.fLA$VCV[,"units"],0.95)
sum(rep2CWT.fLA.rand.VCV$"units"<=posterior.mode(rep2CWT.fLA$VCV[,"units"]))
hist(rep2CWT.fLA.rand.VCV$"units")
#Female heritability LA
h2.fLA <- rep2CWT.fLA$VCV[,"animal"]/(sum(posterior.mode(rep2CWT.fLA$VCV)))
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
h2.fLA.rand <- rep2CWT.fLA.rand.VCV$"animal"/(sum(posterior.mode(as.mcmc(rep2CWT.fLA.rand.VCV[,2:4]))))
hist(h2.fLA.rand)
sum(h2.fLA.rand>=posterior.mode(h2.fLA))
#All parameters p<0.01.