library(MCMCglmm)

#Import data
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

data <- read.csv("data.csv")
str(data)
data<-as.data.frame(data) # 6858 obs. of  14 variables
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

str(FLX1)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1flx.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                        random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1,prior=p2,verbose=FALSE,
                        family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
                        #For shorter run time:
                        #family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1flx.WS,file = "rep1flx.WS")

rep1flx.WS <- readRDS("rep1flx.WS")

#Randomization
rep1flx.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep1flx.WS.rand.VCV) <- dimnames(rep1flx.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. 
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1flx.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                          random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1.rand,prior=p2,verbose=FALSE,
                          family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
                          #For shorter run time:
                          #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1flx.WS.rand.VCV[i,] <- posterior.mode(rep1flx.WS.rand$VCV)
  write.csv(rep1flx.WS.rand.VCV, "rep1flx.WS.rand.VCV.csv")
}

rep1flx.WS.rand.VCV <- read.csv("rep1flx.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep1flx.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep1flx.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep1flx.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep1flx.WS$VCV[,"traitmWS.units"])
HPDinterval(rep1flx.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"]+
    rep1flx.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep1flx.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep1flx.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep1flx.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep1flx.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep1flx.WS$VCV[,"traitfWS.units"])
HPDinterval(rep1flx.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"]+
    rep1flx.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep1flx.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep1flx.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep1flx.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep1flx.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep1flx.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep1flx.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep1flx.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep1flx.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep1flx.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep1flx.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep1flx.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep1flx.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep1flx.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1flx.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
                    #For shorter run time:
                    #family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1flx.DR, file="rep1flx.DR")

rep1flx.DR <- readRDS("rep1flx.DR")

#Randomization
rep1flx.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep1flx.DR.rand.VCV) <- dimnames(rep1flx.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1flx.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
                            #For shorter run time:
                            #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1flx.DR.rand.VCV[i,] <- posterior.mode(rep1flx.DR.rand$VCV)
  write.csv(rep1flx.DR.rand.VCV, "rep1flx.DR.rand.VCV.csv")
}

rep1flx.DR.rand.VCV <- read.csv("rep1flx.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep1flx.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep1flx.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep1flx.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep1flx.DR$VCV[,"traitmDR.units"])
HPDinterval(rep1flx.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep1flx.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep1flx.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep1flx.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep1flx.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep1flx.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep1flx.DR$VCV[,"traitfDR.units"])
HPDinterval(rep1flx.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep1flx.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep1flx.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1flx.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep1flx.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep1flx.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep1flx.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep1flx.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep1flx.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep1flx.DR$VCV[,"traitmDR:traitfDR.animal"]/
   + sqrt(rep1flx.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep1flx.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep1flx.DR.rand.VCV$traitmDR.traitfDR.animal/
   + sqrt(rep1flx.DR.rand.VCV$traitmDR.traitmDR.animal*
           rep1flx.DR.rand.VCV$traitfDR.traitfDR.animal)
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
str(FLX1)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1flx.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
                    #For shorter run time:
                    #family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1flx.LA,file = "rep1flx.LA")

rep1flx.LA <- readRDS("rep1flx.LA")

#Randomization
rep1flx.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep1flx.LA.rand.VCV) <- dimnames(rep1flx.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1flx.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped1,data=FLX1.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
                            #For shorter run time:
                            #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1flx.LA.rand.VCV[i,] <- posterior.mode(rep1flx.LA.rand$VCV)
  write.csv(rep1flx.LA.rand.VCV, "rep1flx.LA.rand.VCV.csv")
}

rep1flx.LA.rand.VCV <- read.csv("rep1flx.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep1flx.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep1flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep1flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep1flx.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep1flx.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep1flx.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep1flx.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep1flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep1flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep1flx.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep1flx.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep1flx.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep1flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep1flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep1flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep1flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep1flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep1flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep1flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep1flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep1flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep1flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep1flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)

#Rep1CFM
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

str(CFM1)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cfm.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
                    #For shorter run time:
                    #family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cfm.WS,file = "rep1cfm.WS")

rep1cfm.WS <- readRDS("rep1cfm.WS")

#Randomization
rep1cfm.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep1cfm.WS.rand.VCV) <- dimnames(rep1cfm.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cfm.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
                            #For shorter run time:
                            #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cfm.WS.rand.VCV[i,] <- posterior.mode(rep1cfm.WS.rand$VCV)
  write.csv(rep1cfm.WS.rand.VCV, "rep1cfm.WS.rand.VCV.csv")
}

rep1cfm.WS.rand.VCV <- read.csv("rep1cfm.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep1cfm.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep1cfm.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep1cfm.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep1cfm.WS$VCV[,"traitmWS.units"])
HPDinterval(rep1cfm.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"]+
     rep1cfm.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep1cfm.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep1cfm.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep1cfm.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep1cfm.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep1cfm.WS$VCV[,"traitfWS.units"])
HPDinterval(rep1cfm.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"]+
     rep1cfm.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep1cfm.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep1cfm.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep1cfm.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep1cfm.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep1cfm.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep1cfm.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep1cfm.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep1cfm.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep1cfm.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep1cfm.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep1cfm.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep1cfm.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep1cfm.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cfm.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
                    #For shorter run time:
                    #family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cfm.DR, file="rep1cfm.DR")

rep1cfm.DR <- readRDS("rep1cfm.DR")

#Randomization
rep1cfm.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep1cfm.DR.rand.VCV) <- dimnames(rep1cfm.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cfm.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cfm.DR.rand.VCV[i,] <- posterior.mode(rep1cfm.DR.rand$VCV)
  write.csv(rep1cfm.DR.rand.VCV, "rep1cfm.DR.rand.VCV.csv")
}

rep1cfm.DR.rand.VCV <- read.csv("rep1cfm.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep1cfm.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep1cfm.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep1cfm.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep1cfm.DR$VCV[,"traitmDR.units"])
HPDinterval(rep1cfm.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep1cfm.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep1cfm.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep1cfm.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep1cfm.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep1cfm.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep1cfm.DR$VCV[,"traitfDR.units"])
HPDinterval(rep1cfm.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep1cfm.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep1cfm.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1cfm.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep1cfm.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep1cfm.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep1cfm.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep1cfm.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep1cfm.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep1cfm.DR$VCV[,"traitmDR:traitfDR.animal"]/
  + sqrt(rep1cfm.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep1cfm.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep1cfm.DR.rand.VCV$"traitmDR.traitfDR.animal"/
  + sqrt(rep1cfm.DR.rand.VCV$"traitmDR.traitmDR.animal"*
           rep1cfm.DR.rand.VCV$"traitfDR.traitfDR.animal")
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cfm.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cfm.LA,file = "rep1cfm.LA")

rep1cfm.LA <- readRDS("rep1cfm.LA")

#Randomization
rep1cfm.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep1cfm.LA.rand.VCV) <- dimnames(rep1cfm.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cfm.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped2,data=CFM1.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cfm.LA.rand.VCV[i,] <- posterior.mode(rep1cfm.LA.rand$VCV)
  write.csv(rep1cfm.LA.rand.VCV, "rep1cfm.LA.rand.VCV.csv")
}

rep1cfm.LA.rand.VCV <- read.csv("rep1cfm.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep1cfm.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep1cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep1cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep1cfm.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep1cfm.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep1cfm.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep1cfm.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep1cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep1cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep1cfm.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep1cfm.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep1cfm.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep1cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep1cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep1cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep1cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep1cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep1cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep1cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep1cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep1cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep1cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep1cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)

#Rep1CWT
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cwt.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cwt.WS,file = "rep1cwt.WS")

rep1cwt.WS <- readRDS("rep1cwt.WS")

#Randomization
rep1cwt.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep1cwt.WS.rand.VCV) <- dimnames(rep1cwt.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cwt.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cwt.WS.rand.VCV[i,] <- posterior.mode(rep1cwt.WS.rand$VCV)
  write.csv(rep1cwt.WS.rand.VCV, "rep1cwt.WS.rand.VCV.csv")
}

rep1cwt.WS.rand.VCV <- read.csv("rep1cwt.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep1cwt.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep1cwt.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep1cwt.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep1cwt.WS$VCV[,"traitmWS.units"])
HPDinterval(rep1cwt.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"]+
     rep1cwt.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep1cwt.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep1cwt.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep1cwt.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep1cwt.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep1cwt.WS$VCV[,"traitfWS.units"])
HPDinterval(rep1cwt.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"]+
     rep1cwt.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep1cwt.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep1cwt.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep1cwt.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep1cwt.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep1cwt.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep1cwt.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep1cwt.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep1cwt.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep1cwt.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep1cwt.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep1cwt.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep1cwt.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep1cwt.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cwt.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cwt.DR, file="rep1cwt.DR")

rep1cwt.DR <- readRDS("rep1cwt.DR")

#Randomization
rep1cwt.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep1cwt.DR.rand.VCV) <- dimnames(rep1cwt.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cwt.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cwt.DR.rand.VCV[i,] <- posterior.mode(rep1cwt.DR.rand$VCV)
  write.csv(rep1cwt.DR.rand.VCV, "rep1cwt.DR.rand.VCV.csv")
}

rep1cwt.DR.rand.VCV <- read.csv("rep1cwt.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep1cwt.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep1cwt.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep1cwt.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep1cwt.DR$VCV[,"traitmDR.units"])
HPDinterval(rep1cwt.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep1cwt.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep1cwt.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep1cwt.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep1cwt.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep1cwt.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep1cwt.DR$VCV[,"traitfDR.units"])
HPDinterval(rep1cwt.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep1cwt.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep1cwt.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep1cwt.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep1cwt.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep1cwt.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep1cwt.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep1cwt.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep1cwt.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep1cwt.DR$VCV[,"traitmDR:traitfDR.animal"]/
  + sqrt(rep1cwt.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep1cwt.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep1cwt.DR.rand.VCV$"traitmDR.traitfDR.animal"/
  + sqrt(rep1cwt.DR.rand.VCV$"traitmDR.traitmDR.animal"*
           rep1cwt.DR.rand.VCV$"traitfDR.traitfDR.animal")
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep1cwt.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cwt.LA,file = "rep1cwt.LA")

rep1cwt.LA <- readRDS("rep1cwt.LA")

#Randomization
rep1cwt.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep1cwt.LA.rand.VCV) <- dimnames(rep1cwt.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT1.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep1cwt.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped3,data=CWT1.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep1cwt.LA.rand.VCV[i,] <- posterior.mode(rep1cwt.LA.rand$VCV)
  write.csv(rep1cwt.LA.rand.VCV, "rep1cwt.LA.rand.VCV.csv")
}

rep1cwt.LA.rand.VCV <- read.csv("rep1cwt.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep1cwt.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep1cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep1cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep1cwt.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep1cwt.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep1cwt.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep1cwt.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep1cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep1cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep1cwt.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep1cwt.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep1cwt.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep1cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep1cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep1cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep1cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep1cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep1cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep1cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep1cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep1cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep1cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep1cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep1cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep1cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)

#Rep2FLX
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

#str(FLX2)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2flx.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2flx.WS,file = "rep2flx.WS")

rep2flx.WS <- readRDS("rep2flx.WS")

#Randomization
rep2flx.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep2flx.WS.rand.VCV) <- dimnames(rep2flx.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2flx.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2flx.WS.rand.VCV[i,] <- posterior.mode(rep2flx.WS.rand$VCV)
  write.csv(rep2flx.WS.rand.VCV, "rep2flx.WS.rand.VCV.csv")
}

rep2flx.WS.rand.VCV <- read.csv("rep2flx.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep2flx.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep2flx.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep2flx.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep2flx.WS$VCV[,"traitmWS.units"])
HPDinterval(rep2flx.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"]+
     rep2flx.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep2flx.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep2flx.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep2flx.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep2flx.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep2flx.WS$VCV[,"traitfWS.units"])
HPDinterval(rep2flx.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"]+
     rep2flx.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep2flx.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep2flx.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep2flx.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep2flx.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep2flx.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep2flx.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep2flx.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep2flx.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep2flx.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep2flx.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep2flx.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep2flx.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep2flx.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2flx.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2flx.DR, file="rep2flx.DR")

rep2flx.DR <- readRDS("rep2flx.DR")

#Randomization
rep2flx.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep2flx.DR.rand.VCV) <- dimnames(rep2flx.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2flx.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2flx.DR.rand.VCV[i,] <- posterior.mode(rep2flx.DR.rand$VCV)
  write.csv(rep2flx.DR.rand.VCV, "rep2flx.DR.rand.VCV.csv")
}

rep2flx.DR.rand.VCV <- read.csv("rep2flx.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep2flx.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep2flx.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep2flx.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep2flx.DR$VCV[,"traitmDR.units"])
HPDinterval(rep2flx.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep2flx.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep2flx.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep2flx.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep2flx.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep2flx.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep2flx.DR$VCV[,"traitfDR.units"])
HPDinterval(rep2flx.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep2flx.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep2flx.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2flx.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep2flx.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep2flx.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep2flx.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep2flx.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep2flx.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep2flx.DR$VCV[,"traitmDR:traitfDR.animal"]/
  + sqrt(rep2flx.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep2flx.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep2flx.DR.rand.VCV$traitmDR.traitfDR.animal/
  + sqrt(rep2flx.DR.rand.VCV$traitmDR.traitmDR.animal*
           rep2flx.DR.rand.VCV$traitfDR.traitfDR.animal)
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
#str(FLX2)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2flx.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2flx.LA,file = "rep2flx.LA")

rep2flx.LA <- readRDS("rep2flx.LA")

#Randomization
rep2flx.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep2flx.LA.rand.VCV) <- dimnames(rep2flx.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  FLX2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(FLX2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2flx.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped4,data=FLX2.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2flx.LA.rand.VCV[i,] <- posterior.mode(rep2flx.LA.rand$VCV)
  write.csv(rep2flx.LA.rand.VCV, "rep2flx.LA.rand.VCV.csv")
}

rep2flx.LA.rand.VCV <- read.csv("rep2flx.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep2flx.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep2flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep2flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep2flx.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep2flx.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep2flx.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep2flx.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep2flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep2flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep2flx.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep2flx.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep2flx.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep2flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep2flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep2flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep2flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep2flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep2flx.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep2flx.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep2flx.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep2flx.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep2flx.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep2flx.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)

#rep2CFM
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

#str(CFM2)
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cfm.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cfm.WS,file = "rep2cfm.WS")

rep2cfm.WS <- readRDS("rep2cfm.WS")

#Randomization
rep2cfm.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep2cfm.WS.rand.VCV) <- dimnames(rep2cfm.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cfm.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cfm.WS.rand.VCV[i,] <- posterior.mode(rep2cfm.WS.rand$VCV)
  write.csv(rep2cfm.WS.rand.VCV, "rep2cfm.WS.rand.VCV.csv")
}

rep2cfm.WS.rand.VCV <- read.csv("rep2cfm.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep2cfm.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep2cfm.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep2cfm.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep2cfm.WS$VCV[,"traitmWS.units"])
HPDinterval(rep2cfm.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"]+
     rep2cfm.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep2cfm.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep2cfm.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep2cfm.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep2cfm.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep2cfm.WS$VCV[,"traitfWS.units"])
HPDinterval(rep2cfm.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"]+
     rep2cfm.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep2cfm.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep2cfm.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep2cfm.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep2cfm.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep2cfm.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep2cfm.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep2cfm.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep2cfm.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep2cfm.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep2cfm.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep2cfm.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep2cfm.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep2cfm.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cfm.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cfm.DR, file="rep2cfm.DR")

rep2cfm.DR <- readRDS("rep2cfm.DR")

#Randomization
rep2cfm.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep2cfm.DR.rand.VCV) <- dimnames(rep2cfm.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cfm.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cfm.DR.rand.VCV[i,] <- posterior.mode(rep2cfm.DR.rand$VCV)
  write.csv(rep2cfm.DR.rand.VCV, "rep2cfm.DR.rand.VCV.csv")
}

rep2cfm.DR.rand.VCV <- read.csv("rep2cfm.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep2cfm.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep2cfm.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep2cfm.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep2cfm.DR$VCV[,"traitmDR.units"])
HPDinterval(rep2cfm.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep2cfm.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep2cfm.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep2cfm.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep2cfm.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep2cfm.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep2cfm.DR$VCV[,"traitfDR.units"])
HPDinterval(rep2cfm.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep2cfm.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep2cfm.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2cfm.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep2cfm.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep2cfm.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep2cfm.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep2cfm.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep2cfm.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep2cfm.DR$VCV[,"traitmDR:traitfDR.animal"]/
  + sqrt(rep2cfm.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep2cfm.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep2cfm.DR.rand.VCV$"traitmDR.traitfDR.animal"/
  + sqrt(rep2cfm.DR.rand.VCV$"traitmDR.traitmDR.animal"*
           rep2cfm.DR.rand.VCV$"traitfDR.traitfDR.animal")
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cfm.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cfm.LA,file = "rep2cfm.LA")

rep2cfm.LA <- readRDS("rep2cfm.LA")

#Randomization
rep2cfm.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep2cfm.LA.rand.VCV) <- dimnames(rep2cfm.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CFM2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CFM2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cfm.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped5,data=CFM2.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cfm.LA.rand.VCV[i,] <- posterior.mode(rep2cfm.LA.rand$VCV)
  write.csv(rep2cfm.LA.rand.VCV, "rep2cfm.LA.rand.VCV.csv")
}

rep2cfm.LA.rand.VCV <- read.csv("rep2cfm.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep2cfm.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep2cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep2cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep2cfm.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep2cfm.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep2cfm.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep2cfm.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep2cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep2cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep2cfm.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep2cfm.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep2cfm.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep2cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep2cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep2cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep2cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep2cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep2cfm.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep2cfm.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep2cfm.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep2cfm.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep2cfm.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep2cfm.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)

#rep2CWT
#WS
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cwt.WS<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cwt.WS,file = "rep2cwt.WS")

rep2cwt.WS <- readRDS("rep2cwt.WS")

#Randomization
rep2cwt.WS.rand.VCV <- matrix(NA,200,10)
colnames(rep2cwt.WS.rand.VCV) <- dimnames(rep2cwt.WS$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2.rand)
  
  #Analysis
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cwt.WS.rand<-MCMCglmm(cbind(mWS,fWS)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cwt.WS.rand.VCV[i,] <- posterior.mode(rep2cwt.WS.rand$VCV)
  write.csv(rep2cwt.WS.rand.VCV, "rep2cwt.WS.rand.VCV.csv")
}

rep2cwt.WS.rand.VCV <- read.csv("rep2cwt.WS.rand.VCV.csv")

#Significance testing
#Male Va wing size
posterior.mode(rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"])
HPDinterval(rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.WS.rand.VCV[,"traitmWS.traitmWS.animal"]))
HPDinterval(as.mcmc(rep2cwt.WS.rand.VCV[,"traitmWS.traitmWS.animal"],0.95))
sum(rep2cwt.WS.rand.VCV$"traitmWS.traitmWS.animal">=posterior.mode(rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"]))
hist(rep2cwt.WS.rand.VCV$"traitmWS.traitmWS.animal")
#Male Vr wing size
posterior.mode(rep2cwt.WS$VCV[,"traitmWS.units"])
HPDinterval(rep2cwt.WS$VCV[,"traitmWS.units"],0.95)
#Male heritability WS
h2.mWS <- rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"]/
  (rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"]+
     rep2cwt.WS$VCV[,"traitmWS:traitmWS.round.sec"]+
     rep2cwt.WS$VCV[,"traitmWS.units"])
posterior.mode(h2.mWS)
HPDinterval(h2.mWS,0.95)
#Female Va wing size
posterior.mode(rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"])
HPDinterval(rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.WS.rand.VCV[,"traitfWS.traitfWS.animal"]))
HPDinterval(as.mcmc(rep2cwt.WS.rand.VCV[,"traitfWS.traitfWS.animal"],0.95))
sum(rep2cwt.WS.rand.VCV$"traitfWS.traitfWS.animal">=posterior.mode(rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"]))
hist(rep2cwt.WS.rand.VCV$"traitfWS.traitfWS.animal")
#Female Vr WS
posterior.mode(rep2cwt.WS$VCV[,"traitfWS.units"])
HPDinterval(rep2cwt.WS$VCV[,"traitfWS.units"],0.95)
#Female heritability WS
h2.fWS <- rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"]/
  (rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"]+
     rep2cwt.WS$VCV[,"traitfWS:traitfWS.round.sec"]+
     rep2cwt.WS$VCV[,"traitfWS.units"])
posterior.mode(h2.fWS)
HPDinterval(h2.fWS,0.95)
#Intersexual genetic covariance WS
posterior.mode(rep2cwt.WS$VCV[,"traitmWS:traitfWS.animal"])
HPDinterval(rep2cwt.WS$VCV[,"traitmWS:traitfWS.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.WS.rand.VCV$"traitmWS.traitfWS.animal"))
HPDinterval(as.mcmc(rep2cwt.WS.rand.VCV$"traitmWS.traitfWS.animal",0.95))
sum(abs(rep2cwt.WS.rand.VCV$"traitmWS.traitfWS.animal")>=abs(posterior.mode(rep2cwt.WS$VCV[,"traitmWS:traitfWS.animal"])))
hist(rep2cwt.WS.rand.VCV$"traitmWS.traitfWS.animal")
#Intersexual genetic correlation WS
rmf.WS<-rep2cwt.WS$VCV[,"traitmWS:traitfWS.animal"]/
  + sqrt(rep2cwt.WS$VCV[,"traitmWS:traitmWS.animal"]*
           rep2cwt.WS$VCV[,"traitfWS:traitfWS.animal"])
posterior.mode(rmf.WS)
HPDinterval(rmf.WS,0.95)
rmf.WS.rand<-rep2cwt.WS.rand.VCV$"traitmWS.traitfWS.animal"/
  + sqrt(rep2cwt.WS.rand.VCV$"traitmWS.traitmWS.animal"*
           rep2cwt.WS.rand.VCV$"traitfWS.traitfWS.animal")
sum(abs(rmf.WS.rand)>=abs(posterior.mode(rmf.WS)))
hist(rmf.WS.rand)

#DR
#Analysis
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")

p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cwt.DR<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2,prior=p2,verbose=FALSE,
                     family=c("gaussian","gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cwt.DR, file="rep2cwt.DR")

rep2cwt.DR <- readRDS("rep2cwt.DR")

#Randomization
rep2cwt.DR.rand.VCV <- matrix(NA,200,10)
colnames(rep2cwt.DR.rand.VCV) <- dimnames(rep2cwt.DR$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cwt.DR.rand<-MCMCglmm(cbind(mDR,fDR)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2.rand,prior=p2,verbose=FALSE,
                            family=c("gaussian","gaussian"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cwt.DR.rand.VCV[i,] <- posterior.mode(rep2cwt.DR.rand$VCV)
  write.csv(rep2cwt.DR.rand.VCV, "rep2cwt.DR.rand.VCV.csv")
}

rep2cwt.DR.rand.VCV <- read.csv("rep2cwt.DR.rand.VCV.csv")

#Significance testing
#Male Va DR
posterior.mode(rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"])
HPDinterval(rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.DR.rand.VCV[,"traitmDR.traitmDR.animal"]))
HPDinterval(as.mcmc(rep2cwt.DR.rand.VCV[,"traitmDR.traitmDR.animal"],0.95))
sum(rep2cwt.DR.rand.VCV$"traitmDR.traitmDR.animal">=posterior.mode(rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"]))
hist(rep2cwt.DR.rand.VCV$"traitmDR.traitmDR.animal")
#Male Vr DR
posterior.mode(rep2cwt.DR$VCV[,"traitmDR.units"])
HPDinterval(rep2cwt.DR$VCV[,"traitmDR.units"],0.95)
#Male heritability DR
h2.mDR <- rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"]/
  (rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"]+
     rep2cwt.DR$VCV[,"traitmDR:traitmDR.round.sec"]+
     rep2cwt.DR$VCV[,"traitmDR.units"])
posterior.mode(h2.mDR)
HPDinterval(h2.mDR,0.95)
#Female Va DR
posterior.mode(rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"])
HPDinterval(rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.DR.rand.VCV[,"traitfDR.traitfDR.animal"]))
HPDinterval(as.mcmc(rep2cwt.DR.rand.VCV[,"traitfDR.traitfDR.animal"],0.95))
sum(rep2cwt.DR.rand.VCV$"traitfDR.traitfDR.animal">=posterior.mode(rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"]))
hist(rep2cwt.DR.rand.VCV$"traitfDR.traitfDR.animal")
#Female Vr DR
posterior.mode(rep2cwt.DR$VCV[,"traitfDR.units"])
HPDinterval(rep2cwt.DR$VCV[,"traitfDR.units"],0.95)
#Female heritability DR
h2.fDR <- rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"]/
  (rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"]+
     rep2cwt.DR$VCV[,"traitfDR:traitfDR.round.sec"]+
     rep2cwt.DR$VCV[,"traitfDR.units"])
posterior.mode(h2.fDR)
HPDinterval(h2.fDR,0.95)
#Intersexual genetic covariance DR
posterior.mode(rep2cwt.DR$VCV[,"traitmDR:traitfDR.animal"])
HPDinterval(rep2cwt.DR$VCV[,"traitmDR:traitfDR.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.DR.rand.VCV$"traitmDR.traitfDR.animal"))
HPDinterval(as.mcmc(rep2cwt.DR.rand.VCV$"traitmDR.traitfDR.animal",0.95))
sum(abs(rep2cwt.DR.rand.VCV$"traitmDR.traitfDR.animal")>=abs(posterior.mode(rep2cwt.DR$VCV[,"traitmDR:traitfDR.animal"])))
hist(rep2cwt.DR.rand.VCV$"traitmDR.traitfDR.animal")
#Intersexual genetic correlation DR
rmf.DR<-rep2cwt.DR$VCV[,"traitmDR:traitfDR.animal"]/
  + sqrt(rep2cwt.DR$VCV[,"traitmDR:traitmDR.animal"]*
           rep2cwt.DR$VCV[,"traitfDR:traitfDR.animal"])
posterior.mode(rmf.DR)
HPDinterval(rmf.DR,0.95)
rmf.DR.rand<-rep2cwt.DR.rand.VCV$"traitmDR.traitfDR.animal"/
  + sqrt(rep2cwt.DR.rand.VCV$"traitmDR.traitmDR.animal"*
           rep2cwt.DR.rand.VCV$"traitfDR.traitfDR.animal")
sum(abs(rmf.DR.rand)>=abs(posterior.mode(rmf.DR)))
hist(rmf.DR.rand)

#Locomotion
p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
rep2cwt.LA<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                     random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2,prior=p2,verbose=FALSE,
                     family=c("poisson","poisson"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("poisson","poisson"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cwt.LA,file = "rep2cwt.LA")

rep2cwt.LA <- readRDS("rep2cwt.LA")

#Randomization
rep2cwt.LA.rand.VCV <- matrix(NA,200,10)
colnames(rep2cwt.LA.rand.VCV) <- dimnames(rep2cwt.LA$VCV)[[2]]
for (i in 1:200) {   
  print(i)
  #4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
  #Round 1 females
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  #str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$mWS <- rep(NA, nrow(r1F))
  r1F.rand$mDR <- rep(NA, nrow(r1F))
  r1F.rand$mLA.active <- rep(NA, nrow(r1F))
  r1F.rand$mLA.passive <- rep(NA, nrow(r1F))
  r1F.rand$fWS <- r1F[rand.order,20]
  r1F.rand$fDR <- r1F[rand.order,22]
  r1F.rand$fLA.active<- r1F[rand.order,17]
  r1F.rand$fLA.passive<- r1F[rand.order,18]
  #Round 2 females
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  #str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$mWS <- rep(NA, nrow(r2F))
  r2F.rand$mDR <- rep(NA, nrow(r2F))
  r2F.rand$mLA.active <- rep(NA, nrow(r2F))
  r2F.rand$mLA.passive <- rep(NA, nrow(r2F))
  r2F.rand$fWS <- r2F[rand.order,20]
  r2F.rand$fDR <- r2F[rand.order,22]
  r2F.rand$fLA.active<- r2F[rand.order,17]
  r2F.rand$fLA.passive<- r2F[rand.order,18]
  #str(r2F.rand)
  #Round 3 females
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  #str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$mWS <- rep(NA, nrow(r3F))
  r3F.rand$mDR <- rep(NA, nrow(r3F))
  r3F.rand$mLA.active <- rep(NA, nrow(r3F))
  r3F.rand$mLA.passive <- rep(NA, nrow(r3F))
  r3F.rand$fWS <- r3F[rand.order,20]
  r3F.rand$fDR <- r3F[rand.order,22]
  r3F.rand$fLA.active<- r3F[rand.order,17]
  r3F.rand$fLA.passive<- r3F[rand.order,18]
  #str(r3F.rand)
  #Round 4 females
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  #str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$mWS <- rep(NA, nrow(r4F))
  r4F.rand$mDR <- rep(NA, nrow(r4F))
  r4F.rand$mLA.active <- rep(NA, nrow(r4F))
  r4F.rand$mLA.passive <- rep(NA, nrow(r4F))
  r4F.rand$fWS <- r4F[rand.order,20]
  r4F.rand$fDR <- r4F[rand.order,22]
  r4F.rand$fLA.active<- r4F[rand.order,17]
  r4F.rand$fLA.passive<- r4F[rand.order,18]
  #str(r4F.rand)
  #Round 1 males
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  #str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$mWS <- r1M[rand.order,19]
  r1M.rand$mDR <- r1M[rand.order,21]
  r1M.rand$mLA.active <- r1M[rand.order,15]
  r1M.rand$mLA.passive <- r1M[rand.order,16]
  r1M.rand$fWS <- rep(NA, nrow(r1M))
  r1M.rand$fDR <- rep(NA, nrow(r1M))
  r1M.rand$fLA.active <- rep(NA, nrow(r1M))
  r1M.rand$fLA.passive <- rep(NA, nrow(r1M))
  #Round 2 males
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  #str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$mWS <- r2M[rand.order,19]
  r2M.rand$mDR <- r2M[rand.order,21]
  r2M.rand$mLA.active <- r2M[rand.order,15]
  r2M.rand$mLA.passive <- r2M[rand.order,16]
  r2M.rand$fWS <- rep(NA, nrow(r2M))
  r2M.rand$fDR <- rep(NA, nrow(r2M))
  r2M.rand$fLA.active <- rep(NA, nrow(r2M))
  r2M.rand$fLA.passive <- rep(NA, nrow(r2M))
  #Round 3 males
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  #str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$mWS <- r3M[rand.order,19]
  r3M.rand$mDR <- r3M[rand.order,21]
  r3M.rand$mLA.active <- r3M[rand.order,15]
  r3M.rand$mLA.passive <- r3M[rand.order,16]
  r3M.rand$fWS <- rep(NA, nrow(r3M))
  r3M.rand$fDR <- rep(NA, nrow(r3M))
  r3M.rand$fLA.active <- rep(NA, nrow(r3M))
  r3M.rand$fLA.passive <- rep(NA, nrow(r3M))
  #Round 4 males
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  #str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$mWS <- r4M[rand.order,19]
  r4M.rand$mDR <- r4M[rand.order,21]
  r4M.rand$mLA.active <- r4M[rand.order,15]
  r4M.rand$mLA.passive <- r4M[rand.order,16]
  r4M.rand$fWS <- rep(NA, nrow(r4M))
  r4M.rand$fDR <- rep(NA, nrow(r4M))
  r4M.rand$fLA.active <- rep(NA, nrow(r4M))
  r4M.rand$fLA.passive <- rep(NA, nrow(r4M))
  #Combine randomized subsets
  CWT2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  #str(CWT2.rand)
  
  #Analysis randomized dataset
  p2<-list(G=list(G1=list(V=diag(2),nu=0.002),G2=list(V=diag(2),nu=0.002)),R=list(V=diag(2),nu=0.002))
  rep2cwt.LA.rand<-MCMCglmm(cbind(mLA.active,fLA.active)~trait-1+trait:round,
                            random=~us(trait):animal+us(trait):round.sec,rcov=~idh(trait):units,pedigree=ped6,data=CWT2.rand,prior=p2,verbose=FALSE,
                            family=c("poisson","poisson"),nitt=550000,thin=2000,burnin=150000) 
  #For shorter run time:
  #family=c("gaussian","gaussian"),nitt=1000,thin=10,burnin=100) 
  rep2cwt.LA.rand.VCV[i,] <- posterior.mode(rep2cwt.LA.rand$VCV)
  write.csv(rep2cwt.LA.rand.VCV, "rep2cwt.LA.rand.VCV.csv")
}

rep2cwt.LA.rand.VCV <- read.csv("rep2cwt.LA.rand.VCV.csv")

#Significance testing
#Male Va LA
posterior.mode(rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"])
HPDinterval(rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"]))
HPDinterval(as.mcmc(rep2cwt.LA.rand.VCV[,"traitmLA.active.traitmLA.active.animal"],0.95))
sum(rep2cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal">=posterior.mode(rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]))
hist(rep2cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal")
#Male Vr LA
posterior.mode(rep2cwt.LA$VCV[,"traitmLA.active.units"])
HPDinterval(rep2cwt.LA$VCV[,"traitmLA.active.units"],0.95)
#Male heritability LA
h2.mLA <- rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]/
  (rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]+
     rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.round.sec"]+
     rep2cwt.LA$VCV[,"traitmLA.active.units"])
posterior.mode(h2.mLA)
HPDinterval(h2.mLA,0.95)
#Female Va LA
posterior.mode(rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
HPDinterval(rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"]))
HPDinterval(as.mcmc(rep2cwt.LA.rand.VCV[,"traitfLA.active.traitfLA.active.animal"],0.95))
sum(rep2cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal">=posterior.mode(rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]))
hist(rep2cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
#Female Vr LA
posterior.mode(rep2cwt.LA$VCV[,"traitfLA.active.units"])
HPDinterval(rep2cwt.LA$VCV[,"traitfLA.active.units"],0.95)
#Female heritability LA
h2.fLA <- rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]/
  (rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"]+
     rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.round.sec"]+
     rep2cwt.LA$VCV[,"traitfLA.active.units"])
posterior.mode(h2.fLA)
HPDinterval(h2.fLA,0.95)
#Intersexual genetic covariance LA
posterior.mode(rep2cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])
HPDinterval(rep2cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"],0.95)
posterior.mode(as.mcmc(rep2cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"))
HPDinterval(as.mcmc(rep2cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal",0.95))
sum(abs(rep2cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")>=abs(posterior.mode(rep2cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"])))
hist(rep2cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal")
#Intersexual genetic correlation LA
rmf.LA<-rep2cwt.LA$VCV[,"traitmLA.active:traitfLA.active.animal"]/
  + sqrt(rep2cwt.LA$VCV[,"traitmLA.active:traitmLA.active.animal"]*
           rep2cwt.LA$VCV[,"traitfLA.active:traitfLA.active.animal"])
posterior.mode(rmf.LA)
HPDinterval(rmf.LA,0.95)
rmf.LA.rand<-rep2cwt.LA.rand.VCV$"traitmLA.active.traitfLA.active.animal"/
  + sqrt(rep2cwt.LA.rand.VCV$"traitmLA.active.traitmLA.active.animal"*
           rep2cwt.LA.rand.VCV$"traitfLA.active.traitfLA.active.animal")
sum(abs(rmf.LA.rand)>=abs(posterior.mode(rmf.LA)))
hist(rmf.LA.rand)
