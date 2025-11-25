#rep2cwt DR

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
data <- read.csv("data.csv")
str(data)
#Rename id column
names(data)[1]<-"animal"
#Make things that should be factors into factors
for(x in c(1:5,12:14))data[,x]<-as.factor(data[,x])

#Subset data
rep2<-data[data$replicat=="rep2",]
cwt2<-rep2[rep2$treatment=="CWT",] 
#  WOMBAT does not allow any missing values, let's remove them
cwt2<-na.omit(cwt2)
#Select wing size for analysis
cwt2DR<-cwt2[c(8,1,5,13:14)]
cwt2DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt2DR$TRAITNO<-ifelse(cwt2DR$sex=="male",1,2)
cwt2DR$sex<-ifelse(cwt2DR$sex=="male",1,2)
#Reorder data
cwt2DR<-cwt2DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt2DR$round.sec <- as.numeric(cwt2DR$round.sec)
#drop unused levels
cwt2DR<-droplevels(cwt2DR) 
str(cwt2DR)
#round: 4 levels
length(unique(cwt2DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2DR$animal<-gsub("M","1",cwt2DR$animal)
cwt2DR$animal<-gsub("F","2",cwt2DR$animal)
cwt2DR$animal<-gsub("m","3",cwt2DR$animal)
cwt2DR$animal<-gsub("f","4",cwt2DR$animal)
cwt2DR$animal<-gsub("-","",cwt2DR$animal)
#Sort by animal ID
cwt2DR <- cwt2DR[order(cwt2DR$animal),]

#Change name of round.sec to NSEC
names(cwt2DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2DR")
write.table(cwt2DR, file="cwt2DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt2 <- read.csv("rep2cwt.p.csv")
head(ped.cwt2)
str(ped.cwt2)
#Change ids into numeric values
ped.cwt2$id<-gsub("M","1",ped.cwt2$id)
ped.cwt2$id<-gsub("F","2",ped.cwt2$id)
ped.cwt2$id<-gsub("m","3",ped.cwt2$id)
ped.cwt2$id<-gsub("f","4",ped.cwt2$id)
ped.cwt2$id<-gsub("-","",ped.cwt2$id)
ped.cwt2$FATHER<-gsub("M","1",ped.cwt2$FATHER)
ped.cwt2$FATHER[is.na(ped.cwt2$FATHER)] <- 0
ped.cwt2$MOTHER<-gsub("F","2",ped.cwt2$MOTHER)
ped.cwt2$MOTHER[is.na(ped.cwt2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2DR")
write.table(ped.cwt2, file="cwt2DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cwt.DR <- readRDS("rep2cwt.DR")
posterior.mode(rep2cwt.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal    traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#1.19553706                  3.76499962                  3.76499962                 10.07541403 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#0.35024075                 -0.17989448                 -0.17989448                  0.01787834 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units     traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#3.91545858                 -7.71844283                 -7.71844283                 13.44201673   

#Approximate CIs
#Va mDR
1.77835-(1.96*0.680988)
1.77835+(1.96*0.680988)
#Va fDR
13.2936-(1.96*4.09256)
13.2936+(1.96*4.09256)
#Cov DR
3.93426-(1.96*1.18862)
3.93426+(1.96*1.18862)
#rmf DR
0.809-(1.96*0.196)
0.809+(1.96*0.196)

#Log-likelihood test for covariance
#logL for model with cov: -1624.503
#logL for model without cov: -1630.628
q <- 2*(-1624.503+1630.628)
q
#q=12.25
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0004652582

