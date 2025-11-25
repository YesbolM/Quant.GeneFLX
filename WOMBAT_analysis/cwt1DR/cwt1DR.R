#rep1cwt DR

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
data <- read.csv("data.csv")
str(data)
#Rename id column
names(data)[1]<-"animal"
#Make things that should be factors into factors
for(x in c(1:5,12:14))data[,x]<-as.factor(data[,x])

#Subset data
rep1<-data[data$replicat=="rep1",]
cwt1<-rep1[rep1$treatment=="CWT",] 
#  WOMBAT does not allow any missing values, let's remove them
cwt1<-na.omit(cwt1)
#Select wing size for analysis
cwt1DR<-cwt1[c(8,1,5,13:14)]
cwt1DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt1DR$TRAITNO<-ifelse(cwt1DR$sex=="male",1,2)
cwt1DR$sex<-ifelse(cwt1DR$sex=="male",1,2)
#Reorder data
cwt1DR<-cwt1DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt1DR$round.sec <- as.numeric(cwt1DR$round.sec)
#drop unused levels
cwt1DR<-droplevels(cwt1DR) 
str(cwt1DR)
#round: 4 levels
length(unique(cwt1DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt1DR$animal<-gsub("M","1",cwt1DR$animal)
cwt1DR$animal<-gsub("F","2",cwt1DR$animal)
cwt1DR$animal<-gsub("m","3",cwt1DR$animal)
cwt1DR$animal<-gsub("f","4",cwt1DR$animal)
cwt1DR$animal<-gsub("-","",cwt1DR$animal)
#Sort by animal ID
cwt1DR <- cwt1DR[order(cwt1DR$animal),]

#Change name of round.sec to NSEC
names(cwt1DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1DR")
write.table(cwt1DR, file="cwt1DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt1 <- read.csv("rep1cwt.p.csv")
head(ped.cwt1)
str(ped.cwt1)
#Change ids into numeric values
ped.cwt1$id<-gsub("M","1",ped.cwt1$id)
ped.cwt1$id<-gsub("F","2",ped.cwt1$id)
ped.cwt1$id<-gsub("m","3",ped.cwt1$id)
ped.cwt1$id<-gsub("f","4",ped.cwt1$id)
ped.cwt1$id<-gsub("-","",ped.cwt1$id)
ped.cwt1$FATHER<-gsub("M","1",ped.cwt1$FATHER)
ped.cwt1$FATHER[is.na(ped.cwt1$FATHER)] <- 0
ped.cwt1$MOTHER<-gsub("F","2",ped.cwt1$MOTHER)
ped.cwt1$MOTHER[is.na(ped.cwt1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1DR")
write.table(ped.cwt1, file="cwt1DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cwt.DR <- readRDS("rep1cwt.DR")
posterior.mode(rep1cwt.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal    traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#2.73344448                  6.31373405                  6.31373405                 13.11918109 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#0.79261428                  1.01495792                  1.01495792                  0.05855447 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units     traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#2.56893112                  5.87713783                  5.87713783                 12.40441356 

#Approximate CIs
#Va mDR
2.71101-(1.96*0.826908)
2.71101+(1.96*0.826908)
#Va fDR
10.9746-(1.96*3.84066)
10.9746+(1.96*3.84066)
#Cov DR
5.45356-(1.96*1.33934)
5.45356+(1.96*1.33934)
#rmf DR
1.000-(1.96*0.166)
1.000+(1.96*0.166)

#Log-likelihood test for covariance
#logL for model with cov: -1512.496
#logL for model without cov: -1524.285
q <- 2*(-1512.496+1524.285)
q
#q=23.578
pchisq(q, df=1, lower.tail = FALSE)
#p=1.199511e-06

