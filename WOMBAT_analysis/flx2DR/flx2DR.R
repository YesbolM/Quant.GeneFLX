#rep2flx DR

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
flx2<-rep2[rep2$treatment=="FLX",] 
#  WOMBAT does not allow any missing values, let's remove them
flx2<-na.omit(flx2)
#Select wing size for analysis
flx2DR<-flx2[c(8,1,5,13:14)]
flx2DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx2DR$TRAITNO<-ifelse(flx2DR$sex=="male",1,2)
flx2DR$sex<-ifelse(flx2DR$sex=="male",1,2)
#Reorder data
flx2DR<-flx2DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx2DR$round.sec <- as.numeric(flx2DR$round.sec)
#drop unused levels
flx2DR<-droplevels(flx2DR) 
str(flx2DR)
#round: 4 levels
length(unique(flx2DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2DR$animal<-gsub("M","1",flx2DR$animal)
flx2DR$animal<-gsub("F","2",flx2DR$animal)
flx2DR$animal<-gsub("m","3",flx2DR$animal)
flx2DR$animal<-gsub("f","4",flx2DR$animal)
flx2DR$animal<-gsub("-","",flx2DR$animal)
#Sort by animal ID
flx2DR <- flx2DR[order(flx2DR$animal),]

#Change name of round.sec to NSEC
names(flx2DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2DR")
write.table(flx2DR, file="flx2DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.flx2 <- read.csv("rep2flx.p.csv")
head(ped.flx2)
str(ped.flx2)
#Change ids into numeric values
ped.flx2$id<-gsub("M","1",ped.flx2$id)
ped.flx2$id<-gsub("F","2",ped.flx2$id)
ped.flx2$id<-gsub("m","3",ped.flx2$id)
ped.flx2$id<-gsub("f","4",ped.flx2$id)
ped.flx2$id<-gsub("-","",ped.flx2$id)
ped.flx2$FATHER<-gsub("M","1",ped.flx2$FATHER)
ped.flx2$FATHER[is.na(ped.flx2$FATHER)] <- 0
ped.flx2$MOTHER<-gsub("F","2",ped.flx2$MOTHER)
ped.flx2$MOTHER[is.na(ped.flx2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2DR")
write.table(ped.flx2, file="flx2DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2flx.DR <- readRDS("rep2flx.DR")
posterior.mode(rep2flx.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal    traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#0.7238784                   2.7116757                   2.7116757                   9.5338945 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#1.2737816                   1.1064754                   1.1064754                   0.9316502 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units     traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#2.3903748                   4.9243025                   4.9243025                  10.8858531 

#Approximate CIs
#Va mDR
1.01878-(1.96*0.383198)
1.01878+(1.96*0.383198)
#Va fDR
8.85594-(1.96*2.86414)
8.85594+(1.96*2.86414)
#Cov DR
2.63559-(1.96*0.746500)
2.63559+(1.96*0.746500)
#rmf DR
0.877-(1.96*0.207)
0.877+(1.96*0.207)

#Log-likelihood test for covariance
#logL for model with cov: -1454.824
#logL for model without cov: -1462.654
q <- 2*(-1454.824+1462.654)
q
#q=15.66
pchisq(q, df=1, lower.tail = FALSE)
#p=7.581066e-05
