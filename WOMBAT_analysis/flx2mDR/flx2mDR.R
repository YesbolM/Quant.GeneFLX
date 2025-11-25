#rep2flx mDR

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
flx2m<-flx2[flx2$sex=="male",]
str(flx2m)
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
flx2mDR<-flx2m[,c(1,8,13:14)]
str(flx2mDR)
#Make round.sec into integers by taking factor codes
flx2mDR$round.sec <- as.numeric(flx2mDR$round.sec)
#Drop unused levels
flx2mDR<-droplevels(flx2mDR) 
str(flx2mDR)
#round: 4 levels
length(unique(flx2mDR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2mDR$animal<-gsub("M","1",flx2mDR$animal)
flx2mDR$animal<-gsub("F","2",flx2mDR$animal)
flx2mDR$animal<-gsub("m","3",flx2mDR$animal)
flx2mDR$animal<-gsub("f","4",flx2mDR$animal)
flx2mDR$animal<-gsub("-","",flx2mDR$animal)

#Change name of round.sec to NSEC
names(flx2mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mDR")
write.table(flx2mDR, file="flx2mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mDR")
write.table(ped.flx2, file="flx2mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2FLXmDR <- readRDS("rep2FLX.mDR")
library(MCMCglmm)
posterior.mode(rep2FLXmDR$VCV)
#animal    round.sec        units 
#1.180969  1.127906  2.048701 

# Testing significance of random effects
#logL for model with animal: -536.386
#logL for model without animal: -543.628
q <- 2*(-536.386+543.628)
q
#q=14.484
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0001413551