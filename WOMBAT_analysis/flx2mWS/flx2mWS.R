#rep2flx mWS

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
#Select WS for analysis
#We need to include animal, WS, round, and round.sec
flx2mWS<-flx2m[,c(1,7,13:14)]
str(flx2mWS)
#Make round.sec into integers by taking factor codes
flx2mWS$round.sec <- as.numeric(flx2mWS$round.sec)
#drop unused levels
flx2mWS<-droplevels(flx2mWS) 
str(flx2mWS)
#round: 4 levels
length(unique(flx2mWS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2mWS$animal<-gsub("M","1",flx2mWS$animal)
flx2mWS$animal<-gsub("F","2",flx2mWS$animal)
flx2mWS$animal<-gsub("m","3",flx2mWS$animal)
flx2mWS$animal<-gsub("f","4",flx2mWS$animal)
flx2mWS$animal<-gsub("-","",flx2mWS$animal)

#Change name of round.sec to NSEC
names(flx2mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mWS")
write.table(flx2mWS, file="flx2mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mWS")
write.table(ped.flx2, file="flx2mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2FLXmWS <- readRDS("rep2FLX.mWS")
library(MCMCglmm)
posterior.mode(rep2FLXmWS$VCV)
#animal    round.sec        units 
#0.0021260612 0.0002274540 0.0005141685

# Testing significance of random effects
#logL for model with animal: 1238.072
#logL for model without animal: 1218.211
q <- 2*(1238.072-1218.211)
q
#q=39.722
pchisq(q, df=1, lower.tail = FALSE)
#p=2.928085e-10
