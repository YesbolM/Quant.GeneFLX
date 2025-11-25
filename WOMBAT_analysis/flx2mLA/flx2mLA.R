#rep2flx mLA

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
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
flx2mLA<-flx2m[,c(1,9,13:14)]
str(flx2mLA)
#Make round.sec into integers by taking factor codes
flx2mLA$round.sec <- as.numeric(flx2mLA$round.sec)
#drop unused levels
flx2mLA<-droplevels(flx2mLA) 
str(flx2mLA)
#round: 4 levels
length(unique(flx2mLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2mLA$animal<-gsub("M","1",flx2mLA$animal)
flx2mLA$animal<-gsub("F","2",flx2mLA$animal)
flx2mLA$animal<-gsub("m","3",flx2mLA$animal)
flx2mLA$animal<-gsub("f","4",flx2mLA$animal)
flx2mLA$animal<-gsub("-","",flx2mLA$animal)

#Change name of round.sec to NSEC
names(flx2mLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mLA")
write.table(flx2mLA, file="flx2mLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2mLA")
write.table(ped.flx2, file="flx2mLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2FLXmLA <- readRDS("rep2FLX.mLA")
library(MCMCglmm)
posterior.mode(rep2FLXmLA$VCV)
#animal    round.sec        units 
#0.002459238 0.133552038 0.909779602

# Testing significance of random effects
#logL for model with animal: -783.144
#logL for model without animal: -783.163
q <- 2*(-783.144+783.163)
q
#q=0.038
pchisq(q, df=1, lower.tail = FALSE)
#p=0.8454431