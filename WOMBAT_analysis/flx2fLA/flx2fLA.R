#rep2flx fLA

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
flx2f<-flx2[flx2$sex=="female",]
str(flx2f)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
flx2fLA<-flx2f[,c(1,9,13:14)]
str(flx2fLA)
#Make round.sec into integers by taking factor codes
flx2fLA$round.sec <- as.numeric(flx2fLA$round.sec)
#drop unused levels
flx2fLA<-droplevels(flx2fLA) 
str(flx2fLA)
#round: 4 levels
length(unique(flx2fLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2fLA$animal<-gsub("M","1",flx2fLA$animal)
flx2fLA$animal<-gsub("F","2",flx2fLA$animal)
flx2fLA$animal<-gsub("m","3",flx2fLA$animal)
flx2fLA$animal<-gsub("f","4",flx2fLA$animal)
flx2fLA$animal<-gsub("-","",flx2fLA$animal)

#Change name of round.sec to NSEC
names(flx2fLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2fLA")
write.table(flx2fLA, file="flx2fLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2fLA")
write.table(ped.flx2, file="flx2fLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2FLXfLA <- readRDS("rep2FLX.fLA")
library(MCMCglmm)
posterior.mode(rep2FLXfLA$VCV)
#animal    round.sec        units 
#0.005146827 0.160110421 1.302812989 

# Testing significance of random effects
#logL for model with animal: -616.889
#logL for model without animal: -618.159
q <- 2*(-616.889+618.159)
q
#q=2.54
pchisq(q, df=1, lower.tail = FALSE)
#p=0.1109948
