#rep1flx mWS

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
flx1<-rep1[rep1$treatment=="FLX",] 
#  WOMBAT does not allow any missing values, let's remove them
flx1<-na.omit(flx1)
flx1m<-flx1[flx1$sex=="male",]
str(flx1m)
#Select WS for analysis
#We need to include animal, WS, round, and round.sec
flx1mWS<-flx1m[,c(1,7,13:14)]
str(flx1mWS)
#Make round.sec into integers by taking factor codes
flx1mWS$round.sec <- as.numeric(flx1mWS$round.sec)
#drop unused levels
flx1mWS<-droplevels(flx1mWS) 
str(flx1mWS)
#round: 4 levels
length(unique(flx1mWS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1mWS$animal<-gsub("M","1",flx1mWS$animal)
flx1mWS$animal<-gsub("F","2",flx1mWS$animal)
flx1mWS$animal<-gsub("m","3",flx1mWS$animal)
flx1mWS$animal<-gsub("f","4",flx1mWS$animal)
flx1mWS$animal<-gsub("-","",flx1mWS$animal)

#Change name of round.sec to NSEC
names(flx1mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mWS")
write.table(flx1mWS, file="flx1mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1mWS_null <- flx1mWS[order(flx1mWS$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mWS\\null model")
write.table(flx1mWS_null, file="flx1mWS_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.flx1 <- read.csv("rep1flx.p.csv")
head(ped.flx1)
str(ped.flx1)
#Change ids into numeric values
ped.flx1$id<-gsub("M","1",ped.flx1$id)
ped.flx1$id<-gsub("F","2",ped.flx1$id)
ped.flx1$id<-gsub("m","3",ped.flx1$id)
ped.flx1$id<-gsub("f","4",ped.flx1$id)
ped.flx1$id<-gsub("-","",ped.flx1$id)
ped.flx1$FATHER<-gsub("M","1",ped.flx1$FATHER)
ped.flx1$FATHER[is.na(ped.flx1$FATHER)] <- 0
ped.flx1$MOTHER<-gsub("F","2",ped.flx1$MOTHER)
ped.flx1$MOTHER[is.na(ped.flx1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mWS")
write.table(ped.flx1, file="flx1mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXmWS <- readRDS("rep1FLX.mWS")
library(MCMCglmm)
posterior.mode(rep1FLXmWS$VCV)
#animal    round.sec        units 
#0.0016804949 0.0002322018 0.0006316800 

# Testing significance of random effects
#logL for model with animal: 1246.059
#logL for model without animal: 1228.942
q <- 2*(1246.059-1228.942)
q
#q=34.234
pchisq(q, df=1, lower.tail = FALSE)
#p=4.886767e-09