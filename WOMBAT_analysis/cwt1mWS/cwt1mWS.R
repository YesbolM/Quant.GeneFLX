#rep1cwt mWS

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
cwt1m<-cwt1[cwt1$sex=="male",]
str(cwt1m)
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cwt1mWS<-cwt1m[,c(1,7,13:14)]
str(cwt1mWS)
#Make round.sec into integers by taking factor codes
cwt1mWS$round.sec <- as.numeric(cwt1mWS$round.sec)
#Drop unused levels
cwt1mWS<-droplevels(cwt1mWS) 
str(cwt1mWS)
#round: 4 levels
length(unique(cwt1mWS$round.sec))
#round.sec: 28 levels

#Change ids into numeric values
cwt1mWS$animal<-gsub("M","1",cwt1mWS$animal)
cwt1mWS$animal<-gsub("F","2",cwt1mWS$animal)
cwt1mWS$animal<-gsub("m","3",cwt1mWS$animal)
cwt1mWS$animal<-gsub("f","4",cwt1mWS$animal)
cwt1mWS$animal<-gsub("-","",cwt1mWS$animal)

#Change name of round.sec to NSEC
names(cwt1mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1mWS")
write.table(cwt1mWS, file="cwt1mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt1 <- read.csv("rep1CWT.p.csv")
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1mWS")
write.table(ped.cwt1, file="cwt1mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CWTmWS <- readRDS("rep1CWT.mWS")
library(MCMCglmm)
posterior.mode(rep1CWTmWS$VCV)
#animal    round.sec        units 
#0.0019850273 0.0003871213 0.0006190984 

# Testing significance of random effects
#logL for model with animal: 1273.325
#logL for model without animal: 1250.782
q <- 2*(1273.325-1250.782)
q
#q=45.086
pchisq(q, df=1, lower.tail = FALSE)
#p=1.885687e-11

