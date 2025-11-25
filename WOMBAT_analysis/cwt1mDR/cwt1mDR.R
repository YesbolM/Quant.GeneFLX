#rep1cwt mDR

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
cwt1mDR<-cwt1m[,c(1,8,13:14)]
str(cwt1mDR)
#Make round.sec into integers by taking factor codes
cwt1mDR$round.sec <- as.numeric(cwt1mDR$round.sec)
#Drop unused levels
cwt1mDR<-droplevels(cwt1mDR) 
str(cwt1mDR)
#round: 4 levels
length(unique(cwt1mDR$round.sec))
#round.sec: 28 levels

#Change ids into numeric values
cwt1mDR$animal<-gsub("M","1",cwt1mDR$animal)
cwt1mDR$animal<-gsub("F","2",cwt1mDR$animal)
cwt1mDR$animal<-gsub("m","3",cwt1mDR$animal)
cwt1mDR$animal<-gsub("f","4",cwt1mDR$animal)
cwt1mDR$animal<-gsub("-","",cwt1mDR$animal)

#Change name of round.sec to NSEC
names(cwt1mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1mDR")
write.table(cwt1mDR, file="cwt1mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1mDR")
write.table(ped.cwt1, file="cwt1mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CWTmDR <- readRDS("rep1CWT.mDR")
library(MCMCglmm)
posterior.mode(rep1CWTmDR$VCV)
#animal    round.sec        units 
#2.4570931 0.9757411 3.1689966

# Testing significance of random effects
#logL for model with animal: -697.934
#logL for model without animal: -704.984
q <- 2*(-697.934+704.984)
q
#q=14.1
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0001733438
