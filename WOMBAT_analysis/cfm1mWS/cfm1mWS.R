#rep1cfm mWS

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
cfm1<-rep1[rep1$treatment=="CFM",] 
#  WOMBAT does not allow any missing values, let's remove them
cfm1<-na.omit(cfm1) 
cfm1m<-cfm1[cfm1$sex=="male",]
str(cfm1m)
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cfm1mWS<-cfm1m[,c(1,7,13:14)]
str(cfm1mWS)
#Make round.sec into integers by taking factor codes
cfm1mWS$round.sec <- as.numeric(cfm1mWS$round.sec)
#Drop unused levels
cfm1mWS<-droplevels(cfm1mWS) 
str(cfm1mWS)
#round: 4 levels
length(unique(cfm1mWS$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm1mWS$animal<-gsub("M","1",cfm1mWS$animal)
cfm1mWS$animal<-gsub("F","2",cfm1mWS$animal)
cfm1mWS$animal<-gsub("m","3",cfm1mWS$animal)
cfm1mWS$animal<-gsub("f","4",cfm1mWS$animal)
cfm1mWS$animal<-gsub("-","",cfm1mWS$animal)

#Change name of round.sec to NSEC
names(cfm1mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mWS")
write.table(cfm1mWS, file="cfm1mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cfm1 <- read.csv("rep1cfm.p.csv")
head(ped.cfm1)
str(ped.cfm1)
#Change ids into numeric values
ped.cfm1$id<-gsub("M","1",ped.cfm1$id)
ped.cfm1$id<-gsub("F","2",ped.cfm1$id)
ped.cfm1$id<-gsub("m","3",ped.cfm1$id)
ped.cfm1$id<-gsub("f","4",ped.cfm1$id)
ped.cfm1$id<-gsub("-","",ped.cfm1$id)
ped.cfm1$FATHER<-gsub("M","1",ped.cfm1$FATHER)
ped.cfm1$FATHER[is.na(ped.cfm1$FATHER)] <- 0
ped.cfm1$MOTHER<-gsub("F","2",ped.cfm1$MOTHER)
ped.cfm1$MOTHER[is.na(ped.cfm1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mWS")
write.table(ped.cfm1, file="cfm1mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CFMmWS <- readRDS("rep1CFM.mWS")
library(MCMCglmm)
posterior.mode(rep1CFMmWS$VCV)
#animal    round.sec        units 
#0.0018875054 0.0001875669 0.0007920928 

# Testing significance of random effects
#logL for model with animal: 1220.583
#logL for model without animal: 1204.850
q <- 2*(1220.583-1204.850)
q
#q= 31.466
pchisq(q, df=1, lower.tail = FALSE)
#p=2.029632e-08