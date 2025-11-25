#rep1cfm fLA

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
cfm1f<-cfm1[cfm1$sex=="female",]
str(cfm1f)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cfm1fLA<-cfm1f[,c(1,9,13:14)]
str(cfm1fLA)
#Make round.sec into integers by taking factor codes
cfm1fLA$round.sec <- as.numeric(cfm1fLA$round.sec)
#drop unused levels
cfm1fLA<-droplevels(cfm1fLA) 
str(cfm1fLA)
#round: 4 levels
length(unique(cfm1fLA$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm1fLA$animal<-gsub("M","1",cfm1fLA$animal)
cfm1fLA$animal<-gsub("F","2",cfm1fLA$animal)
cfm1fLA$animal<-gsub("m","3",cfm1fLA$animal)
cfm1fLA$animal<-gsub("f","4",cfm1fLA$animal)
cfm1fLA$animal<-gsub("-","",cfm1fLA$animal)

#Change name of round.sec to NSEC
names(cfm1fLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1fLA")
write.table(cfm1fLA, file="cfm1fLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1fLA")
write.table(ped.cfm1, file="cfm1fLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CFMfLA <- readRDS("rep1CFM.fLA")
library(MCMCglmm)
posterior.mode(rep1CFMfLA$VCV)
#animal    round.sec        units 
#0.009311440 0.003077282 1.220311464

# Testing significance of random effects
#logL for model with animal: -589.624
#logL for model without animal: -592.951
q <- 2*(-589.624+592.951)
q
#q=6.654
pchisq(q, df=1, lower.tail = FALSE)
#p=0.009893348