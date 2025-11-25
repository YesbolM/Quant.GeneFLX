#rep2cfm mWS

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
cfm2<-rep2[rep2$treatment=="CFM",]
#  WOMBAT does not allow any missing values, let's remove them
cfm2<-na.omit(cfm2) 
cfm2m<-cfm2[cfm2$sex=="male",]
str(cfm2m)
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cfm2mWS<-cfm2m[,c(1,7,13:14)]
str(cfm2mWS)
#Make round.sec into integers by taking factor codes
cfm2mWS$round.sec <- as.numeric(cfm2mWS$round.sec)
#Drop unused levels
cfm2mWS<-droplevels(cfm2mWS) 
str(cfm2mWS)
#round: 4 levels
length(unique(cfm2mWS$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm2mWS$animal<-gsub("M","1",cfm2mWS$animal)
cfm2mWS$animal<-gsub("F","2",cfm2mWS$animal)
cfm2mWS$animal<-gsub("m","3",cfm2mWS$animal)
cfm2mWS$animal<-gsub("f","4",cfm2mWS$animal)
cfm2mWS$animal<-gsub("-","",cfm2mWS$animal)

#Change name of round.sec to NSEC
names(cfm2mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mWS")
write.table(cfm2mWS, file="cfm2mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cfm2 <- read.csv("rep2cfm.p.csv")
head(ped.cfm2)
str(ped.cfm2)
#Change ids into numeric values
ped.cfm2$id<-gsub("M","1",ped.cfm2$id)
ped.cfm2$id<-gsub("F","2",ped.cfm2$id)
ped.cfm2$id<-gsub("m","3",ped.cfm2$id)
ped.cfm2$id<-gsub("f","4",ped.cfm2$id)
ped.cfm2$id<-gsub("-","",ped.cfm2$id)
ped.cfm2$FATHER<-gsub("M","1",ped.cfm2$FATHER)
ped.cfm2$FATHER[is.na(ped.cfm2$FATHER)] <- 0
ped.cfm2$MOTHER<-gsub("F","2",ped.cfm2$MOTHER)
ped.cfm2$MOTHER[is.na(ped.cfm2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mWS")
write.table(ped.cfm2, file="cfm2mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CFMmWS <- readRDS("rep2CFM.mWS")
library(MCMCglmm)
posterior.mode(rep2CFMmWS$VCV)
#animal    round.sec        units 
#0.0018622586 0.0003164373 0.0010363658

# Testing significance of random effects
#logL for model with animal: 1233.827
#logL for model without animal: 1215.332
q <- 2*(1233.827-1215.332)
q
#q=36.99
pchisq(q, df=1, lower.tail = FALSE)
#p=1.187366e-09