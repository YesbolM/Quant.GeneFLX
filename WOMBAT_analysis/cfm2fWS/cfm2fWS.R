#rep2cfm fWS

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
cfm2f<-cfm2[cfm2$sex=="female",]
str(cfm2f)
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cfm2fWS<-cfm2f[,c(1,7,13:14)]
str(cfm2fWS)
#Make round.sec into integers by taking factor codes
cfm2fWS$round.sec <- as.numeric(cfm2fWS$round.sec)
#Drop unused levels
cfm2fWS<-droplevels(cfm2fWS) 
str(cfm2fWS)
#round: 4 levels
length(unique(cfm2fWS$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm2fWS$animal<-gsub("M","1",cfm2fWS$animal)
cfm2fWS$animal<-gsub("F","2",cfm2fWS$animal)
cfm2fWS$animal<-gsub("m","3",cfm2fWS$animal)
cfm2fWS$animal<-gsub("f","4",cfm2fWS$animal)
cfm2fWS$animal<-gsub("-","",cfm2fWS$animal)

#Change name of round.sec to NSEC
names(cfm2fWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2fWS")
write.table(cfm2fWS, file="cfm2fWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2fWS")
write.table(ped.cfm2, file="cfm2fWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CFMfWS <- readRDS("rep2CFM.fWS")
library(MCMCglmm)
posterior.mode(rep2CFMfWS$VCV)
#animal    round.sec        units 
#0.0017793122 0.0003894269 0.0017439748 

# Testing significance of random effects
#logL for model with animal: 1013.566
#logL for model without animal: 1005.309
q <- 2*(1013.566-1005.309)
q
#q=16.514
pchisq(q, df=1, lower.tail = FALSE)
#p=4.829218e-05