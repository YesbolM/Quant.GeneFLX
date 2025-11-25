#rep2cfm mLA

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
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cfm2mLA<-cfm2m[,c(1,9,13:14)]
str(cfm2mLA)
#Make round.sec into integers by taking factor codes
cfm2mLA$round.sec <- as.numeric(cfm2mLA$round.sec)
#drop unused levels
cfm2mLA<-droplevels(cfm2mLA) 
str(cfm2mLA)
#round: 4 levels
length(unique(cfm2mLA$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm2mLA$animal<-gsub("M","1",cfm2mLA$animal)
cfm2mLA$animal<-gsub("F","2",cfm2mLA$animal)
cfm2mLA$animal<-gsub("m","3",cfm2mLA$animal)
cfm2mLA$animal<-gsub("f","4",cfm2mLA$animal)
cfm2mLA$animal<-gsub("-","",cfm2mLA$animal)

#Change name of round.sec to NSEC
names(cfm2mLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mLA")
write.table(cfm2mLA, file="cfm2mLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mLA")
write.table(ped.cfm2, file="cfm2mLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CFMmLA <- readRDS("rep2CFM.mLA")
library(MCMCglmm)
posterior.mode(rep2CFMmLA$VCV)
#animal    round.sec        units 
#0.2350769 0.1254656 0.5940939

# Testing significance of random effects
#logL for model with animal: -841.819
#logL for model without animal: -845.985
q <- 2*(-841.819+845.985)
q
#q=8.332
pchisq(q, df=1, lower.tail = FALSE)
#p=0.003895275