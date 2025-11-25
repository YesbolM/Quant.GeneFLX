#rep1cfm mLA

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
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cfm1mLA<-cfm1m[,c(1,9,13:14)]
str(cfm1mLA)
#Make round.sec into integers by taking factor codes
cfm1mLA$round.sec <- as.numeric(cfm1mLA$round.sec)
#drop unused levels
cfm1mLA<-droplevels(cfm1mLA) 
str(cfm1mLA)
#round: 4 levels
length(unique(cfm1mLA$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm1mLA$animal<-gsub("M","1",cfm1mLA$animal)
cfm1mLA$animal<-gsub("F","2",cfm1mLA$animal)
cfm1mLA$animal<-gsub("m","3",cfm1mLA$animal)
cfm1mLA$animal<-gsub("f","4",cfm1mLA$animal)
cfm1mLA$animal<-gsub("-","",cfm1mLA$animal)

#Change name of round.sec to NSEC
names(cfm1mLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mLA")
write.table(cfm1mLA, file="cfm1mLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mLA")
write.table(ped.cfm1, file="cfm1mLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CFMmLA <- readRDS("rep1CFM.mLA")
library(MCMCglmm)
posterior.mode(rep1CFMmLA$VCV)
#animal    round.sec        units 
#0.686755651 0.002932446 0.612958039 

# Testing significance of random effects
#logL for model with animal: -801.581
#logL for model without animal: -813.733
q <- 2*(-801.581+813.733)
q
#q=24.304
pchisq(q, df=1, lower.tail = FALSE)
#p=8.226805e-07