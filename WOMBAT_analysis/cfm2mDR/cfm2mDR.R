#rep2cfm mDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
cfm2mDR<-cfm2m[,c(1,8,13:14)]
str(cfm2mDR)
#Make round.sec into integers by taking factor codes
cfm2mDR$round.sec <- as.numeric(cfm2mDR$round.sec)
#Drop unused levels
cfm2mDR<-droplevels(cfm2mDR) 
str(cfm2mDR)
#round: 4 levels
length(unique(cfm2mDR$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm2mDR$animal<-gsub("M","1",cfm2mDR$animal)
cfm2mDR$animal<-gsub("F","2",cfm2mDR$animal)
cfm2mDR$animal<-gsub("m","3",cfm2mDR$animal)
cfm2mDR$animal<-gsub("f","4",cfm2mDR$animal)
cfm2mDR$animal<-gsub("-","",cfm2mDR$animal)

#Change name of round.sec to NSEC
names(cfm2mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mDR")
write.table(cfm2mDR, file="cfm2mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2mDR")
write.table(ped.cfm2, file="cfm2mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CFMmDR <- readRDS("rep2CFM.mDR")
library(MCMCglmm)
posterior.mode(rep2CFMmDR$VCV)
#animal    round.sec        units 
#1.235598073 0.003891989 2.173844383 

# Testing significance of random effects
#logL for model with animal: -545.011
#logL for model without animal: -550.937
q <- 2*(-545.011+550.937)
q
#q=11.852
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0005759945
