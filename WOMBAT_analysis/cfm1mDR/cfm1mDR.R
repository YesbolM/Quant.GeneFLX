#rep1cfm mDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
cfm1mDR<-cfm1m[,c(1,8,13:14)]
str(cfm1mDR)
#Make round.sec into integers by taking factor codes
cfm1mDR$round.sec <- as.numeric(cfm1mDR$round.sec)
#Drop unused levels
cfm1mDR<-droplevels(cfm1mDR) 
str(cfm1mDR)
#round: 4 levels
length(unique(cfm1mDR$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm1mDR$animal<-gsub("M","1",cfm1mDR$animal)
cfm1mDR$animal<-gsub("F","2",cfm1mDR$animal)
cfm1mDR$animal<-gsub("m","3",cfm1mDR$animal)
cfm1mDR$animal<-gsub("f","4",cfm1mDR$animal)
cfm1mDR$animal<-gsub("-","",cfm1mDR$animal)

#Change name of round.sec to NSEC
names(cfm1mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mDR")
write.table(cfm1mDR, file="cfm1mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1mDR")
write.table(ped.cfm1, file="cfm1mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CFMmDR <- readRDS("rep1CFM.mDR")
library(MCMCglmm)
posterior.mode(rep1CFMmDR$VCV)
#animal    round.sec        units 
#2.28341899 0.01101003 1.59089363 

# Testing significance of random effects
#logL for model with animal: -583.691
#logL for model without animal: -597.872
q <- 2*(-583.691+597.872)
q
#q=28.362
pchisq(q, df=1, lower.tail = FALSE)
#p=1.006211e-07
