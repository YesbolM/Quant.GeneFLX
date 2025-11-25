#rep1cfm fDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
cfm1fDR<-cfm1f[,c(1,8,13:14)]
str(cfm1fDR)
#Make round.sec into integers by taking factor codes
cfm1fDR$round.sec <- as.numeric(cfm1fDR$round.sec)
#Drop unused levels
cfm1fDR<-droplevels(cfm1fDR) 
str(cfm1fDR)
#round: 4 levels
length(unique(cfm1fDR$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm1fDR$animal<-gsub("M","1",cfm1fDR$animal)
cfm1fDR$animal<-gsub("F","2",cfm1fDR$animal)
cfm1fDR$animal<-gsub("m","3",cfm1fDR$animal)
cfm1fDR$animal<-gsub("f","4",cfm1fDR$animal)
cfm1fDR$animal<-gsub("-","",cfm1fDR$animal)

#Change name of round.sec to NSEC
names(cfm1fDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1fDR")
write.table(cfm1fDR, file="cfm1fDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1fDR")
write.table(ped.cfm1, file="cfm1fDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CFMfDR <- readRDS("rep1CFM.fDR")
library(MCMCglmm)
posterior.mode(rep1CFMfDR$VCV)
#animal    round.sec        units 
#10.804217  4.273915 16.194070 

# Testing significance of random effects
#logL for model with animal: -974.878
#logL for model without animal: -982.410
q <- 2*(-974.878+982.410)
q
#q=15.064
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0001039265
