#rep2cfm fDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
cfm2fDR<-cfm2f[,c(1,8,13:14)]
str(cfm2fDR)
#Make round.sec into integers by taking factor codes
cfm2fDR$round.sec <- as.numeric(cfm2fDR$round.sec)
#Drop unused levels
cfm2fDR<-droplevels(cfm2fDR) 
str(cfm2fDR)
#round: 4 levels
length(unique(cfm2fDR$round.sec))
#round.sec: 27 levels

#Change ids into numeric values
cfm2fDR$animal<-gsub("M","1",cfm2fDR$animal)
cfm2fDR$animal<-gsub("F","2",cfm2fDR$animal)
cfm2fDR$animal<-gsub("m","3",cfm2fDR$animal)
cfm2fDR$animal<-gsub("f","4",cfm2fDR$animal)
cfm2fDR$animal<-gsub("-","",cfm2fDR$animal)

#Change name of round.sec to NSEC
names(cfm2fDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2fDR")
write.table(cfm2fDR, file="cfm2fDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2fDR")
write.table(ped.cfm2, file="cfm2fDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CFMfDR <- readRDS("rep2CFM.fDR")
library(MCMCglmm)
posterior.mode(rep2CFMfDR$VCV)
#animal    round.sec        units 
#8.610229616  0.005502964 16.861336038 

# Testing significance of random effects
#logL for model with animal: -934.393
#logL for model without animal: -937.940
q <- 2*(-934.393+937.940)
q
#q=7.094
pchisq(q, df=1, lower.tail = FALSE)
#p=0.007734243
