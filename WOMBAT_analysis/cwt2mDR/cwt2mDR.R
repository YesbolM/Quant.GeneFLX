#rep2cwt mDR

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
cwt2<-rep2[rep2$treatment=="CWT",] 
#  WOMBAT does not allow any missing values, let's remove them
cwt2<-na.omit(cwt2)
cwt2m<-cwt2[cwt2$sex=="male",]
str(cwt2m)
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cwt2mDR<-cwt2m[,c(1,8,13:14)]
str(cwt2mDR)
#Make round.sec into integers by taking factor codes
cwt2mDR$round.sec <- as.numeric(cwt2mDR$round.sec)
#Drop unused levels
cwt2mDR<-droplevels(cwt2mDR) 
str(cwt2mDR)
#round: 4 levels
length(unique(cwt2mDR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2mDR$animal<-gsub("M","1",cwt2mDR$animal)
cwt2mDR$animal<-gsub("F","2",cwt2mDR$animal)
cwt2mDR$animal<-gsub("m","3",cwt2mDR$animal)
cwt2mDR$animal<-gsub("f","4",cwt2mDR$animal)
cwt2mDR$animal<-gsub("-","",cwt2mDR$animal)

#Change name of round.sec to NSEC
names(cwt2mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mDR")
write.table(cwt2mDR, file="cwt2mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt2 <- read.csv("rep2CWT.p.csv")
head(ped.cwt2)
str(ped.cwt2)
#Change ids into numeric values
ped.cwt2$id<-gsub("M","1",ped.cwt2$id)
ped.cwt2$id<-gsub("F","2",ped.cwt2$id)
ped.cwt2$id<-gsub("m","3",ped.cwt2$id)
ped.cwt2$id<-gsub("f","4",ped.cwt2$id)
ped.cwt2$id<-gsub("-","",ped.cwt2$id)
ped.cwt2$FATHER<-gsub("M","1",ped.cwt2$FATHER)
ped.cwt2$FATHER[is.na(ped.cwt2$FATHER)] <- 0
ped.cwt2$MOTHER<-gsub("F","2",ped.cwt2$MOTHER)
ped.cwt2$MOTHER[is.na(ped.cwt2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mDR")
write.table(ped.cwt2, file="cwt2mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CWTmDR <- readRDS("rep2CWT.mDR")
library(MCMCglmm)
posterior.mode(rep2CWTmDR$VCV)
#animal    round.sec        units 
#1.88191017 0.01114971 3.25290647 

# Testing significance of random effects
#logL for model with animal: -668.919
#logL for model without animal: -673.826
q <- 2*(-668.919+673.826)
q
#q=9.814
pchisq(q, df=1, lower.tail = FALSE)
#p=0.001731884