#rep2cwt mLA

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
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cwt2mLA<-cwt2m[,c(1,9,13:14)]
str(cwt2mLA)
#Make round.sec into integers by taking factor codes
cwt2mLA$round.sec <- as.numeric(cwt2mLA$round.sec)
#drop unused levels
cwt2mLA<-droplevels(cwt2mLA) 
str(cwt2mLA)
#round: 4 levels
length(unique(cwt2mLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2mLA$animal<-gsub("M","1",cwt2mLA$animal)
cwt2mLA$animal<-gsub("F","2",cwt2mLA$animal)
cwt2mLA$animal<-gsub("m","3",cwt2mLA$animal)
cwt2mLA$animal<-gsub("f","4",cwt2mLA$animal)
cwt2mLA$animal<-gsub("-","",cwt2mLA$animal)

#Change name of round.sec to NSEC
names(cwt2mLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mLA")
write.table(cwt2mLA, file="cwt2mLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mLA")
write.table(ped.cwt2, file="cwt2mLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CWTmLA <- readRDS("rep2CWT.mLA")
library(MCMCglmm)
posterior.mode(rep2CWTmLA$VCV)
#animal    round.sec        units 
#0.837254549 0.003102857 0.005346287 

# Testing significance of random effects
#logL for model with animal: -663.357
#logL for model without animal: -668.017
q <- 2*(-663.357+668.017)
q
#q=9.32
pchisq(q, df=1, lower.tail = FALSE)
#p=0.00226666