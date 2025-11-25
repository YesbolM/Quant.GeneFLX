#rep1cwt fLA

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
cwt1<-rep1[rep1$treatment=="CWT",] 
#  WOMBAT does not allow any missing values, let's remove them
cwt1<-na.omit(cwt1) 
cwt1f<-cwt1[cwt1$sex=="female",]
str(cwt1f)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cwt1fLA<-cwt1f[,c(1,9,13:14)]
str(cwt1fLA)
#Make round.sec into integers by taking factor codes
cwt1fLA$round.sec <- as.numeric(cwt1fLA$round.sec)
#drop unused levels
cwt1fLA<-droplevels(cwt1fLA) 
str(cwt1fLA)
#round: 4 levels
length(unique(cwt1fLA$round.sec))
#round.sec: 28 levels

#Change ids into numeric values
cwt1fLA$animal<-gsub("M","1",cwt1fLA$animal)
cwt1fLA$animal<-gsub("F","2",cwt1fLA$animal)
cwt1fLA$animal<-gsub("m","3",cwt1fLA$animal)
cwt1fLA$animal<-gsub("f","4",cwt1fLA$animal)
cwt1fLA$animal<-gsub("-","",cwt1fLA$animal)

#Change name of round.sec to NSEC
names(cwt1fLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1fLA")
write.table(cwt1fLA, file="cwt1fLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt1 <- read.csv("rep1CWT.p.csv")
head(ped.cwt1)
str(ped.cwt1)
#Change ids into numeric values
ped.cwt1$id<-gsub("M","1",ped.cwt1$id)
ped.cwt1$id<-gsub("F","2",ped.cwt1$id)
ped.cwt1$id<-gsub("m","3",ped.cwt1$id)
ped.cwt1$id<-gsub("f","4",ped.cwt1$id)
ped.cwt1$id<-gsub("-","",ped.cwt1$id)
ped.cwt1$FATHER<-gsub("M","1",ped.cwt1$FATHER)
ped.cwt1$FATHER[is.na(ped.cwt1$FATHER)] <- 0
ped.cwt1$MOTHER<-gsub("F","2",ped.cwt1$MOTHER)
ped.cwt1$MOTHER[is.na(ped.cwt1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1fLA")
write.table(ped.cwt1, file="cwt1fLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CWTfLA <- readRDS("rep1CWT.fLA")
library(MCMCglmm)
posterior.mode(rep1CWTfLA$VCV)
#animal    round.sec        units 
#0.005812646 0.002393337 1.310340230 

# Testing significance of random effects
#logL for model with animal: -578.271
#logL for model without animal: -579.569
q <- 2*(-578.271+579.569)
q
#q=2.596
pchisq(q, df=1, lower.tail = FALSE)
#p=0.1071338