#rep1cwt fWS

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
#Select wing size for analysis
#We need to include animal, WS, round, and round.sec
cwt1fWS<-cwt1f[,c(1,7,13:14)]
str(cwt1fWS)
#Make round.sec into integers by taking factor codes
cwt1fWS$round.sec <- as.numeric(cwt1fWS$round.sec)
#Drop unused levels
cwt1fWS<-droplevels(cwt1fWS) 
str(cwt1fWS)
#round: 4 levels
length(unique(cwt1fWS$round.sec))
#round.sec: 28 levels

#Change ids into numeric values
cwt1fWS$animal<-gsub("M","1",cwt1fWS$animal)
cwt1fWS$animal<-gsub("F","2",cwt1fWS$animal)
cwt1fWS$animal<-gsub("m","3",cwt1fWS$animal)
cwt1fWS$animal<-gsub("f","4",cwt1fWS$animal)
cwt1fWS$animal<-gsub("-","",cwt1fWS$animal)

#Change name of round.sec to NSEC
names(cwt1fWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1fWS")
write.table(cwt1fWS, file="cwt1fWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1fWS")
write.table(ped.cwt1, file="cwt1fWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1CWTfWS <- readRDS("rep1CWT.fWS")
library(MCMCglmm)
posterior.mode(rep1CWTfWS$VCV)
#animal    round.sec        units 
#0.0020532729 0.0005070573 0.0009641562 

# Testing significance of random effects
#logL for model with animal: 943.173
#logL for model without animal: 929.052
q <- 2*(943.173-929.052)
q
#q=28.242
pchisq(q, df=1, lower.tail = FALSE)
#p=1.070562e-07
