#rep2cwt fLA

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
cwt2f<-cwt2[cwt2$sex=="female",]
str(cwt2f)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
cwt2fLA<-cwt2f[,c(1,9,13:14)]
str(cwt2fLA)
#Make round.sec into integers by taking factor codes
cwt2fLA$round.sec <- as.numeric(cwt2fLA$round.sec)
#drop unused levels
cwt2fLA<-droplevels(cwt2fLA) 
str(cwt2fLA)
#round: 4 levels
length(unique(cwt2fLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2fLA$animal<-gsub("M","1",cwt2fLA$animal)
cwt2fLA$animal<-gsub("F","2",cwt2fLA$animal)
cwt2fLA$animal<-gsub("m","3",cwt2fLA$animal)
cwt2fLA$animal<-gsub("f","4",cwt2fLA$animal)
cwt2fLA$animal<-gsub("-","",cwt2fLA$animal)

#Change name of round.sec to NSEC
names(cwt2fLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2fLA")
write.table(cwt2fLA, file="cwt2fLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2fLA")
write.table(ped.cwt2, file="cwt2fLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CWTfLA <- readRDS("rep2CWT.fLA")
library(MCMCglmm)
posterior.mode(rep2CWTfLA$VCV)
#animal    round.sec        units 
#0.009843630 0.294866254 0.005356234 

# Testing significance of random effects
#logL for model with animal: -518.340
#logL for model without animal: -520.132
q <- 2*(-518.340+520.132)
q
#q=3.584
pchisq(q, df=1, lower.tail = FALSE)
#p=0.05833852

