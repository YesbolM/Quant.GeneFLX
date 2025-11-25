#rep2cwt mWS

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
cwt2mWS<-cwt2m[,c(1,7,13:14)]
str(cwt2mWS)
#Make round.sec into integers by taking factor codes
cwt2mWS$round.sec <- as.numeric(cwt2mWS$round.sec)
#Drop unused levels
cwt2mWS<-droplevels(cwt2mWS) 
str(cwt2mWS)
#round: 4 levels
length(unique(cwt2mWS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2mWS$animal<-gsub("M","1",cwt2mWS$animal)
cwt2mWS$animal<-gsub("F","2",cwt2mWS$animal)
cwt2mWS$animal<-gsub("m","3",cwt2mWS$animal)
cwt2mWS$animal<-gsub("f","4",cwt2mWS$animal)
cwt2mWS$animal<-gsub("-","",cwt2mWS$animal)

#Change name of round.sec to NSEC
names(cwt2mWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mWS")
write.table(cwt2mWS, file="cwt2mWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2mWS")
write.table(ped.cwt2, file="cwt2mWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CWTmWS <- readRDS("rep2CWT.mWS")
library(MCMCglmm)
posterior.mode(rep2CWTmWS$VCV)
#animal    round.sec        units 
#0.0032082498 0.0005990073 0.0005721071 

# Testing significance of random effects
#logL for model with animal: 1174.744
#logL for model without animal: 1144.840
q <- 2*(1174.744-1144.840)
q
#q=59.808
pchisq(q, df=1, lower.tail = FALSE)
#p=1.045774e-14
