#rep1cwt LA

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
#Select wing size for analysis
cwt1LA<-cwt1[c(9,1,5,13:14)]
cwt1LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt1LA$TRAITNO<-ifelse(cwt1LA$sex=="male",1,2)
cwt1LA$sex<-ifelse(cwt1LA$sex=="male",1,2)
#Reorder data
cwt1LA<-cwt1LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt1LA$round.sec <- as.numeric(cwt1LA$round.sec)
#Drop unused levels
cwt1LA<-droplevels(cwt1LA) 
str(cwt1LA)
#round: 4 levels
length(unique(cwt1LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt1LA$animal<-gsub("M","1",cwt1LA$animal)
cwt1LA$animal<-gsub("F","2",cwt1LA$animal)
cwt1LA$animal<-gsub("m","3",cwt1LA$animal)
cwt1LA$animal<-gsub("f","4",cwt1LA$animal)
cwt1LA$animal<-gsub("-","",cwt1LA$animal)
#Sort by animal ID
cwt1LA <- cwt1LA[order(cwt1LA$animal),]

#Change name of round.sec to NSEC
names(cwt1LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1LA")
write.table(cwt1LA, file="cwt1LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt1 <- read.csv("rep1cwt.p.csv")
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1LA")
write.table(ped.cwt1, file="cwt1LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cwt.LA <- readRDS("rep1cwt.LA")
posterior.mode(rep1cwt.LA$VCV)
#traitmLA.active:traitmLA.active.animal    traitfLA.active:traitmLA.active.animal 
#0.6389641318                              0.5388375348 
#traitmLA.active:traitfLA.active.animal    traitfLA.active:traitfLA.active.animal 
#0.5388375348                              0.5541888111 
#traitmLA.active:traitmLA.active.round.sec traitfLA.active:traitmLA.active.round.sec 
#0.0033353788                              0.0005784652 
#traitmLA.active:traitfLA.active.round.sec traitfLA.active:traitfLA.active.round.sec 
#0.0005784652                              0.0055531489 
#traitmLA.active:traitmLA.active.units     traitfLA.active:traitmLA.active.units 
#0.5474514706                             -0.8799710711 
#traitmLA.active:traitfLA.active.units     traitfLA.active:traitfLA.active.units 
#-0.8799710711                              1.2916614589 

#Approximate CIs
#Va mLA
2.83206-(1.96*1.28245)
2.83206+(1.96*1.28245)
#Va fLA
1.46874-(1.96*0.959819)
1.46874+(1.96*0.959819)
#Cov LA
1.42614-(1.96*0.744966)
1.42614+(1.96*0.744966)
#rmf LA
0.699-(1.96*0.358)
0.699+(1.96*0.358)

#Log-likelihood test for covariance
#logL for model with cov: -1453.296
#logL for model without cov: -1456.354
q <- 2*(-1453.296+1456.354)
q
#q=6.116
pchisq(q, df=1, lower.tail = FALSE)
#p=0.01339636

