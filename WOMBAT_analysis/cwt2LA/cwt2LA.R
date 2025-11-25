#rep2cwt LA

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
#Select wing size for analysis
cwt2LA<-cwt2[c(9,1,5,13:14)]
cwt2LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt2LA$TRAITNO<-ifelse(cwt2LA$sex=="male",1,2)
cwt2LA$sex<-ifelse(cwt2LA$sex=="male",1,2)
#Reorder data
cwt2LA<-cwt2LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt2LA$round.sec <- as.numeric(cwt2LA$round.sec)
#Drop unused levels
cwt2LA<-droplevels(cwt2LA) 
str(cwt2LA)
#round: 4 levels
length(unique(cwt2LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2LA$animal<-gsub("M","1",cwt2LA$animal)
cwt2LA$animal<-gsub("F","2",cwt2LA$animal)
cwt2LA$animal<-gsub("m","3",cwt2LA$animal)
cwt2LA$animal<-gsub("f","4",cwt2LA$animal)
cwt2LA$animal<-gsub("-","",cwt2LA$animal)
#Sort by animal ID
cwt2LA <- cwt2LA[order(cwt2LA$animal),]

#Change name of round.sec to NSEC
names(cwt2LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2LA")
write.table(cwt2LA, file="cwt2LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt2 <- read.csv("rep2cwt.p.csv")
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2LA")
write.table(ped.cwt2, file="cwt2LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cwt.LA <- readRDS("rep2cwt.LA")
posterior.mode(rep2cwt.LA$VCV)
#traitmLA.active:traitmLA.active.animal    traitfLA.active:traitmLA.active.animal 
#0.952931132                               0.575690348 
#traitmLA.active:traitfLA.active.animal    traitfLA.active:traitfLA.active.animal 
#0.575690348                               0.602646934 
#traitmLA.active:traitmLA.active.round.sec traitfLA.active:traitmLA.active.round.sec 
#0.005058656                               0.105894886 
#traitmLA.active:traitfLA.active.round.sec traitfLA.active:traitfLA.active.round.sec 
#0.105894886                               0.226650751 
#traitmLA.active:traitmLA.active.units     traitfLA.active:traitmLA.active.units 
#0.007367713                               0.298059765 
#traitmLA.active:traitfLA.active.units     traitfLA.active:traitfLA.active.units 
#0.298059765                               0.560551775 

#Approximate CIs
#Va mLA
1.72995-(1.96*0.658901)
1.72995+(1.96*0.658901)
#Va fLA
0.666965-(1.96*0.381841)
0.666965+(1.96*0.381841)
#Cov LA
1.07335-(1.96*0.360916)
1.07335+(1.96*0.360916)
#rmf LA
0.999-(1.96*0.304)
0.999+(1.96*0.304)

#Log-likelihood test for covariance
#logL for model with cov: -1170.495
#logL for model without cov: -1176.979
q <- 2*(-1170.495+1176.979)
q
#q=12.968
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0003168603
