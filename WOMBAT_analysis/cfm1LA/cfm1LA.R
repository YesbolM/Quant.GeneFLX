#rep1cfm LA

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
#Select wing size for analysis
cfm1LA<-cfm1[c(9,1,5,13:14)]
cfm1LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm1LA$TRAITNO<-ifelse(cfm1LA$sex=="male",1,2)
cfm1LA$sex<-ifelse(cfm1LA$sex=="male",1,2)
#Reorder data
cfm1LA<-cfm1LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm1LA$round.sec <- as.numeric(cfm1LA$round.sec)
#Drop unused levels
cfm1LA<-droplevels(cfm1LA) 
str(cfm1LA)
#round: 4 levels
length(unique(cfm1LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm1LA$animal<-gsub("M","1",cfm1LA$animal)
cfm1LA$animal<-gsub("F","2",cfm1LA$animal)
cfm1LA$animal<-gsub("m","3",cfm1LA$animal)
cfm1LA$animal<-gsub("f","4",cfm1LA$animal)
cfm1LA$animal<-gsub("-","",cfm1LA$animal)
#Sort by animal ID
cfm1LA <- cfm1LA[order(cfm1LA$animal),]

#Change name of round.sec to NSEC
names(cfm1LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1LA")
write.table(cfm1LA, file="cfm1LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1LA")
write.table(ped.cfm1, file="cfm1LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cfm.LA <- readRDS("rep1cfm.LA")
posterior.mode(rep1cfm.LA$VCV)
#traitmLA.active:traitmLA.active.animal 
#0.805625352 
#traitfLA.active:traitmLA.active.animal 
#0.359061573 
#traitmLA.active:traitfLA.active.animal 
#0.359061573 
#traitfLA.active:traitfLA.active.animal 
#0.136327939 
#traitmLA.active:traitmLA.active.round.sec 
#0.003002279 
#traitfLA.active:traitmLA.active.round.sec 
#-0.000374850 
#traitmLA.active:traitfLA.active.round.sec 
#-0.000374850 
#traitfLA.active:traitfLA.active.round.sec 
#0.003348680 
#traitmLA.active:traitmLA.active.units 
#0.549776963 
#traitfLA.active:traitmLA.active.units 
#0.855796635 
#traitmLA.active:traitfLA.active.units 
#0.855796635 
#traitfLA.active:traitfLA.active.units 
#1.429147100 

#Approximate CIs
#Va mLA
4.70135-(1.96*1.36243)
4.70135+(1.96*1.36243)
#Va fLA
1.13812-(1.96*0.550150)
1.13812+(1.96*0.550150)
#Cov LA
1.16606-(1.96*0.604650)
1.16606+(1.96*0.604650)
#rmf LA
0.504-(1.96*0.229)
0.504+(1.96*0.229)

#Log-likelihood test for covariance
#logL for model with cov: -1387.478
#logL for model without cov: -1392.432
q <- 2*(-1387.478+1392.432)
q
#q=9.908
pchisq(q, df=1, lower.tail = FALSE)
#p=0.001645619
