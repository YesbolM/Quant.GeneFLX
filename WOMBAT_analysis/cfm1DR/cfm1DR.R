#rep1cfm DR

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
cfm1DR<-cfm1[c(8,1,5,13:14)]
cfm1DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm1DR$TRAITNO<-ifelse(cfm1DR$sex=="male",1,2)
cfm1DR$sex<-ifelse(cfm1DR$sex=="male",1,2)
#Reorder data
cfm1DR<-cfm1DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm1DR$round.sec <- as.numeric(cfm1DR$round.sec)
#drop unused levels
cfm1DR<-droplevels(cfm1DR) 
str(cfm1DR)
#round: 4 levels
length(unique(cfm1DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm1DR$animal<-gsub("M","1",cfm1DR$animal)
cfm1DR$animal<-gsub("F","2",cfm1DR$animal)
cfm1DR$animal<-gsub("m","3",cfm1DR$animal)
cfm1DR$animal<-gsub("f","4",cfm1DR$animal)
cfm1DR$animal<-gsub("-","",cfm1DR$animal)
#Sort by animal ID
cfm1DR <- cfm1DR[order(cfm1DR$animal),]

#Change name of round.sec to NSEC
names(cfm1DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1DR")
write.table(cfm1DR, file="cfm1DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1DR")
write.table(ped.cfm1, file="cfm1DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cfm.DR <- readRDS("rep1cfm.DR")
posterior.mode(rep1cfm.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal 
#2.2507932                   3.0196881 
#traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#3.0196881                   7.8518145 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec 
#0.6570223                   1.7685777 
#traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#1.7685777                   4.9787394 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units 
#1.9546576                  -6.1628371 
#traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#-6.1628371                  17.7248881 

#Approximate CIs
#Va mDR
2.35096-(1.96*0.613035)
2.35096+(1.96*0.613035)
#Va fDR
10.9300-(1.96*3.69482)
10.9300+(1.96*3.69482)
#Cov DR
2.49691-(1.96*1.05114)
2.49691+(1.96*1.05114)
#rmf DR
0.493-(1.96*0.181)
0.493+(1.96*0.181)

#Log-likelihood test for covariance
#logL for model with cov: -1545.757
#logL for model without cov: -1550.231
q <- 2*(-1545.757+1550.231)
q
#q=8.948
pchisq(q, df=1, lower.tail = FALSE)
#p=0.002777735
