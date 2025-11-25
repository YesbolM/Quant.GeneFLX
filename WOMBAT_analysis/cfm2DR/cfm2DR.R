#rep2cfm DR

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
cfm2<-rep2[rep2$treatment=="CFM",] 
#  WOMBAT does not allow any missing values, let's remove them
cfm2<-na.omit(cfm2)
#Select wing size for analysis
cfm2DR<-cfm2[c(8,1,5,13:14)]
cfm2DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm2DR$TRAITNO<-ifelse(cfm2DR$sex=="male",1,2)
cfm2DR$sex<-ifelse(cfm2DR$sex=="male",1,2)
#Reorder data
cfm2DR<-cfm2DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm2DR$round.sec <- as.numeric(cfm2DR$round.sec)
#drop unused levels
cfm2DR<-droplevels(cfm2DR) 
str(cfm2DR)
#round: 4 levels
length(unique(cfm2DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm2DR$animal<-gsub("M","1",cfm2DR$animal)
cfm2DR$animal<-gsub("F","2",cfm2DR$animal)
cfm2DR$animal<-gsub("m","3",cfm2DR$animal)
cfm2DR$animal<-gsub("f","4",cfm2DR$animal)
cfm2DR$animal<-gsub("-","",cfm2DR$animal)
#Sort by animal ID
cfm2DR <- cfm2DR[order(cfm2DR$animal),]

#Change name of round.sec to NSEC
names(cfm2DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2DR")
write.table(cfm2DR, file="cfm2DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cfm2 <- read.csv("rep2cfm.p.csv")
head(ped.cfm2)
str(ped.cfm2)
#Change ids into numeric values
ped.cfm2$id<-gsub("M","1",ped.cfm2$id)
ped.cfm2$id<-gsub("F","2",ped.cfm2$id)
ped.cfm2$id<-gsub("m","3",ped.cfm2$id)
ped.cfm2$id<-gsub("f","4",ped.cfm2$id)
ped.cfm2$id<-gsub("-","",ped.cfm2$id)
ped.cfm2$FATHER<-gsub("M","1",ped.cfm2$FATHER)
ped.cfm2$FATHER[is.na(ped.cfm2$FATHER)] <- 0
ped.cfm2$MOTHER<-gsub("F","2",ped.cfm2$MOTHER)
ped.cfm2$MOTHER[is.na(ped.cfm2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2DR")
write.table(ped.cfm2, file="cfm2DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cfm.DR <- readRDS("rep2cfm.DR")
posterior.mode(rep2cfm.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal    traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#0.939894074                 2.520872262                 2.520872262                 6.603425386 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#0.006792709                 0.004379413                 0.004379413                 0.022124255 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units     traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#2.227236934                -6.248283392                -6.248283392                17.901002933 

#Approximate CIs
#Va mDR
1.15507-(1.96*0.402033)
1.15507+(1.96*0.402033)
#Va fDR
6.53003-(1.96*2.91361)
6.53003+(1.96*2.91361)
#Cov DR
2.48181-(1.96*0.780031)
2.48181+(1.96*0.780031)
#rmf DR
0.904-(1.96*0.229)
0.904+(1.96*0.229)

#Log-likelihood test for covariance
#logL for model with cov: -1467.348
#logL for model without cov: -1474.915
q <- 2*(-1467.348+1474.915)
q
#q=15.134
pchisq(q, df=1, lower.tail = FALSE)
#p=0.0001001434
