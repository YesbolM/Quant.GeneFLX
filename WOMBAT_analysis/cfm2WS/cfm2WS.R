#rep2cfm WS

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
cfm2WS<-cfm2[c(7,1,5,13:14)]
cfm2WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm2WS$TRAITNO<-ifelse(cfm2WS$sex=="male",1,2)
cfm2WS$sex<-ifelse(cfm2WS$sex=="male",1,2)
#Reorder data
cfm2WS<-cfm2WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm2WS$round.sec <- as.numeric(cfm2WS$round.sec)
#drop unused levels
cfm2WS<-droplevels(cfm2WS) 
str(cfm2WS)
#round: 4 levels
length(unique(cfm2WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm2WS$animal<-gsub("M","1",cfm2WS$animal)
cfm2WS$animal<-gsub("F","2",cfm2WS$animal)
cfm2WS$animal<-gsub("m","3",cfm2WS$animal)
cfm2WS$animal<-gsub("f","4",cfm2WS$animal)
cfm2WS$animal<-gsub("-","",cfm2WS$animal)
#Sort by animal ID
cfm2WS <- cfm2WS[order(cfm2WS$animal),]

#Change name of round.sec to NSEC
names(cfm2WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2WS")
write.table(cfm2WS, file="cfm2WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2WS")
write.table(ped.cfm2, file="cfm2WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cfm.WS <- readRDS("rep2cfm.WS")
posterior.mode(rep2cfm.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal    traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#0.0017196119                0.0011640266                0.0011640266                0.0017689080 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#0.0002786700                0.0001004975                0.0001004975                0.0003294403 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units     traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#0.0012785949               -0.0001131584               -0.0001131584                0.0022053953 

#Approximate CIs
#Va mWS
0.174372E-02-(1.96*0.409462E-03)
0.174372E-02+(1.96*0.409462E-03)
#Va fWS
0.189622E-02-(1.96*0.522023E-03)
0.189622E-02+(1.96*0.522023E-03)
#Cov WS
0.154614E-02-(1.96*0.346273E-03)
0.154614E-02+(1.96*0.346273E-03)
#rmf WS
0.850-(1.96*0.126)
0.850+(1.96*0.126)

#Log-likelihood test for covariance
#logL for model with cov: 2278.778
#logL for model without cov: 2258.552
q <- 2*(2278.778-2258.552)
q
#q=40.452
pchisq(q, df=1, lower.tail = FALSE)
#p=2.015061e-10