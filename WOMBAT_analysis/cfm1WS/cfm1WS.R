#rep1cfm WS

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
cfm1WS<-cfm1[c(7,1,5,13:14)]
cfm1WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm1WS$TRAITNO<-ifelse(cfm1WS$sex=="male",1,2)
cfm1WS$sex<-ifelse(cfm1WS$sex=="male",1,2)
#Reorder data
cfm1WS<-cfm1WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm1WS$round.sec <- as.numeric(cfm1WS$round.sec)
#drop unused levels
cfm1WS<-droplevels(cfm1WS) 
str(cfm1WS)
#round: 4 levels
length(unique(cfm1WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm1WS$animal<-gsub("M","1",cfm1WS$animal)
cfm1WS$animal<-gsub("F","2",cfm1WS$animal)
cfm1WS$animal<-gsub("m","3",cfm1WS$animal)
cfm1WS$animal<-gsub("f","4",cfm1WS$animal)
cfm1WS$animal<-gsub("-","",cfm1WS$animal)
#Sort by animal ID
cfm1WS <- cfm1WS[order(cfm1WS$animal),]

#Change name of round.sec to NSEC
names(cfm1WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1WS")
write.table(cfm1WS, file="cfm1WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1WS")
write.table(ped.cfm1, file="cfm1WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cfm.WS <- readRDS("rep1cfm.WS")
posterior.mode(rep1cfm.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal 
#1.740797e-03                1.323846e-03 
#traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#1.323846e-03                2.149574e-03 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec 
#2.237250e-04                5.965911e-05 
#traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#5.965911e-05                3.475242e-04 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units 
#9.321869e-04                2.781808e-04 
#traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#2.781808e-04                7.794911e-04 

#Approximate CIs
#Va mWS
0.154026E-02-(1.96*0.371820E-03)
0.154026E-02+(1.96*0.371820E-03)
#Va fWS
0.200811E-02-(1.96*0.465042E-03)
0.200811E-02+(1.96*0.465042E-03)
#Cov WS
0.153248E-02-(1.96*0.315612E-03)
0.153248E-02+(1.96*0.315612E-03)
#rmf WS
0.871-(1.96*0.107)
0.871+(1.96*0.107)

#Log-likelihood test for covariance
#logL for model with cov: 2359.033
#logL for model without cov: 2340.809
q <- 2*(2359.033-2340.809)
q
#q=36.448
pchisq(q, df=1, lower.tail = FALSE)
#p=1.567939e-09
