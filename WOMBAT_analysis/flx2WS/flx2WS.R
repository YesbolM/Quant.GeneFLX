#rep2flx WS

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
flx2<-rep2[rep2$treatment=="FLX",] 
#  WOMBAT does not allow any missing values, let's remove them
flx2<-na.omit(flx2)
#Select wing size for analysis
flx2WS<-flx2[c(7,1,5,13:14)]
flx2WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx2WS$TRAITNO<-ifelse(flx2WS$sex=="male",1,2)
flx2WS$sex<-ifelse(flx2WS$sex=="male",1,2)
#Reorder data
flx2WS<-flx2WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx2WS$round.sec <- as.numeric(flx2WS$round.sec)
#drop unused levels
flx2WS<-droplevels(flx2WS) 
str(flx2WS)
#round: 4 levels
length(unique(flx2WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2WS$animal<-gsub("M","1",flx2WS$animal)
flx2WS$animal<-gsub("F","2",flx2WS$animal)
flx2WS$animal<-gsub("m","3",flx2WS$animal)
flx2WS$animal<-gsub("f","4",flx2WS$animal)
flx2WS$animal<-gsub("-","",flx2WS$animal)
#Sort by animal ID
flx2WS <- flx2WS[order(flx2WS$animal),]

#Change name of round.sec to NSEC
names(flx2WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2WS")
write.table(flx2WS, file="flx2WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.flx2 <- read.csv("rep2flx.p.csv")
head(ped.flx2)
str(ped.flx2)
#Change ids into numeric values
ped.flx2$id<-gsub("M","1",ped.flx2$id)
ped.flx2$id<-gsub("F","2",ped.flx2$id)
ped.flx2$id<-gsub("m","3",ped.flx2$id)
ped.flx2$id<-gsub("f","4",ped.flx2$id)
ped.flx2$id<-gsub("-","",ped.flx2$id)
ped.flx2$FATHER<-gsub("M","1",ped.flx2$FATHER)
ped.flx2$FATHER[is.na(ped.flx2$FATHER)] <- 0
ped.flx2$MOTHER<-gsub("F","2",ped.flx2$MOTHER)
ped.flx2$MOTHER[is.na(ped.flx2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2WS")
write.table(ped.flx2, file="flx2WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2flx.WS <- readRDS("rep2flx.WS")
posterior.mode(rep2flx.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal    traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#1.740457e-03                1.268251e-03                1.268251e-03                2.358883e-03 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#2.778483e-04                6.473318e-05                6.473318e-05                3.605006e-04 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units     traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#7.156357e-04                1.217686e-05                1.217686e-05                9.603044e-04 

#Approximate CIs
#Va mWS
0.220236E-02-(1.96*0.432085E-03)
0.220236E-02+(1.96*0.432085E-03)
#Va fWS
0.258153E-02-(1.96*0.539726E-03)
0.258153E-02+(1.96*0.539726E-03)
#Cov WS
0.180290E-02-(1.96*0.370948E-03)
0.180290E-02+(1.96*0.370948E-03)
#rmf WS
0.756-(1.96*0.096)
0.756+(1.96*0.096)

#Log-likelihood test for covariance
#logL for model with cov: 2381.044
#logL for model without cov: 2360.541
q <- 2*(2381.044-2360.541)
q
#q=41.006
pchisq(q, df=1, lower.tail = FALSE)
#p=1.517626e-10
