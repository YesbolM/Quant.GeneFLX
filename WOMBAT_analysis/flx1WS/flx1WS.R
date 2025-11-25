#rep1flx WS

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
flx1<-rep1[rep1$treatment=="FLX",] 
#  WOMBAT does not allow any missing values, let's remove them
flx1<-na.omit(flx1)
#Select wing size for analysis
flx1WS<-flx1[c(7,1,5,13:14)]
flx1WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx1WS$TRAITNO<-ifelse(flx1WS$sex=="male",1,2)
flx1WS$sex<-ifelse(flx1WS$sex=="male",1,2)
#Reorder data
flx1WS<-flx1WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx1WS$round.sec <- as.numeric(flx1WS$round.sec)
#drop unused levels
flx1WS<-droplevels(flx1WS) 
str(flx1WS)
#round: 4 levels
length(unique(flx1WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1WS$animal<-gsub("M","1",flx1WS$animal)
flx1WS$animal<-gsub("F","2",flx1WS$animal)
flx1WS$animal<-gsub("m","3",flx1WS$animal)
flx1WS$animal<-gsub("f","4",flx1WS$animal)
flx1WS$animal<-gsub("-","",flx1WS$animal)
#Sort by animal ID
flx1WS <- flx1WS[order(flx1WS$animal),]

#Change name of round.sec to NSEC
names(flx1WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1WS")
write.table(flx1WS, file="flx1WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.flx1 <- read.csv("rep1flx.p.csv")
head(ped.flx1)
str(ped.flx1)
#Change ids into numeric values
ped.flx1$id<-gsub("M","1",ped.flx1$id)
ped.flx1$id<-gsub("F","2",ped.flx1$id)
ped.flx1$id<-gsub("m","3",ped.flx1$id)
ped.flx1$id<-gsub("f","4",ped.flx1$id)
ped.flx1$id<-gsub("-","",ped.flx1$id)
ped.flx1$FATHER<-gsub("M","1",ped.flx1$FATHER)
ped.flx1$FATHER[is.na(ped.flx1$FATHER)] <- 0
ped.flx1$MOTHER<-gsub("F","2",ped.flx1$MOTHER)
ped.flx1$MOTHER[is.na(ped.flx1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1WS")
write.table(ped.flx1, file="flx1WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1flx.WS <- readRDS("rep1flx.WS")
posterior.mode(rep1flx.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal 
#1.614706e-03                1.762379e-03 
#traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#1.762379e-03                2.232470e-03 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec 
#2.460493e-04                7.056426e-05 
#traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#7.056426e-05                3.975539e-04 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units 
#4.620382e-04                4.205792e-05 
#traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#4.205792e-05                3.586765e-04 

#Approximate CIs
#Va mWS
0.18273E-02-(1.96*0.366403E-03)
0.18273E-02+(1.96*0.366403E-03)
#Va fWS
0.257129E-02-(1.96*0.478395E-03)
0.257129E-02+(1.96*0.478395E-03)
#Cov WS
0.215900E-02-(1.96*0.346179E-03)
0.215900E-02+(1.96*0.346179E-03)
#rmf WS
0.996-(1.96*0.060)
0.996+(1.96*0.060)

#Log-likelihood test for covariance
#logL for model with cov: 2447.045
#logL for model without cov: 2407.657
q <- 2*(2447.045-2407.657)
q
#q=78.776
pchisq(q, df=1, lower.tail = FALSE)
#p=6.95671e-19
