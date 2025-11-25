#rep1cwt WS

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
cwt1WS<-cwt1[c(7,1,5,13:14)]
cwt1WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt1WS$TRAITNO<-ifelse(cwt1WS$sex=="male",1,2)
cwt1WS$sex<-ifelse(cwt1WS$sex=="male",1,2)
#Reorder data
cwt1WS<-cwt1WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt1WS$round.sec <- as.numeric(cwt1WS$round.sec)
#drop unused levels
cwt1WS<-droplevels(cwt1WS) 
str(cwt1WS)
#round: 4 levels
length(unique(cwt1WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt1WS$animal<-gsub("M","1",cwt1WS$animal)
cwt1WS$animal<-gsub("F","2",cwt1WS$animal)
cwt1WS$animal<-gsub("m","3",cwt1WS$animal)
cwt1WS$animal<-gsub("f","4",cwt1WS$animal)
cwt1WS$animal<-gsub("-","",cwt1WS$animal)
#Sort by animal ID
cwt1WS <- cwt1WS[order(cwt1WS$animal),]

#Change name of round.sec to NSEC
names(cwt1WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1WS")
write.table(cwt1WS, file="cwt1WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1WS")
write.table(ped.cwt1, file="cwt1WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1cwt.WS <- readRDS("rep1cwt.WS")
posterior.mode(rep1cwt.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal 
#1.850073e-03                1.617926e-03 
#traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#1.617926e-03                2.152812e-03 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec 
#4.114821e-04                2.261106e-04 
#traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#2.261106e-04                5.920097e-04 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units 
#6.520133e-04               -2.813732e-05 
#traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#-2.813732e-05                6.454667e-04 

#Approximate CIs
#Va mWS
0.214594E-02-(1.96*0.422407E-03)
0.214594E-02+(1.96*0.422407E-03)
#Va fWS
0.229516E-02-(1.96*0.558021E-03)
0.229516E-02+(1.96*0.558021E-03)
#Cov WS
0.179470E-02-(1.96*0.376741E-03)
0.179470E-02+(1.96*0.376741E-03)
#rmf WS
0.809-(1.96*0.100)
0.809+(1.96*0.100)

#Log-likelihood test for covariance
#logL for model with cov: 2249.696
#logL for model without cov: 2232.627
q <- 2*(2249.696-2232.627)
q
#q=34.138
pchisq(q, df=1, lower.tail = FALSE)
#p=5.133886e-09
