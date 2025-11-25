#rep2cwt WS

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
cwt2WS<-cwt2[c(7,1,5,13:14)]
cwt2WS$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cwt2WS$TRAITNO<-ifelse(cwt2WS$sex=="male",1,2)
cwt2WS$sex<-ifelse(cwt2WS$sex=="male",1,2)
#Reorder data
cwt2WS<-cwt2WS[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cwt2WS$round.sec <- as.numeric(cwt2WS$round.sec)
#drop unused levels
cwt2WS<-droplevels(cwt2WS) 
str(cwt2WS)
#round: 4 levels
length(unique(cwt2WS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2WS$animal<-gsub("M","1",cwt2WS$animal)
cwt2WS$animal<-gsub("F","2",cwt2WS$animal)
cwt2WS$animal<-gsub("m","3",cwt2WS$animal)
cwt2WS$animal<-gsub("f","4",cwt2WS$animal)
cwt2WS$animal<-gsub("-","",cwt2WS$animal)
#Sort by animal ID
cwt2WS <- cwt2WS[order(cwt2WS$animal),]

#Change name of round.sec to NSEC
names(cwt2WS)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2WS")
write.table(cwt2WS, file="cwt2WS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2WS")
write.table(ped.cwt2, file="cwt2WS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cwt.WS <- readRDS("rep2cwt.WS")
posterior.mode(rep2cwt.WS$VCV)
#traitmWS:traitmWS.animal    traitfWS:traitmWS.animal    traitmWS:traitfWS.animal    traitfWS:traitfWS.animal 
#2.801250e-03                2.169691e-03                2.169691e-03                3.374578e-03 
#traitmWS:traitmWS.round.sec traitfWS:traitmWS.round.sec traitmWS:traitfWS.round.sec traitfWS:traitfWS.round.sec 
#6.967135e-04               -1.592182e-04               -1.592182e-04                6.210159e-04 
#traitmWS:traitmWS.units     traitfWS:traitmWS.units     traitmWS:traitfWS.units     traitfWS:traitfWS.units 
#5.529878e-04                5.349630e-06                5.349630e-06                9.166396e-04 

#Approximate CIs
#Va mWS
0.326004E-02-(1.96*0.613600E-03)
0.326004E-02+(1.96*0.613600E-03)
#Va fWS
0.485488E-02-(1.96*0.882881E-03)
0.485488E-02+(1.96*0.882881E-03)
#Cov WS
0.252721E-02-(1.96*0.543697E-03)
0.252721E-02+(1.96*0.543697E-03)
#rmf WS
0.635-(1.96*0.096)
0.635+(1.96*0.096)

#Log-likelihood test for covariance
#logL for model with cov: 2209.055
#logL for model without cov: 2195.129
q <- 2*(2209.055-2195.129)
q
#q=27.852
pchisq(q, df=1, lower.tail = FALSE)
#p=1.309586e-07
