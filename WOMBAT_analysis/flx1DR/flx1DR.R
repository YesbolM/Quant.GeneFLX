#rep1flx DR

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
flx1DR<-flx1[c(8,1,5,13:14)]
flx1DR$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx1DR$TRAITNO<-ifelse(flx1DR$sex=="male",1,2)
flx1DR$sex<-ifelse(flx1DR$sex=="male",1,2)
#Reorder data
flx1DR<-flx1DR[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx1DR$round.sec <- as.numeric(flx1DR$round.sec)
#drop unused levels
flx1DR<-droplevels(flx1DR) 
str(flx1DR)
#round: 4 levels
length(unique(flx1DR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1DR$animal<-gsub("M","1",flx1DR$animal)
flx1DR$animal<-gsub("F","2",flx1DR$animal)
flx1DR$animal<-gsub("m","3",flx1DR$animal)
flx1DR$animal<-gsub("f","4",flx1DR$animal)
flx1DR$animal<-gsub("-","",flx1DR$animal)
#Sort by animal ID
flx1DR <- flx1DR[order(flx1DR$animal),]

#Change name of round.sec to NSEC
names(flx1DR)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1DR")
write.table(flx1DR, file="flx1DR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1DR")
write.table(ped.flx1, file="flx1DR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1flx.DR <- readRDS("rep1flx.DR")
posterior.mode(rep1flx.DR$VCV)
#traitmDR:traitmDR.animal    traitfDR:traitmDR.animal 
#2.9596462                   3.5819607 
#traitmDR:traitfDR.animal    traitfDR:traitfDR.animal 
#3.5819607                   7.4753616 
#traitmDR:traitmDR.round.sec traitfDR:traitmDR.round.sec 
#0.5685580                   0.5414712 
#traitmDR:traitfDR.round.sec traitfDR:traitfDR.round.sec 
#0.5414712                   0.5021076 
#traitmDR:traitmDR.units     traitfDR:traitmDR.units 
#1.1206512                  -3.2782057 
#traitmDR:traitfDR.units     traitfDR:traitfDR.units 
#-3.2782057                  10.5388849 

#Approximate CIs
#Va mDR
2.79035-(1.96*0.631951)
2.79035+(1.96*0.631951)
#Va fDR
9.45321-(1.96*2.85056)
9.45321+(1.96*2.85056)
#Cov DR
3.21481-(1.96*0.977542)
3.21481+(1.96*0.977542)
#rmf DR
0.626-(1.96*0.151)
0.626+(1.96*0.151)

#Log-likelihood test for covariance
#logL for model with cov: -1442.526
#logL for model without cov: -1452.236
q <- 2*(-1442.526+1452.236)
q
#q=19.42
pchisq(q, df=1, lower.tail = FALSE)
#p=1.049027e-05