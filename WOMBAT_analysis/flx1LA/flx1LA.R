#rep1flx LA

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
#Select trait for analysis
flx1LA<-flx1[c(9,1,5,13:14)]
flx1LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx1LA$TRAITNO<-ifelse(flx1LA$sex=="male",1,2)
flx1LA$sex<-ifelse(flx1LA$sex=="male",1,2)
#Reorder data
flx1LA<-flx1LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx1LA$round.sec <- as.numeric(flx1LA$round.sec)
#Drop unused levels
flx1LA<-droplevels(flx1LA) 
str(flx1LA)
#round: 4 levels
length(unique(flx1LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1LA$animal<-gsub("M","1",flx1LA$animal)
flx1LA$animal<-gsub("F","2",flx1LA$animal)
flx1LA$animal<-gsub("m","3",flx1LA$animal)
flx1LA$animal<-gsub("f","4",flx1LA$animal)
flx1LA$animal<-gsub("-","",flx1LA$animal)
#Sort by animal ID
flx1LA <- flx1LA[order(flx1LA$animal),]

#Change name of round.sec to NSEC
names(flx1LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1LA")
write.table(flx1LA, file="flx1LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1LA")
write.table(ped.flx1, file="flx1LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep1flx.LA <- readRDS("rep1flx.LA")
posterior.mode(rep1flx.LA$VCV)
#traitmLA.active:traitmLA.active.animal 
#0.800193022 
#traitfLA.active:traitmLA.active.animal 
#0.118447803 
#traitmLA.active:traitfLA.active.animal 
#0.118447803 
#traitfLA.active:traitfLA.active.animal 
#0.004214209 
#traitmLA.active:traitmLA.active.round.sec 
#0.004422699 
#traitfLA.active:traitmLA.active.round.sec 
#0.045459067 
#traitmLA.active:traitfLA.active.round.sec 
#0.045459067 
#traitfLA.active:traitfLA.active.round.sec 
#0.061985855 
#traitmLA.active:traitmLA.active.units 
#0.624806844 
#traitfLA.active:traitmLA.active.units 
#0.655607224 
#traitmLA.active:traitfLA.active.units 
#0.655607224 
#traitfLA.active:traitfLA.active.units 
#0.992207674 

#Approximate CIs
#Va mLA
3.23604-(1.96*1.13711)
3.23604+(1.96*1.13711)
#Va fLA
1.27573-(1.96*0.782768)
1.27573+(1.96*0.782768)
#Cov LA
0.135360-(1.96*0.643569)
0.135360+(1.96*0.643569)
#rmf LA
0.067-(1.96*0.316)
0.067+(1.96*0.316)

#Log-likelihood test for covariance
#logL for model with cov: -1445.288
#logL for model without cov: -1447.191
q <- 2*(-1445.288+1447.191)
q
#q=3.806
pchisq(q, df=1, lower.tail = FALSE)
#p=0.05106927
