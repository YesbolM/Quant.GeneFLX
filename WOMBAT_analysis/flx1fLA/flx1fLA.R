#rep1flx fLA

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
flx1f<-flx1[flx1$sex=="female",]
str(flx1f)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
flx1fLA<-flx1f[,c(1,9,13:14)]
str(flx1fLA)
#Make round.sec into integers by taking factor codes
flx1fLA$round.sec <- as.numeric(flx1fLA$round.sec)
#drop unused levels
flx1fLA<-droplevels(flx1fLA) 
str(flx1fLA)
#round: 4 levels
length(unique(flx1fLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1fLA$animal<-gsub("M","1",flx1fLA$animal)
flx1fLA$animal<-gsub("F","2",flx1fLA$animal)
flx1fLA$animal<-gsub("m","3",flx1fLA$animal)
flx1fLA$animal<-gsub("f","4",flx1fLA$animal)
flx1fLA$animal<-gsub("-","",flx1fLA$animal)

#Change name of round.sec to NSEC
names(flx1fLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fLA")
write.table(flx1fLA, file="flx1fLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1fLA_null <- flx1fLA[order(flx1fLA$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fLA\\null model")
write.table(flx1fLA_null, file="flx1fLA_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fLA")
write.table(ped.flx1, file="flx1fLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXfLA <- readRDS("rep1FLX.fLA")
library(MCMCglmm)
posterior.mode(rep1FLXfLA$VCV)
#animal    round.sec        units 
#0.004430139 0.002441968 0.710829756

# Testing significance of random effects
#logL for model with animal: -677.033
#logL for model without animal: -679.064
q <- 2*(-677.033+679.064)
q
#q=4.062
pchisq(q, df=1, lower.tail = FALSE)
#p=0.04385853