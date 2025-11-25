#rep1flx mLA

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
flx1m<-flx1[flx1$sex=="male",]
str(flx1m)
#Select LA for analysis
#We need to include animal, LA, round, and round.sec
flx1mLA<-flx1m[,c(1,9,13:14)]
str(flx1mLA)
#Make round.sec into integers by taking factor codes
flx1mLA$round.sec <- as.numeric(flx1mLA$round.sec)
#drop unused levels
flx1mLA<-droplevels(flx1mLA) 
str(flx1mLA)
#round: 4 levels
length(unique(flx1mLA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1mLA$animal<-gsub("M","1",flx1mLA$animal)
flx1mLA$animal<-gsub("F","2",flx1mLA$animal)
flx1mLA$animal<-gsub("m","3",flx1mLA$animal)
flx1mLA$animal<-gsub("f","4",flx1mLA$animal)
flx1mLA$animal<-gsub("-","",flx1mLA$animal)

#Change name of round.sec to NSEC
names(flx1mLA)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mLA")
write.table(flx1mLA, file="flx1mLA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1mLA_null <- flx1mLA[order(flx1mLA$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mLA\\null model")
write.table(flx1mLA_null, file="flx1mLA_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mLA")
write.table(ped.flx1, file="flx1mLA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXmLA <- readRDS("rep1FLX.mLA")
library(MCMCglmm)
posterior.mode(rep1FLXmLA$VCV)
#animal    round.sec        units 
#0.705715007 0.001879763 0.526334617

# Testing significance of random effects
#logL for model with animal: -768.889
#logL for model without animal: -776.781
q <- 2*(-768.889+776.781)
q
#q=15.784
pchisq(q, df=1, lower.tail = FALSE)
#p=7.10004e-05
