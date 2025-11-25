#rep1flx mDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
flx1mDR<-flx1m[,c(1,8,13:14)]
str(flx1mDR)
#Make round.sec into integers by taking factor codes
flx1mDR$round.sec <- as.numeric(flx1mDR$round.sec)
#Drop unused levels
flx1mDR<-droplevels(flx1mDR) 
str(flx1mDR)
#round: 4 levels
length(unique(flx1mDR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1mDR$animal<-gsub("M","1",flx1mDR$animal)
flx1mDR$animal<-gsub("F","2",flx1mDR$animal)
flx1mDR$animal<-gsub("m","3",flx1mDR$animal)
flx1mDR$animal<-gsub("f","4",flx1mDR$animal)
flx1mDR$animal<-gsub("-","",flx1mDR$animal)

#Change name of round.sec to NSEC
names(flx1mDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mDR")
write.table(flx1mDR, file="flx1mDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1mDR_null <- flx1mDR[order(flx1mDR$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mDR\\null model")
write.table(flx1mDR_null, file="flx1mDR_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1mDR")
write.table(ped.flx1, file="flx1mDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXmDR <- readRDS("rep1FLX.mDR")
library(MCMCglmm)
posterior.mode(rep1FLXmDR$VCV)
#animal    round.sec        units 
#4.295603856 0.007182980 0.008435656

# Testing significance of random effects
#logL for model with animal: -554.321
#logL for model without animal: -576.047
q <- 2*(-554.321+576.047)
q
#q=43.452
pchisq(q, df=1, lower.tail = FALSE)
#p=4.344873e-11
