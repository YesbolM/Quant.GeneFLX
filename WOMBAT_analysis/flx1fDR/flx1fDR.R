#rep1flx fDR

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
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
flx1fDR<-flx1f[,c(1,8,13:14)]
str(flx1fDR)
#Make round.sec into integers by taking factor codes
flx1fDR$round.sec <- as.numeric(flx1fDR$round.sec)
#Drop unused levels
flx1fDR<-droplevels(flx1fDR) 
str(flx1fDR)
#round: 4 levels
length(unique(flx1fDR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1fDR$animal<-gsub("M","1",flx1fDR$animal)
flx1fDR$animal<-gsub("F","2",flx1fDR$animal)
flx1fDR$animal<-gsub("m","3",flx1fDR$animal)
flx1fDR$animal<-gsub("f","4",flx1fDR$animal)
flx1fDR$animal<-gsub("-","",flx1fDR$animal)

#Change name of round.sec to NSEC
names(flx1fDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fDR")
write.table(flx1fDR, file="flx1fDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1fDR_null <- flx1fDR[order(flx1fDR$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fDR\\null model")
write.table(flx1fDR_null, file="flx1fDR_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fDR")
write.table(ped.flx1, file="flx1fDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXfDR <- readRDS("rep1FLX.fDR")
library(MCMCglmm)
posterior.mode(rep1FLXfDR$VCV)
#animal    round.sec        units 
#12.045197552  0.004507341  9.451230947 

# Testing significance of random effects
#logL for model with animal: -901.293
#logL for model without animal: -911.275
q <- 2*(-901.293+911.275)
q
#q=19.964
pchisq(q, df=1, lower.tail = FALSE)
#p=7.891401e-06
