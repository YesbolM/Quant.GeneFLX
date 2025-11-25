#rep1flx fWS

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
flx1fWS<-flx1f[,c(1,7,13:14)]
str(flx1fWS)
#Make round.sec into integers by taking factor codes
flx1fWS$round.sec <- as.numeric(flx1fWS$round.sec)
#drop unused levels
flx1fWS<-droplevels(flx1fWS) 
str(flx1fWS)
#round: 4 levels
length(unique(flx1fWS$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx1fWS$animal<-gsub("M","1",flx1fWS$animal)
flx1fWS$animal<-gsub("F","2",flx1fWS$animal)
flx1fWS$animal<-gsub("m","3",flx1fWS$animal)
flx1fWS$animal<-gsub("f","4",flx1fWS$animal)
flx1fWS$animal<-gsub("-","",flx1fWS$animal)

#Change name of round.sec to NSEC
names(flx1fWS)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fWS")
write.table(flx1fWS, file="flx1fWS_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#sort by NSEC for null model analysis
flx1fWS_null <- flx1fWS[order(flx1fWS$NSEC),]
#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fWS\\null model")
write.table(flx1fWS_null, file="flx1fWS_datnull.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1fWS")
write.table(ped.flx1, file="flx1fWS_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep1FLXfWS <- readRDS("rep1FLX.fWS")
library(MCMCglmm)
posterior.mode(rep1FLXfWS$VCV)
#animal    round.sec        units 
#0.0022678066 0.0003423411 0.0004275114

# Testing significance of random effects
#logL for model with animal: 1135.268
#logL for model without animal: 1111.656
q <- 2*(1135.268-1111.656)
q
#q=47.224
pchisq(q, df=1, lower.tail = FALSE)
#p=6.332036e-12