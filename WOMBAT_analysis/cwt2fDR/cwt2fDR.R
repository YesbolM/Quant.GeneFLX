#rep2cwt fDR

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
cwt2f<-cwt2[cwt2$sex=="female",]
str(cwt2f)
#Select DR for analysis
#We need to include animal, DR, round, and round.sec
cwt2fDR<-cwt2f[,c(1,8,13:14)]
str(cwt2fDR)
#Make round.sec into integers by taking factor codes
cwt2fDR$round.sec <- as.numeric(cwt2fDR$round.sec)
#Drop unused levels
cwt2fDR<-droplevels(cwt2fDR) 
str(cwt2fDR)
#round: 4 levels
length(unique(cwt2fDR$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cwt2fDR$animal<-gsub("M","1",cwt2fDR$animal)
cwt2fDR$animal<-gsub("F","2",cwt2fDR$animal)
cwt2fDR$animal<-gsub("m","3",cwt2fDR$animal)
cwt2fDR$animal<-gsub("f","4",cwt2fDR$animal)
cwt2fDR$animal<-gsub("-","",cwt2fDR$animal)

#Change name of round.sec to NSEC
names(cwt2fDR)[4]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2fDR")
write.table(cwt2fDR, file="cwt2fDR_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt2 <- read.csv("rep2CWT.p.csv")
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2fDR")
write.table(ped.cwt2, file="cwt2fDR_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\univariate models")
rep2CWTfDR <- readRDS("rep2CWT.fDR")
library(MCMCglmm)
posterior.mode(rep2CWTfDR$VCV)
#animal    round.sec        units 
#13.86198661  0.01162422 13.73410839

# Testing significance of random effects
#logL for model with animal: -962.753
#logL for model without animal: -972.502
q <- 2*(-962.753+972.502)
q
#q=19.498
pchisq(q, df=1, lower.tail = FALSE)
#p=1.007051e-05
