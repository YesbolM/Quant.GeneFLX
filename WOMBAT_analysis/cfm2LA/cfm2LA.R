#rep2cfm LA

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
cfm2<-rep2[rep2$treatment=="CFM",] 
#  WOMBAT does not allow any missing values, let's remove them
cfm2<-na.omit(cfm2)
#Select wing size for analysis
cfm2LA<-cfm2[c(9,1,5,13:14)]
cfm2LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
cfm2LA$TRAITNO<-ifelse(cfm2LA$sex=="male",1,2)
cfm2LA$sex<-ifelse(cfm2LA$sex=="male",1,2)
#Reorder data
cfm2LA<-cfm2LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
cfm2LA$round.sec <- as.numeric(cfm2LA$round.sec)
#Drop unused levels
cfm2LA<-droplevels(cfm2LA) 
str(cfm2LA)
#round: 4 levels
length(unique(cfm2LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
cfm2LA$animal<-gsub("M","1",cfm2LA$animal)
cfm2LA$animal<-gsub("F","2",cfm2LA$animal)
cfm2LA$animal<-gsub("m","3",cfm2LA$animal)
cfm2LA$animal<-gsub("f","4",cfm2LA$animal)
cfm2LA$animal<-gsub("-","",cfm2LA$animal)
#Sort by animal ID
cfm2LA <- cfm2LA[order(cfm2LA$animal),]

#Change name of round.sec to NSEC
names(cfm2LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2LA")
write.table(cfm2LA, file="cfm2LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cfm2 <- read.csv("rep2cfm.p.csv")
head(ped.cfm2)
str(ped.cfm2)
#Change ids into numeric values
ped.cfm2$id<-gsub("M","1",ped.cfm2$id)
ped.cfm2$id<-gsub("F","2",ped.cfm2$id)
ped.cfm2$id<-gsub("m","3",ped.cfm2$id)
ped.cfm2$id<-gsub("f","4",ped.cfm2$id)
ped.cfm2$id<-gsub("-","",ped.cfm2$id)
ped.cfm2$FATHER<-gsub("M","1",ped.cfm2$FATHER)
ped.cfm2$FATHER[is.na(ped.cfm2$FATHER)] <- 0
ped.cfm2$MOTHER<-gsub("F","2",ped.cfm2$MOTHER)
ped.cfm2$MOTHER[is.na(ped.cfm2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2LA")
write.table(ped.cfm2, file="cfm2LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2cfm.LA <- readRDS("rep2cfm.LA")
posterior.mode(rep2cfm.LA$VCV)
#traitmLA.active:traitmLA.active.animal    traitfLA.active:traitmLA.active.animal 
#0.37019569                                0.31716549 
#traitmLA.active:traitfLA.active.animal    traitfLA.active:traitfLA.active.animal 
#0.31716549                                0.20883009 
#traitmLA.active:traitmLA.active.round.sec traitfLA.active:traitmLA.active.round.sec 
#0.11970081                                0.03062276 
#traitmLA.active:traitfLA.active.round.sec traitfLA.active:traitfLA.active.round.sec 
#0.03062276                                0.24382472 
#traitmLA.active:traitmLA.active.units     traitfLA.active:traitmLA.active.units 
#0.66630383                               -0.67373984 
#traitmLA.active:traitfLA.active.units     traitfLA.active:traitfLA.active.units 
#-0.67373984                                1.08795363 

#Approximate CIs
#Va mLA
2.54377-(1.96*1.10483)
2.54377+(1.96*1.10483)
#Va fLA
1.45261-(1.96*0.799420)
1.45261+(1.96*0.799420)
#Cov LA
1.43402-(1.96*0.638585)
1.43402+(1.96*0.638585)
#rmf LA
0.746-(1.96*0.305)
0.746+(1.96*0.305)

#Log-likelihood test for covariance
#logL for model with cov: -1474.100
#logL for model without cov: -1479.454
q <- 2*(-1474.100+1479.454)
q
#q=10.708
pchisq(q, df=1, lower.tail = FALSE)
#p=0.001066733
