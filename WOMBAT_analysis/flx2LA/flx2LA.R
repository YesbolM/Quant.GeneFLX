#rep2flx LA

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
flx2<-rep2[rep2$treatment=="FLX",] 
#  WOMBAT does not allow any missing values, let's remove them
flx2<-na.omit(flx2)
#Select wing size for analysis
flx2LA<-flx2[c(9,1,5,13:14)]
flx2LA$TRAITNO<-NA
#Code trait by sex, Males = 1, Females = 2
flx2LA$TRAITNO<-ifelse(flx2LA$sex=="male",1,2)
flx2LA$sex<-ifelse(flx2LA$sex=="male",1,2)
#Reorder data
flx2LA<-flx2LA[c(6,1:5)]
#Make round.sec into integers by taking factor codes
flx2LA$round.sec <- as.numeric(flx2LA$round.sec)
#Drop unused levels
flx2LA<-droplevels(flx2LA) 
str(flx2LA)
#round: 4 levels
length(unique(flx2LA$round.sec))
#round.sec: 26 levels

#Change ids into numeric values
flx2LA$animal<-gsub("M","1",flx2LA$animal)
flx2LA$animal<-gsub("F","2",flx2LA$animal)
flx2LA$animal<-gsub("m","3",flx2LA$animal)
flx2LA$animal<-gsub("f","4",flx2LA$animal)
flx2LA$animal<-gsub("-","",flx2LA$animal)
#Sort by animal ID
flx2LA <- flx2LA[order(flx2LA$animal),]

#Change name of round.sec to NSEC
names(flx2LA)[5]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2LA")
write.table(flx2LA, file="flx2LA_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.flx2 <- read.csv("rep2flx.p.csv")
head(ped.flx2)
str(ped.flx2)
#Change ids into numeric values
ped.flx2$id<-gsub("M","1",ped.flx2$id)
ped.flx2$id<-gsub("F","2",ped.flx2$id)
ped.flx2$id<-gsub("m","3",ped.flx2$id)
ped.flx2$id<-gsub("f","4",ped.flx2$id)
ped.flx2$id<-gsub("-","",ped.flx2$id)
ped.flx2$FATHER<-gsub("M","1",ped.flx2$FATHER)
ped.flx2$FATHER[is.na(ped.flx2$FATHER)] <- 0
ped.flx2$MOTHER<-gsub("F","2",ped.flx2$MOTHER)
ped.flx2$MOTHER[is.na(ped.flx2$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2LA")
write.table(ped.flx2, file="flx2LA_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\bivariate models")
rep2flx.LA <- readRDS("rep2flx.LA")
posterior.mode(rep2flx.LA$VCV)
#traitmLA.active:traitmLA.active.animal    traitfLA.active:traitmLA.active.animal 
#0.0031088328                              0.0007583186 
#traitmLA.active:traitfLA.active.animal    traitfLA.active:traitfLA.active.animal 
#0.0007583186                              0.0066203771 
#traitmLA.active:traitmLA.active.round.sec traitfLA.active:traitmLA.active.round.sec 
#0.1170349467                              0.0925386125 
#traitmLA.active:traitfLA.active.round.sec traitfLA.active:traitfLA.active.round.sec 
#0.0925386125                              0.2299885250 
#traitmLA.active:traitmLA.active.units     traitfLA.active:traitmLA.active.units 
#0.9386583792                              1.1175800662 
#traitmLA.active:traitfLA.active.units     traitfLA.active:traitfLA.active.units 
#1.1175800662                              1.3442842929 

#Approximate CIs
#Va mLA
0.142849-(1.96*0.744552)
0.142849+(1.96*0.744552)
#Va fLA
0.714140-(1.96*0.536977)
0.714140+(1.96*0.536977)
#Cov LA
0.152659-(1.96*0.419527)
0.152659+(1.96*0.419527)
#rmf LA
0.478-(1.96*0.720)
0.478+(1.96*0.720)

#Log-likelihood test for covariance
#logL for model with cov: -1399.921
#logL for model without cov: -1401.148
q <- 2*(-1399.921+1401.148)
q
#q=2.454
pchisq(q, df=1, lower.tail = FALSE)
#p=0.1172258