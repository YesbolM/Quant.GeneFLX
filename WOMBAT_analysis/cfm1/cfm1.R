#rep1cfm

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
cfm1<-rep1[rep1$treatment=="CFM",] 
#Select traits for analysis
cfm1<-cfm1[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
cfm1<-na.omit(cfm1)

#Mean standardization of values
cfm1$mLA<-ifelse(cfm1$sex=="male",cfm1$LA_active,NA)
meanLA.m<-mean(cfm1$mLA, na.rm=TRUE)
cfm1$LA.m<-cfm1$mLA/meanLA.m

cfm1$fLA.<-ifelse(cfm1$sex=="female",cfm1$LA_active,NA)
meanLA.f<-mean(cfm1$fLA, na.rm=TRUE)
cfm1$LA.f<-cfm1$fLA/meanLA.f

cfm1$mWS<-ifelse(cfm1$sex=="male",cfm1$WS,NA)
meanWs.m<-mean(cfm1$mWS,na.rm=TRUE)
cfm1$meanWs.m<-cfm1$mWS/meanWs.m

cfm1$fWS<-ifelse(cfm1$sex=="female",cfm1$WS,NA)
meanWS.f<-mean(cfm1$fWS,na.rm=TRUE)
cfm1$meanWs.f<-cfm1$fWS/meanWS.f

cfm1$mDR<-ifelse(cfm1$sex=="male",cfm1$DR ,NA)
meanDR.m<-mean(cfm1$mDR,na.rm=TRUE)
cfm1$meanDr.m<-cfm1$mDR/meanDR.m

cfm1$fDR<-ifelse(cfm1$sex=="female",cfm1$DR,NA)
meanDR.f<-mean(cfm1$fDR,na.rm=TRUE)
cfm1$meanDr.f<-cfm1$fDR/meanDR.f

cfm1$WS<-ifelse(cfm1$sex=="male",cfm1$meanWs.m,cfm1$meanWs.f)
cfm1$DR<-ifelse(cfm1$sex=="male",cfm1$meanDr.m,cfm1$meanDr.f)
cfm1$LA_active<-ifelse(cfm1$sex=="male",cfm1$LA.m,cfm1$LA.f)

cfm1<-cfm1[,1:7]

#Set up trait numbers
cfm1$TRAITNO<-NA
#Select wing size data
cfm1.WS<-cfm1[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
cfm1.WS$TRAITNO<-ifelse(cfm1.WS$sex=="male",1,4)
#Rename trait column
names(cfm1.WS)[2] = "trait"
#Select DR data
cfm1.DR<-cfm1[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
cfm1.DR$TRAITNO<-ifelse(cfm1.DR$sex=="male",2,5)
#Rename trait column
names(cfm1.DR)[2] = "trait"
#Select LA data
cfm1.LA<-cfm1[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
cfm1.LA$TRAITNO<-ifelse(cfm1.LA$sex=="male",3,6)
#Rename trait column
names(cfm1.LA)[2] = "trait"
#Put together into single dataset
cfm1.dat<-rbind(cfm1.WS,cfm1.DR,cfm1.LA)
#Code trait by sex, Males = 1, Females = 2
cfm1.dat$sex<-ifelse(cfm1.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
cfm1.dat$round.sec <- as.numeric(cfm1.dat$round.sec)
#drop unused levels
cfm1.dat<-droplevels(cfm1.dat) 
str(cfm1.dat)
#round: 4 levels
length(unique(cfm1.dat$round.sec))
#round.sec: 27 levels
length(unique(cfm1.dat$animal))
#animal: 947 levels

#Change ids into numeric values
cfm1.dat$animal<-gsub("M","1",cfm1.dat$animal)
cfm1.dat$animal<-gsub("F","2",cfm1.dat$animal)
cfm1.dat$animal<-gsub("m","3",cfm1.dat$animal)
cfm1.dat$animal<-gsub("f","4",cfm1.dat$animal)
cfm1.dat$animal<-gsub("-","",cfm1.dat$animal)
#Sort by animal ID and then trait number
cfm1.dat <- cfm1.dat[order(cfm1.dat$animal, cfm1.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(cfm1.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1")
write.table(cfm1.dat, file="cfm1_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cfm1 <- read.csv("rep1cfm.p.csv")
head(ped.cfm1)
str(ped.cfm1)
#Change ids into numeric values
ped.cfm1$id<-gsub("M","1",ped.cfm1$id)
ped.cfm1$id<-gsub("F","2",ped.cfm1$id)
ped.cfm1$id<-gsub("m","3",ped.cfm1$id)
ped.cfm1$id<-gsub("f","4",ped.cfm1$id)
ped.cfm1$id<-gsub("-","",ped.cfm1$id)
ped.cfm1$FATHER<-gsub("M","1",ped.cfm1$FATHER)
ped.cfm1$FATHER[is.na(ped.cfm1$FATHER)] <- 0
ped.cfm1$MOTHER<-gsub("F","2",ped.cfm1$MOTHER)
ped.cfm1$MOTHER[is.na(ped.cfm1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm1")
write.table(ped.cfm1, file="cfm1_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep1cfm <- readRDS("rep1cfm")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep1cfm$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#0.0009973198                               0.0004690675 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#0.0023926187                               0.0006452605 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#0.0012393469                               0.0006616731 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#0.0680139607                              -0.0320732536 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#0.0005434323                               0.0341247278 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#-0.0341828167                               0.4437949512 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#0.0007671578                               0.0582813672 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#0.2392783852                               0.0009260752 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#0.0001205472                               0.0010288189 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#0.0556632163                               0.0108103220 
#traitmeanActive.f:traitmeanActive.f.animal 
#0.1852414166 

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#2.466591e-04                               -3.718968e-04 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#1.603535e-04                                3.669512e-05 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-6.441845e-04                               -1.822410e-04 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#2.699051e-02                                6.252646e-03 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#-4.329081e-04                                2.676312e-02 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#1.283780e-02                                8.004794e-03 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#-1.312689e-04                               -6.421362e-03 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#6.165614e-03                                2.509296e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#-1.611736e-04                               -2.220275e-04 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#3.976145e-02                                1.294439e-02 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#1.096133e-02 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
#traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
#5.015979e-04                              4.540509e-04 
#traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
#9.969142e-04                              2.562255e-05 
#traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
#-3.071637e-04                             -1.372138e-05 
#traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
#5.450666e-02                              1.758122e-03 
#traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
#-5.067446e-04                             -3.111102e-02 
#traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
#8.180119e-02                              7.186831e-01 
#traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
#-1.526627e-05                             -2.200724e-01 
#traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
#9.135398e-01                              3.791899e-04 
#traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
#2.905022e-04                              4.897957e-03 
#traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
#1.128801e-01                             -1.034826e-01 
#traitmeanActive.f:traitmeanActive.f.units 
#1.169746e+00 