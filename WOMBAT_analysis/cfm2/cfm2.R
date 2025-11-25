#rep2cfm

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
#Select traits for analysis
cfm2<-cfm2[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
cfm2<-na.omit(cfm2)

#Mean standardization of values
cfm2$mLA<-ifelse(cfm2$sex=="male",cfm2$LA_active,NA)
meanLA.m<-mean(cfm2$mLA, na.rm=TRUE)
cfm2$LA.m<-cfm2$mLA/meanLA.m

cfm2$fLA.<-ifelse(cfm2$sex=="female",cfm2$LA_active,NA)
meanLA.f<-mean(cfm2$fLA, na.rm=TRUE)
cfm2$LA.f<-cfm2$fLA/meanLA.f

cfm2$mWS<-ifelse(cfm2$sex=="male",cfm2$WS,NA)
meanWs.m<-mean(cfm2$mWS,na.rm=TRUE)
cfm2$meanWs.m<-cfm2$mWS/meanWs.m

cfm2$fWS<-ifelse(cfm2$sex=="female",cfm2$WS,NA)
meanWS.f<-mean(cfm2$fWS,na.rm=TRUE)
cfm2$meanWs.f<-cfm2$fWS/meanWS.f

cfm2$mDR<-ifelse(cfm2$sex=="male",cfm2$DR ,NA)
meanDR.m<-mean(cfm2$mDR,na.rm=TRUE)
cfm2$meanDr.m<-cfm2$mDR/meanDR.m

cfm2$fDR<-ifelse(cfm2$sex=="female",cfm2$DR,NA)
meanDR.f<-mean(cfm2$fDR,na.rm=TRUE)
cfm2$meanDr.f<-cfm2$fDR/meanDR.f

cfm2$WS<-ifelse(cfm2$sex=="male",cfm2$meanWs.m,cfm2$meanWs.f)
cfm2$DR<-ifelse(cfm2$sex=="male",cfm2$meanDr.m,cfm2$meanDr.f)
cfm2$LA_active<-ifelse(cfm2$sex=="male",cfm2$LA.m,cfm2$LA.f)

cfm2<-cfm2[,1:7]

#Set up trait numbers
cfm2$TRAITNO<-NA
#Select wing size data
cfm2.WS<-cfm2[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
cfm2.WS$TRAITNO<-ifelse(cfm2.WS$sex=="male",1,4)
#Rename trait column
names(cfm2.WS)[2] = "trait"
#Select DR data
cfm2.DR<-cfm2[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
cfm2.DR$TRAITNO<-ifelse(cfm2.DR$sex=="male",2,5)
#Rename trait column
names(cfm2.DR)[2] = "trait"
#Select LA data
cfm2.LA<-cfm2[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
cfm2.LA$TRAITNO<-ifelse(cfm2.LA$sex=="male",3,6)
#Rename trait column
names(cfm2.LA)[2] = "trait"
#Put together into single dataset
cfm2.dat<-rbind(cfm2.WS,cfm2.DR,cfm2.LA)
#Code trait by sex, Males = 1, Females = 2
cfm2.dat$sex<-ifelse(cfm2.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
cfm2.dat$round.sec <- as.numeric(cfm2.dat$round.sec)
#drop unused levels
cfm2.dat<-droplevels(cfm2.dat) 
str(cfm2.dat)
#round: 4 levels
length(unique(cfm2.dat$round.sec))
#round.sec: 27 levels
length(unique(cfm2.dat$animal))
#animal: 953 levels

#Change ids into numeric values
cfm2.dat$animal<-gsub("M","1",cfm2.dat$animal)
cfm2.dat$animal<-gsub("F","2",cfm2.dat$animal)
cfm2.dat$animal<-gsub("m","3",cfm2.dat$animal)
cfm2.dat$animal<-gsub("f","4",cfm2.dat$animal)
cfm2.dat$animal<-gsub("-","",cfm2.dat$animal)
#Sort by animal ID and then trait number
cfm2.dat <- cfm2.dat[order(cfm2.dat$animal, cfm2.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(cfm2.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2")
write.table(cfm2.dat, file="cfm2_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cfm2")
write.table(ped.cfm2, file="cfm2_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep2cfm <- readRDS("rep2cfm")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep2cfm$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#1.069756e-03                              -8.430176e-04 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#-1.170656e-04                               5.900434e-04 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#1.219509e-04                              -1.355515e-03 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#3.132515e-02                              -1.834932e-02 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#6.730295e-05                               3.766106e-02 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#2.096537e-03                               2.412923e-01 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#5.154019e-03                              -1.231539e-02 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#3.006118e-01                               9.579861e-04 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#1.151806e-03                               8.469808e-03 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#6.455471e-02                              -2.492414e-02 
#traitmeanActive.f:traitmeanActive.f.animal 
#5.786106e-01 

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#2.638866e-04                                5.030625e-05 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#1.030024e-03                                2.898538e-05 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-1.804079e-05                                5.410344e-04 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#6.067153e-03                               -1.035513e-02 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#9.819428e-05                                5.188029e-03 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#9.899943e-03                                8.444778e-02 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#5.908923e-04                               -2.137355e-02 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#-2.825546e-02                                2.440399e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#2.664967e-04                                4.383792e-04 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#1.342121e-02                                2.138112e-02 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#1.861731e-01 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
#traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
#6.019280e-04                              6.974320e-04 
#traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
#7.419532e-04                             -5.310313e-05 
#traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
#4.283215e-04                             -2.303415e-03 
#traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
#5.789557e-02                              2.372642e-02 
#traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
#3.915837e-05                             -5.442241e-02 
#traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
#2.189794e-02                              1.006441e+00 
#traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
#6.165889e-03                              4.607183e-02 
#traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
#-1.003721e+00                              6.806827e-04 
#traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
#5.979348e-04                             -4.733967e-03 
#traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
#1.056645e-01                             -3.318093e-02 
#traitmeanActive.f:traitmeanActive.f.units 
#1.335135e+00 