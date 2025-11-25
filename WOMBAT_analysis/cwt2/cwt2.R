#rep2cwt

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
#Select traits for analysis
cwt2<-cwt2[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
cwt2<-na.omit(cwt2)

#Mean standardization of values
cwt2$mLA<-ifelse(cwt2$sex=="male",cwt2$LA_active,NA)
meanLA.m<-mean(cwt2$mLA, na.rm=TRUE)
cwt2$LA.m<-cwt2$mLA/meanLA.m

cwt2$fLA.<-ifelse(cwt2$sex=="female",cwt2$LA_active,NA)
meanLA.f<-mean(cwt2$fLA, na.rm=TRUE)
cwt2$LA.f<-cwt2$fLA/meanLA.f

cwt2$mWS<-ifelse(cwt2$sex=="male",cwt2$WS,NA)
meanWs.m<-mean(cwt2$mWS,na.rm=TRUE)
cwt2$meanWs.m<-cwt2$mWS/meanWs.m

cwt2$fWS<-ifelse(cwt2$sex=="female",cwt2$WS,NA)
meanWS.f<-mean(cwt2$fWS,na.rm=TRUE)
cwt2$meanWs.f<-cwt2$fWS/meanWS.f

cwt2$mDR<-ifelse(cwt2$sex=="male",cwt2$DR ,NA)
meanDR.m<-mean(cwt2$mDR,na.rm=TRUE)
cwt2$meanDr.m<-cwt2$mDR/meanDR.m

cwt2$fDR<-ifelse(cwt2$sex=="female",cwt2$DR,NA)
meanDR.f<-mean(cwt2$fDR,na.rm=TRUE)
cwt2$meanDr.f<-cwt2$fDR/meanDR.f

cwt2$WS<-ifelse(cwt2$sex=="male",cwt2$meanWs.m,cwt2$meanWs.f)
cwt2$DR<-ifelse(cwt2$sex=="male",cwt2$meanDr.m,cwt2$meanDr.f)
cwt2$LA_active<-ifelse(cwt2$sex=="male",cwt2$LA.m,cwt2$LA.f)

cwt2<-cwt2[,1:7]

#Set up trait numbers
cwt2$TRAITNO<-NA
#Select wing size data
cwt2.WS<-cwt2[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
cwt2.WS$TRAITNO<-ifelse(cwt2.WS$sex=="male",1,4)
#Rename trait column
names(cwt2.WS)[2] = "trait"
#Select DR data
cwt2.DR<-cwt2[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
cwt2.DR$TRAITNO<-ifelse(cwt2.DR$sex=="male",2,5)
#Rename trait column
names(cwt2.DR)[2] = "trait"
#Select LA data
cwt2.LA<-cwt2[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
cwt2.LA$TRAITNO<-ifelse(cwt2.LA$sex=="male",3,6)
#Rename trait column
names(cwt2.LA)[2] = "trait"
#Put together into single dataset
cwt2.dat<-rbind(cwt2.WS,cwt2.DR,cwt2.LA)
#Code trait by sex, Males = 1, Females = 2
cwt2.dat$sex<-ifelse(cwt2.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
cwt2.dat$round.sec <- as.numeric(cwt2.dat$round.sec)
#drop unused levels
cwt2.dat<-droplevels(cwt2.dat) 
str(cwt2.dat)
#round: 4 levels
length(unique(cwt2.dat$round.sec))
#round.sec: 26 levels
length(unique(cwt2.dat$animal))
#animal: 950 levels

#Change ids into numeric values
cwt2.dat$animal<-gsub("M","1",cwt2.dat$animal)
cwt2.dat$animal<-gsub("F","2",cwt2.dat$animal)
cwt2.dat$animal<-gsub("m","3",cwt2.dat$animal)
cwt2.dat$animal<-gsub("f","4",cwt2.dat$animal)
cwt2.dat$animal<-gsub("-","",cwt2.dat$animal)
#Sort by animal ID and then trait number
cwt2.dat <- cwt2.dat[order(cwt2.dat$animal, cwt2.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(cwt2.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2")
write.table(cwt2.dat, file="cwt2_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt2 <- read.csv("rep2cwt.p.csv")
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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt2")
write.table(ped.cwt2, file="cwt2_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep2cwt <- readRDS("rep2cwt")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep2cwt$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#0.0015520546                               0.0007655915 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#0.0045973883                               0.0008979807 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#-0.0017876461                              -0.0011861102 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#0.0473277903                              -0.0077838942 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#0.0017520805                               0.0484448274 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#-0.0425454541                               0.1864862481 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#0.0038963738                              -0.0230731653 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#0.1018426593                               0.0015680607 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#0.0032564412                              -0.0014341119 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#0.0850614451                              -0.0316081403 
#traitmeanActive.f:traitmeanActive.f.animal 
#0.0727359100 

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#5.341646e-04                               -6.186609e-04 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#3.767749e-03                               -1.018960e-04 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-6.130743e-05                                3.718164e-03 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#1.455358e-02                                3.546568e-03 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#-3.725988e-04                                7.457031e-04 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#-1.086728e-02                                5.738630e-02 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#-1.921569e-03                               -1.008954e-03 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#7.705456e-02                                4.406747e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#5.620309e-05                               -1.918620e-03 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#4.836603e-03                                9.233851e-03 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#2.561160e-01 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
3.846099e-04                              1.424101e-04 
traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
-1.490724e-03                             -1.241547e-05 
traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
-7.219962e-04                              1.481181e-03 
traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
1.132006e-01                              9.829397e-03 
traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
6.008595e-04                              7.089292e-02 
traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
1.365599e-01                              4.662575e-01 
traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
-8.460435e-04                              1.497366e-02 
traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
-5.247137e-01                              2.897487e-04 
traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
-3.571012e-04                              6.938585e-04 
traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
9.193411e-02                             -3.512146e-03 
traitmeanActive.f:traitmeanActive.f.units 
8.214362e-01 