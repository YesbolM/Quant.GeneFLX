#rep1flx

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
#Select traits for analysis
flx1<-flx1[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
flx1<-na.omit(flx1)

#Mean standardization of values
flx1$mLA<-ifelse(flx1$sex=="male",flx1$LA_active,NA)
meanLA.m<-mean(flx1$mLA, na.rm=TRUE)
flx1$LA.m<-flx1$mLA/meanLA.m

flx1$fLA.<-ifelse(flx1$sex=="female",flx1$LA_active,NA)
meanLA.f<-mean(flx1$fLA, na.rm=TRUE)
flx1$LA.f<-flx1$fLA/meanLA.f

flx1$mWS<-ifelse(flx1$sex=="male",flx1$WS,NA)
meanWs.m<-mean(flx1$mWS,na.rm=TRUE)
flx1$meanWs.m<-flx1$mWS/meanWs.m

flx1$fWS<-ifelse(flx1$sex=="female",flx1$WS,NA)
meanWS.f<-mean(flx1$fWS,na.rm=TRUE)
flx1$meanWs.f<-flx1$fWS/meanWS.f

flx1$mDR<-ifelse(flx1$sex=="male",flx1$DR ,NA)
meanDR.m<-mean(flx1$mDR,na.rm=TRUE)
flx1$meanDr.m<-flx1$mDR/meanDR.m

flx1$fDR<-ifelse(flx1$sex=="female",flx1$DR,NA)
meanDR.f<-mean(flx1$fDR,na.rm=TRUE)
flx1$meanDr.f<-flx1$fDR/meanDR.f

flx1$WS<-ifelse(flx1$sex=="male",flx1$meanWs.m,flx1$meanWs.f)
flx1$DR<-ifelse(flx1$sex=="male",flx1$meanDr.m,flx1$meanDr.f)
flx1$LA_active<-ifelse(flx1$sex=="male",flx1$LA.m,flx1$LA.f)

flx1<-flx1[,1:7]

#Set up trait numbers
flx1$TRAITNO<-NA
#Select wing size data
flx1.WS<-flx1[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
flx1.WS$TRAITNO<-ifelse(flx1.WS$sex=="male",1,4)
#Rename trait column
names(flx1.WS)[2] = "trait"
#Select DR data
flx1.DR<-flx1[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
flx1.DR$TRAITNO<-ifelse(flx1.DR$sex=="male",2,5)
#Rename trait column
names(flx1.DR)[2] = "trait"
#Select LA data
flx1.LA<-flx1[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
flx1.LA$TRAITNO<-ifelse(flx1.LA$sex=="male",3,6)
#Rename trait column
names(flx1.LA)[2] = "trait"
#Put together into single dataset
flx1.dat<-rbind(flx1.WS,flx1.DR,flx1.LA)
#Code trait by sex, Males = 1, Females = 2
flx1.dat$sex<-ifelse(flx1.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
flx1.dat$round.sec <- as.numeric(flx1.dat$round.sec)
#drop unused levels
flx1.dat<-droplevels(flx1.dat) 
str(flx1.dat)
#round: 4 levels
length(unique(flx1.dat$round.sec))
#round.sec: 26 levels
length(unique(flx1.dat$animal))
#animal: 945 levels

#Change ids into numeric values
flx1.dat$animal<-gsub("M","1",flx1.dat$animal)
flx1.dat$animal<-gsub("F","2",flx1.dat$animal)
flx1.dat$animal<-gsub("m","3",flx1.dat$animal)
flx1.dat$animal<-gsub("f","4",flx1.dat$animal)
flx1.dat$animal<-gsub("-","",flx1.dat$animal)
#Sort by animal ID and then trait number
flx1.dat <- flx1.dat[order(flx1.dat$animal, flx1.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(flx1.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1")
write.table(flx1.dat, file="flx1_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx1")
write.table(ped.flx1, file="flx1_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep1flx <- readRDS("rep1flx")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep1flx$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#0.0009714517                               0.0021352285 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#0.0038383550                               0.0007137548 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#0.0015635071                               0.0021532318 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#0.0675469158                               0.0234716991 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#0.0007763721                               0.0400512211 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#0.0661091755                               0.3765542478 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#0.0015782511                              -0.0137369018 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#-0.0167293208                               0.0009277351 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#0.0008681601                               0.0016023233 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#0.0666567452                               0.0476071158 
#traitmeanActive.f:traitmeanActive.f.animal 
#0.1348304246 

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#2.273910e-04                               -5.144454e-04 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#7.345749e-04                                6.880101e-05 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-5.082861e-05                                4.740871e-04 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#2.698518e-02                               -4.444999e-02 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#-6.760750e-04                                1.233885e-02 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#-2.927563e-02                                6.814820e-02 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#1.361452e-03                               -1.705246e-02 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#1.257523e-01                                2.676005e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#-3.507076e-04                                1.235554e-03 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#9.970048e-03                               -1.323562e-02 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#1.469418e-01 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
#traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
#2.739703e-04                             -4.261326e-04 
#traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
#-9.986395e-04                             -9.222264e-06 
#traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
#-3.684390e-04                             -1.703470e-03 
#traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
#3.916256e-02                             -2.881221e-02 
#traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
#1.510667e-04                              4.813272e-02 
#traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
#5.640863e-02                              7.022670e-01 
#traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
#-4.690786e-05                             -6.425878e-02 
#traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
#-1.117370e+00                              2.386179e-04 
#traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
#3.713422e-04                              1.086034e-03 
#traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
#7.889975e-02                             -5.760093e-02 
#traitmeanActive.f:traitmeanActive.f.units 
#1.789611e+00 
