#rep2flx

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
#Select traits for analysis
flx2<-flx2[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
flx2<-na.omit(flx2)

#Mean standardization of values
flx2$mLA<-ifelse(flx2$sex=="male",flx2$LA_active,NA)
meanLA.m<-mean(flx2$mLA, na.rm=TRUE)
flx2$LA.m<-flx2$mLA/meanLA.m

flx2$fLA.<-ifelse(flx2$sex=="female",flx2$LA_active,NA)
meanLA.f<-mean(flx2$fLA, na.rm=TRUE)
flx2$LA.f<-flx2$fLA/meanLA.f

flx2$mWS<-ifelse(flx2$sex=="male",flx2$WS,NA)
meanWs.m<-mean(flx2$mWS,na.rm=TRUE)
flx2$meanWs.m<-flx2$mWS/meanWs.m

flx2$fWS<-ifelse(flx2$sex=="female",flx2$WS,NA)
meanWS.f<-mean(flx2$fWS,na.rm=TRUE)
flx2$meanWs.f<-flx2$fWS/meanWS.f

flx2$mDR<-ifelse(flx2$sex=="male",flx2$DR ,NA)
meanDR.m<-mean(flx2$mDR,na.rm=TRUE)
flx2$meanDr.m<-flx2$mDR/meanDR.m

flx2$fDR<-ifelse(flx2$sex=="female",flx2$DR,NA)
meanDR.f<-mean(flx2$fDR,na.rm=TRUE)
flx2$meanDr.f<-flx2$fDR/meanDR.f

flx2$WS<-ifelse(flx2$sex=="male",flx2$meanWs.m,flx2$meanWs.f)
flx2$DR<-ifelse(flx2$sex=="male",flx2$meanDr.m,flx2$meanDr.f)
flx2$LA_active<-ifelse(flx2$sex=="male",flx2$LA.m,flx2$LA.f)

flx2<-flx2[,1:7]

#Set up trait numbers
flx2$TRAITNO<-NA
#Select wing size data
flx2.WS<-flx2[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
flx2.WS$TRAITNO<-ifelse(flx2.WS$sex=="male",1,4)
#Rename trait column
names(flx2.WS)[2] = "trait"
#Select DR data
flx2.DR<-flx2[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
flx2.DR$TRAITNO<-ifelse(flx2.DR$sex=="male",2,5)
#Rename trait column
names(flx2.DR)[2] = "trait"
#Select LA data
flx2.LA<-flx2[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
flx2.LA$TRAITNO<-ifelse(flx2.LA$sex=="male",3,6)
#Rename trait column
names(flx2.LA)[2] = "trait"
#Put together into single dataset
flx2.dat<-rbind(flx2.WS,flx2.DR,flx2.LA)
#Code trait by sex, Males = 1, Females = 2
flx2.dat$sex<-ifelse(flx2.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
flx2.dat$round.sec <- as.numeric(flx2.dat$round.sec)
#drop unused levels
flx2.dat<-droplevels(flx2.dat) 
str(flx2.dat)
#round: 4 levels
length(unique(flx2.dat$round.sec))
#round.sec: 26 levels
length(unique(flx2.dat$animal))
#animal: 954 levels

#Change ids into numeric values
flx2.dat$animal<-gsub("M","1",flx2.dat$animal)
flx2.dat$animal<-gsub("F","2",flx2.dat$animal)
flx2.dat$animal<-gsub("m","3",flx2.dat$animal)
flx2.dat$animal<-gsub("f","4",flx2.dat$animal)
flx2.dat$animal<-gsub("-","",flx2.dat$animal)
#Sort by animal ID and then trait number
flx2.dat <- flx2.dat[order(flx2.dat$animal, flx2.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(flx2.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2")
write.table(flx2.dat, file="flx2_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

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
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\flx2")
write.table(ped.flx2, file="flx2_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep2flx <- readRDS("rep2flx")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep2flx$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#1.131490e-03                              -4.921477e-05 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#-2.041107e-04                               5.282327e-04 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#3.471941e-04                               7.252594e-04 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#3.661580e-02                              -1.310448e-03 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#-8.983748e-05                               3.215700e-02 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#-3.365084e-02                               9.476748e-03 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#-1.411851e-04                              -4.804578e-03 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#9.346420e-03                               9.649350e-04 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#2.881537e-04                               6.779235e-04 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#5.746240e-02                              -3.997272e-02 
#traitmeanActive.f:traitmeanActive.f.animal 
#2.207929e-02 

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#2.403992e-04                               -6.051327e-04 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#2.696640e-04                                4.178072e-05 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-3.287310e-04                               -2.443616e-04 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#5.665976e-02                                1.591976e-02 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#-7.087867e-04                                2.133931e-02 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#2.754768e-02                                7.727601e-02 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#3.969094e-04                               -2.539442e-03 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#1.494483e-02                                3.613153e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#-6.025614e-04                               -2.860629e-04 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#2.888585e-02                                3.147813e-02 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#9.649710e-02 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
#traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
#3.319258e-04                              8.808054e-04 
#traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
#2.361668e-03                              3.081148e-05 
#traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
#1.088053e-03                             -9.500723e-05 
#traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
#5.600682e-02                             -4.953405e-03 
#traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
#-1.036265e-04                              5.667775e-02 
#traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
#1.516546e-01                              1.014288e+00 
#traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
#2.435702e-03                              1.363283e-01 
#traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
#-1.136502e+00                              4.307519e-04 
#traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
#5.457380e-04                             -2.333717e-03 
#traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
#7.577411e-02                              1.056784e-02 
#traitmeanActive.f:traitmeanActive.f.units 
#1.390230e+00 
