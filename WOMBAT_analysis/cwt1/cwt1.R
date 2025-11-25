#rep1cwt

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
cwt1<-rep1[rep1$treatment=="CWT",] 
#Select traits for analysis
cwt1<-cwt1[c(7,8,9,1,5,13:14)]
#  WOMBAT does not allow any missing values, let's remove them
cwt1<-na.omit(cwt1)

#Mean standardization of values
cwt1$mLA<-ifelse(cwt1$sex=="male",cwt1$LA_active,NA)
meanLA.m<-mean(cwt1$mLA, na.rm=TRUE)
cwt1$LA.m<-cwt1$mLA/meanLA.m

cwt1$fLA.<-ifelse(cwt1$sex=="female",cwt1$LA_active,NA)
meanLA.f<-mean(cwt1$fLA, na.rm=TRUE)
cwt1$LA.f<-cwt1$fLA/meanLA.f

cwt1$mWS<-ifelse(cwt1$sex=="male",cwt1$WS,NA)
meanWs.m<-mean(cwt1$mWS,na.rm=TRUE)
cwt1$meanWs.m<-cwt1$mWS/meanWs.m

cwt1$fWS<-ifelse(cwt1$sex=="female",cwt1$WS,NA)
meanWS.f<-mean(cwt1$fWS,na.rm=TRUE)
cwt1$meanWs.f<-cwt1$fWS/meanWS.f

cwt1$mDR<-ifelse(cwt1$sex=="male",cwt1$DR ,NA)
meanDR.m<-mean(cwt1$mDR,na.rm=TRUE)
cwt1$meanDr.m<-cwt1$mDR/meanDR.m

cwt1$fDR<-ifelse(cwt1$sex=="female",cwt1$DR,NA)
meanDR.f<-mean(cwt1$fDR,na.rm=TRUE)
cwt1$meanDr.f<-cwt1$fDR/meanDR.f

cwt1$WS<-ifelse(cwt1$sex=="male",cwt1$meanWs.m,cwt1$meanWs.f)
cwt1$DR<-ifelse(cwt1$sex=="male",cwt1$meanDr.m,cwt1$meanDr.f)
cwt1$LA_active<-ifelse(cwt1$sex=="male",cwt1$LA.m,cwt1$LA.f)

cwt1<-cwt1[,1:7]

#Set up trait numbers
cwt1$TRAITNO<-NA
#Select wing size data
cwt1.WS<-cwt1[c(8,1,4:7)]
#Code traits, 1 = male wing size, 4 = female wing size
cwt1.WS$TRAITNO<-ifelse(cwt1.WS$sex=="male",1,4)
#Rename trait column
names(cwt1.WS)[2] = "trait"
#Select DR data
cwt1.DR<-cwt1[c(8,2,4:7)]
#Code traits, 2 = male DR, 5 = female DR
cwt1.DR$TRAITNO<-ifelse(cwt1.DR$sex=="male",2,5)
#Rename trait column
names(cwt1.DR)[2] = "trait"
#Select LA data
cwt1.LA<-cwt1[c(8,3:7)]
#Code traits, 3 = male LA, 6 = female LA
cwt1.LA$TRAITNO<-ifelse(cwt1.LA$sex=="male",3,6)
#Rename trait column
names(cwt1.LA)[2] = "trait"
#Put together into single dataset
cwt1.dat<-rbind(cwt1.WS,cwt1.DR,cwt1.LA)
#Code trait by sex, Males = 1, Females = 2
cwt1.dat$sex<-ifelse(cwt1.dat$sex=="male",1,2)
#Make round.sec into integers by taking factor codes
cwt1.dat$round.sec <- as.numeric(cwt1.dat$round.sec)
#drop unused levels
cwt1.dat<-droplevels(cwt1.dat) 
str(cwt1.dat)
#round: 4 levels
length(unique(cwt1.dat$round.sec))
#round.sec: 28 levels
length(unique(cwt1.dat$animal))
#animal: 914 levels

#Change ids into numeric values
cwt1.dat$animal<-gsub("M","1",cwt1.dat$animal)
cwt1.dat$animal<-gsub("F","2",cwt1.dat$animal)
cwt1.dat$animal<-gsub("m","3",cwt1.dat$animal)
cwt1.dat$animal<-gsub("f","4",cwt1.dat$animal)
cwt1.dat$animal<-gsub("-","",cwt1.dat$animal)
#Sort by animal ID and then trait number
cwt1.dat <- cwt1.dat[order(cwt1.dat$animal, cwt1.dat$TRAITNO),]

#Change name of round.sec to NSEC
names(cwt1.dat)[6]<-"NSEC"

#write data file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1")
write.table(cwt1.dat, file="cwt1_dat.dat", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Create pedigree.
#M=1 (sire), F=2 (mother), m=3 (male offspring), f=4 (female offspring)
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

#Import data
ped.cwt1 <- read.csv("rep1cwt.p.csv")
head(ped.cwt1)
str(ped.cwt1)
#Change ids into numeric values
ped.cwt1$id<-gsub("M","1",ped.cwt1$id)
ped.cwt1$id<-gsub("F","2",ped.cwt1$id)
ped.cwt1$id<-gsub("m","3",ped.cwt1$id)
ped.cwt1$id<-gsub("f","4",ped.cwt1$id)
ped.cwt1$id<-gsub("-","",ped.cwt1$id)
ped.cwt1$FATHER<-gsub("M","1",ped.cwt1$FATHER)
ped.cwt1$FATHER[is.na(ped.cwt1$FATHER)] <- 0
ped.cwt1$MOTHER<-gsub("F","2",ped.cwt1$MOTHER)
ped.cwt1$MOTHER[is.na(ped.cwt1$MOTHER)] <- 0

#write pedigree file
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\Wombat analysis\\cwt1")
write.table(ped.cwt1, file="cwt1_ped.ped", row.names = FALSE, quote = FALSE, col.names = FALSE)

#Get starting values for par file from MCMCglmm results
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
rep1cwt <- readRDS("rep1cwt")
library(MCMCglmm)
vcvmatrix<-posterior.mode(rep1cwt$VCV)
#VAR ANIMAL
vcvmatrix[c(1:6,8:12,15:18,22:24,29:30,36)]
#traitmeanWs.m:traitmeanWs.m.animal         traitmeanDr.m:traitmeanWs.m.animal 
#0.0010135619                               0.0017243495 
#traitmeanActive.m:traitmeanWs.m.animal         traitmeanWs.f:traitmeanWs.m.animal 
#0.0062189367                               0.0007081119 
#traitmeanDr.f:traitmeanWs.m.animal     traitmeanActive.f:traitmeanWs.m.animal 
#0.0014235847                               0.0070785664 
#traitmeanDr.m:traitmeanDr.m.animal     traitmeanActive.m:traitmeanDr.m.animal 
#0.0780701071                              -0.0636919679 
#traitmeanWs.f:traitmeanDr.m.animal         traitmeanDr.f:traitmeanDr.m.animal 
#0.0020567686                               0.0766817182 
#traitmeanActive.f:traitmeanDr.m.animal traitmeanActive.m:traitmeanActive.m.animal 
#-0.0119690798                               0.4604724963 
#traitmeanWs.f:traitmeanActive.m.animal     traitmeanDr.f:traitmeanActive.m.animal 
#0.0077738434                              -0.0481697554 
#traitmeanActive.f:traitmeanActive.m.animal         traitmeanWs.f:traitmeanWs.f.animal 
#0.3800519872                               0.0008778647 
#traitmeanDr.f:traitmeanWs.f.animal     traitmeanActive.f:traitmeanWs.f.animal 
#0.0015773463                               0.0077539090 
#traitmeanDr.f:traitmeanDr.f.animal     traitmeanActive.f:traitmeanDr.f.animal 
#0.0856468329                               0.0407111088 
#traitmeanActive.f:traitmeanActive.f.animal 
#0.5826502884

#VAR NSEC
vcvmatrix[c(37:42,44:48,51:54,58:60,65:66,72)]
#traitmeanWs.m:traitmeanWs.m.new.sec         traitmeanDr.m:traitmeanWs.m.new.sec 
#4.382112e-04                               -1.529295e-04 
#traitmeanActive.m:traitmeanWs.m.new.sec         traitmeanWs.f:traitmeanWs.m.new.sec 
#-2.207460e-04                                1.328045e-04 
#traitmeanDr.f:traitmeanWs.m.new.sec     traitmeanActive.f:traitmeanWs.m.new.sec 
#-2.223878e-05                               -1.057756e-04 
#traitmeanDr.m:traitmeanDr.m.new.sec     traitmeanActive.m:traitmeanDr.m.new.sec 
#3.430672e-02                                1.126425e-02 
#traitmeanWs.f:traitmeanDr.m.new.sec         traitmeanDr.f:traitmeanDr.m.new.sec 
#-9.949411e-05                                2.140604e-02 
#traitmeanActive.f:traitmeanDr.m.new.sec traitmeanActive.m:traitmeanActive.m.new.sec 
#-2.584849e-03                                1.003059e-02 
#traitmeanWs.f:traitmeanActive.m.new.sec     traitmeanDr.f:traitmeanActive.m.new.sec 
#-2.897188e-04                                2.939581e-04 
#traitmeanActive.f:traitmeanActive.m.new.sec         traitmeanWs.f:traitmeanWs.f.new.sec 
#2.810117e-03                                4.696593e-04 
#traitmeanDr.f:traitmeanWs.f.new.sec     traitmeanActive.f:traitmeanWs.f.new.sec 
#-1.311094e-04                               -2.665869e-05 
#traitmeanDr.f:traitmeanDr.f.new.sec     traitmeanActive.f:traitmeanDr.f.new.sec 
#1.922956e-02                                2.481056e-03 
#traitmeanActive.f:traitmeanActive.f.new.sec 
#1.019808e-02 

#VAR RESIDUAL
vcvmatrix[c(73:78,80:84,87:90,94:96,101:102,108)]
#traitmeanWs.m:traitmeanWs.m.units         traitmeanDr.m:traitmeanWs.m.units 
#4.011628e-04                             -9.543139e-04 
#traitmeanActive.m:traitmeanWs.m.units         traitmeanWs.f:traitmeanWs.m.units 
#-2.186521e-03                              1.594795e-05 
#traitmeanDr.f:traitmeanWs.m.units     traitmeanActive.f:traitmeanWs.m.units 
#8.540481e-04                              1.482159e-03 
#traitmeanDr.m:traitmeanDr.m.units     traitmeanActive.m:traitmeanDr.m.units 
#7.403146e-02                              6.884633e-02 
#traitmeanWs.f:traitmeanDr.m.units         traitmeanDr.f:traitmeanDr.m.units 
#-7.676350e-04                              7.719561e-02 
#traitmeanActive.f:traitmeanDr.m.units traitmeanActive.m:traitmeanActive.m.units 
#1.924053e-01                              9.682296e-01 
#traitmeanWs.f:traitmeanActive.m.units     traitmeanDr.f:traitmeanActive.m.units 
#-1.238393e-03                              9.282551e-02 
#traitmeanActive.f:traitmeanActive.m.units         traitmeanWs.f:traitmeanWs.f.units 
#1.153693e+00                              4.266777e-04 
#traitmeanDr.f:traitmeanWs.f.units     traitmeanActive.f:traitmeanWs.f.units 
#8.872743e-04                             -3.075394e-03 
#traitmeanDr.f:traitmeanDr.f.units     traitmeanActive.f:traitmeanDr.f.units 
#8.432051e-02                             -1.090081e-02 
#traitmeanActive.f:traitmeanActive.f.units 
#1.653955e+00 