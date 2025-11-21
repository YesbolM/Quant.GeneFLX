#Covariance tensor matrix comparisons

#Import data
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")

model.flx1<-readRDS("rep1flx")
VC_rep1flx <- model.flx1$VCV
model.cfm1<-readRDS("rep1cfm")
VC_rep1cfm <- model.cfm1$VCV
model.cwt1<-readRDS("rep1cwt")
VC_rep1cwt <- model.cwt1$VCV
model.flx2<-readRDS("rep2flx")
VC_rep2flx <- model.flx2$VCV
model.cfm2<-readRDS("rep2cfm")
VC_rep2cfm <- model.cfm2$VCV
model.cwt2<-readRDS("rep2cwt")
VC_rep2cwt <- model.cwt2$VCV

#Set random number seed to make analysis repeatable
set.seed(2112)

#number of MCMC samples
MCMCsamp <- 1000 
#number of traits 
n <- 6 
#number of matrices to compare
m <- 6
#number of random effects specified in the model. In our analyses these were animal, section and residual effects.
r <- 3
#trait names
traitnames <- c("Ws.m","Dr.m","Active.m","Ws.f","Dr.f","Active.f") 
#matrix labels 
Gnames <- c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2")

#empty array 
MCMCarray <- array(NA,c(MCMCsamp,(n^2)*r,m)) 
#G1 stored as the 1st element of dim[3] 
MCMCarray[,,1] <- as.matrix(VC_rep1flx)
#G2 stored as the 2nd element of dim[3]
MCMCarray[,,2] <- as.matrix(VC_rep1cfm)
#G3 stored as the 3rd element of dim[3]
MCMCarray[,,3] <- as.matrix(VC_rep1cwt)
#G4 stored as the 4th element of dim[3]
MCMCarray[,,4] <- as.matrix(VC_rep2flx) 
#G5 stored as the 5th element of dim[3]
MCMCarray[,,5] <- as.matrix(VC_rep2cfm) 
#G6 stored as the 6th element of dim[3]
MCMCarray[,,6] <- as.matrix(VC_rep2cwt) 

#Reshaping the array and standardizing G
Garray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)
Parray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames,Gnames)
for (i in 1:m){
  for (j in 1:MCMCsamp){
    G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
    CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    R <- matrix(MCMCarray[j,(((n^2)*2)+1):((n^2)*3),i],ncol= n)
    Garray[,,i,j] <- G
    Parray[,,i,j] <- G + CE + R
  }
}

#Define standardization function c.f. Hansen and Houle (2008)
#Should be used if traits are not on the same scale. The HHGarray object will contain the standardized G array. 
#I assume that HHGarray should be used anywhere below where the Garray is called.
inv.rootP <- function (P){
  rootP <- matrix(0,n, n)  
  for (i in 1:n){
    val <- eigen(P)$values
    vec <- eigen(P)$vectors
    rootP <- rootP + (vec[,i] %*% t(vec[,i]))*sqrt(val[i])
  }
  solve(rootP)
}

HHGarray <- array(NA,c(n,n,m,MCMCsamp))
for (k in 1:MCMCsamp){
  for (j in 1:m){
    P <- inv.rootP(Parray[,,j,k])
    HHGarray[,,j,k] <- P %*% Garray[,,j,k] %*% P
  }
}

#Generating randomised G matrices for hypothesis tests. 
#NB - Morissey et al. 2019 suggest a more complicated randomization procedure where phenotypes for each individual are simulated from the output of the MCMCglmm model and then reanalysed in the same way as the observed data. This is not practical to implement for these datasets due to runtime issues, but we apply the Walters et al 2018 randomization procedure below.
#Import pedigree data
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")
Ped.flx1 <- as.data.frame(read.csv(file="rep1FLX.p.csv",header =T))
Ped.cfm1 <- as.data.frame(read.csv(file="rep1CFM.p.csv",header =T))
Ped.cwt1 <- as.data.frame(read.csv(file="rep1CWT.p.csv",header =T))
Ped.flx2 <- as.data.frame(read.csv(file="rep2FLX.p.csv",header =T))
Ped.cfm2 <- as.data.frame(read.csv(file="rep2CFM.p.csv",header =T))
Ped.cwt2 <- as.data.frame(read.csv(file="rep2CWT.p.csv",header =T))

library(MCMCglmm)

#Simulate a set of G that have been sampled from the same population, so the only dissimilarity among them is random sampling error
rand.Garray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(rand.Garray) <- list(traitnames,traitnames,Gnames)
for (i in 1:MCMCsamp){
  flx1.bv<-rbv(Ped.flx1,HHGarray[,,1,i])
  cfm1.bv<-rbv(Ped.cfm1,HHGarray[,,2,i])
  cwt1.bv<-rbv(Ped.cwt1,HHGarray[,,3,i])
  flx2.bv<-rbv(Ped.flx2,HHGarray[,,4,i])
  cfm2.bv<-rbv(Ped.cfm2,HHGarray[,,5,i])
  cwt2.bv<-rbv(Ped.cwt2,HHGarray[,,6,i])
  a.pop <- cumsum(c(dim(Ped.flx1)[1],dim(Ped.cfm1)[1],dim(Ped.cwt1)[1],dim(Ped.flx2)[1],dim(Ped.cfm2)[1],dim(Ped.cwt2)[1]))
  pop.bv <- rbind(flx1.bv,cfm1.bv,cwt1.bv,flx2.bv,cfm2.bv,cwt2.bv)
  rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
  rand.Garray[,,1,i] <- cov(rand.pop.bv[1:a.pop[1],])
  rand.Garray[,,2,i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2],])
  rand.Garray[,,3,i] <- cov(rand.pop.bv[(a.pop[2] + 1):a.pop[3],])
  rand.Garray[,,4,i] <- cov(rand.pop.bv[(a.pop[3] + 1):a.pop[4],])
  rand.Garray[,,5,i] <- cov(rand.pop.bv[(a.pop[4] + 1):a.pop[5],])
  rand.Garray[,,6,i] <- cov(rand.pop.bv[(a.pop[5] + 1):a.pop[6],])
}

#required packages
library(gdata)
library(matrixcalc)
library(MCMCglmm)

#Tensor method c.f. Hine (2009)
#Create function to carry out analysis
#START
covtensor <- function(Gs){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  neigten <- n*(n+1)/2 
  #Number of eigentensors
  MCMC.S <- array(NA,c(neigten, neigten, MCMCsamp))
  dimnames(MCMC.S) <- list(paste("e", 1:neigten, sep=""), paste("e", 1:neigten, sep=""))
  for (k in 1:MCMCsamp){
    MCMCG <- Gs[,,,k] 
    MCMCvarmat <- t(apply(MCMCG, 3, diag)) 
    #find the variances of the kth G and store them 
    MCMCcovmat <- t(apply(MCMCG, 3, lowerTriangle)) 
    #find the covariances of the kth G and store them
    MCMC.S[1:n,1:n, k] <- cov(MCMCvarmat, MCMCvarmat) 
    #fill the upper left quadrant of the kth S
    MCMC.S[(n+1):neigten,(n+1):neigten, k] <- 2*cov(MCMCcovmat, MCMCcovmat)
    #fill the lower right quadrant of the kth S
    MCMC.S[1:n,(n+1):neigten, k] <- sqrt(2)*cov(MCMCvarmat, MCMCcovmat)
    #fill the upper right quadrant of the kth S
    MCMC.S[(n+1):neigten,1:n, k] <- sqrt(2)*cov(MCMCcovmat, MCMCvarmat)
    #fill the lower left quadrant of the kthS
  }  
  av.S <- apply(MCMC.S, 1:2, mean)
  #posterior mean S
  av.S.val <- eigen(av.S)$values
  #eigenvalues of posterior mean S 
  av.S.vec <- eigen(av.S)$vectors
  #eigenvalues of posterior mean S
  eTmat <- array(NA, c(n, n, neigten))
  dimnames(eTmat) <- list(traitnames, traitnames, paste("E", 1:neigten, sep=""))  
  for (i in 1:neigten){
    emat <- matrix(0, n, n) 
    lowerTriangle(emat) <- 1/sqrt(2)*av.S.vec[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- av.S.vec[1:n,i]
    eTmat[,,i] <- emat 
  }
  #construct the second-order eigentensors of posterior mean S
  eT.eigen <- array(NA, c(n+1, n, neigten))
  for (i in 1:neigten){
    eT.eigen[1,,i] <- t(eigen(eTmat[,,i])$values) 
    #Eigenvalues of the ith eigentensor
    eT.eigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors 
    #Eigenvectors of the ith eigentensor
    eT.eigen[,,i] <- eT.eigen[,order(abs(eT.eigen[1,,i]), decreasing = T), i]
  }
  MCMC.S.val <- matrix(NA, MCMCsamp, neigten)
  colnames(MCMC.S.val) <- paste("E", 1:neigten, sep="")
  for (i in 1:MCMCsamp){
    for(j in 1:neigten){
      MCMC.S.val[i,j] <- t(av.S.vec[,j]) %*% MCMC.S[,,i] %*% av.S.vec[,j]
    }
  }
  #posterior distribution of the genetic variance for the eigenvectors of posterior mean S
  av.G.coord <- array(NA, c(m, neigten, 1))
  dimnames(av.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    av.G.coord[,i,] <- apply((apply(Gs, 1:3, mean)) , 3, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the jth avG for the eigentensors of posterior mean S
  MCMC.G.coord <- array(NA, c(m, neigten, MCMCsamp))
  dimnames(MCMC.G.coord) <- list(Gnames, paste("E", 1:neigten, sep=""))
  for (i in 1:neigten){
    MCMC.G.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
  }
  #Coordinates of the kth MCMC sample of the jth G for the eigentensors of posterior mean S
  tensor.summary <- data.frame(rep(av.S.val,each=n), t(data.frame(eT.eigen)))
  colnames(tensor.summary) <- c("S.eigval", "eT.val", traitnames)
  rownames(tensor.summary)<- paste(paste("e", rep(1:neigten, each=n), sep=""), rep(1:n,neigten), sep=".")
  list(tensor.summary = tensor.summary, av.S = av.S, eTmat = eTmat, av.G.coord = av.G.coord, MCMC.S = MCMC.S, MCMC.S.val = MCMC.S.val, MCMC.G.coord = MCMC.G.coord)
}
#END

#Apply covtensor to the observed G array
MCMC.covtensor <- covtensor(Garray)

#calculate the maximum number of nonzero eigentensors
nnonzero <- min(n*(n+1)/2,m-1)
nnonzero
#In this case the result is 5

#Apply the covtensor function to the randomised G array for significance testing
MCMC.covtensor.rand <- covtensor(rand.Garray)

#Compare observed and randomized posterior means and 95% HPD intervals of the eigenvalues for the nonzero eigentensors
HPD.eT.val <- cbind(mean.obs=posterior.mode(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero])), HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), mean.rand=posterior.mode(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero])), HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]), prob=0.95))
round(HPD.eT.val, 3)
#Observed on the left half of the table, randomized on the right. We want to find the ones that don't overlap.
#It seems like there could be some evidence of a difference in E1 and E2.

#Plot equivalent to figure 5 in the tutorial:
bp<-barplot(t(as.matrix(HPD.eT.val[,c(1,4)])), legend=c("Observed","Randomized"), xlab="Alpha", ylim=c(0, 0.6), beside=TRUE)
arrows(bp, t(as.matrix(HPD.eT.val[,c(2,5)])), bp, t(as.matrix(HPD.eT.val[,c(3,6)])), code=3, angle=90)

#To examine trait combinations for E1 and E2 we can examine the first 12 (i.e. n*2) rows of the $tensor.summary: 
round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)
#Variation explained by each eigenvector of the eigentensors:
x<-MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]]
E1<-abs(x[1:n,1])/sum(abs(x[1:n,1]))
E1
#e11 = 87% of variation explained
E2<-abs(x[(n+1):(n*2),1])/sum(abs(x[(n+1):(n*2),1]))
E2
#e21 = 57% of variation explained
#e22 = 38% of variation explained
#These three eigenvectors seem useful to work with.

#e11 loads heavily on male and female activity, in the same direction (positive)
#e21 loads heavily on male activity (negative)
#e22 loads heavily on female activity (positive) and to some extent female dessication resistance (negative)

#The contribution of each population to the eigenvalues of the tensor can be visualized by plotting the coordinates of each G in the space of eigentensors E1 and E2. Plot equivalent to figure 6 in the tutorial:
str(MCMC.covtensor$MCMC.G.coord)
#dim 1 = populations
#dim 2 = Eigentensors (E1-E...)
#dim 3 = MCMC values

#For E1
flx1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,1,]))
flx1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,1,]))
cfm1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,1,]))
cfm1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,1,]))
cwt1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,1,]))
cwt1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,1,]))
flx2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,1,]))
flx2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,1,]))
cfm2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,1,]))
cfm2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,1,]))
cwt2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,1,]))
cwt2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,1,]))
#Plot
bp <- barplot(c(flx1.mean.E1, cfm1.mean.E1, cwt1.mean.E1, flx2.mean.E1, cfm2.mean.E1, cwt2.mean.E1), ylab="E1", xlab="Treatment", names=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), ylim=c(-2,0.25))
arrows(bp, c(flx1.CI.E1[1], cfm1.CI.E1[1], cwt1.CI.E1[1], flx2.CI.E1[1], cfm2.CI.E1[1], cwt2.CI.E1[1]), bp, c(flx1.CI.E1[2], cfm1.CI.E1[2], cwt1.CI.E1[2], flx2.CI.E1[2], cfm2.CI.E1[2], cwt2.CI.E1[2]), code=3, angle=90)
#The main difference along E1 is between CWT and FLX in rep1, but between FLX and CFM in rep2. This dimension seems to represent variation in the amount of additive genetic variance in male and female locomotory activity - most in CWT and least in FLX for both sexes in rep 1. (I don't think it is capturing changes in the between sex genetic correlation for this trait because CFM is the odd one out in the point estimates.)

#For E2
flx1.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,2,]))
flx1.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,2,]))
cfm1.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,2,]))
cfm1.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,2,]))
cwt1.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,2,]))
cwt1.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,2,]))
flx2.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,2,]))
flx2.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,2,]))
cfm2.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,2,]))
cfm2.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,2,]))
cwt2.mean.E2 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,2,]))
cwt2.CI.E2 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,2,]))
#Plot
bp <- barplot(c(flx1.mean.E2, cfm1.mean.E2, cwt1.mean.E2, flx2.mean.E2, cfm2.mean.E2, cwt2.mean.E2), ylab="E2", xlab="Treatment", names=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), ylim=c(-0.8,0.5))
arrows(bp, c(flx1.CI.E2[1], cfm1.CI.E2[1], cwt1.CI.E2[1], flx2.CI.E2[1], cfm2.CI.E2[1], cwt2.CI.E2[1]), bp, c(flx1.CI.E2[2], cfm1.CI.E2[2], cwt1.CI.E2[2], flx2.CI.E2[2], cfm2.CI.E2[2], cwt2.CI.E2[2]), code=3, angle=90)
#The main difference along E2 is between FLX and CFM in both replicates. It is mainly explained by male activity in e21, and female activity (and to some extent male desiccation resistance) in e22. This dimension seems to represent changes in the intersexual genetic correlation for locomotory activity, as well as maybe the cross-sex genetic correlation for locomotory acitivity and desiccation resistance.

#Examine the contribution of specific trait combinations to coordinated changes among G by projecting the eigenvectors of eigentensors on the observed G array. 
#e11 vector
e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
#e21 vector
e21 <- c(as.numeric(MCMC.covtensor$tensor.summary[(n+1),3:dim(MCMC.covtensor$tensor.summary)[2]]))
#e22 vector
e22 <- c(as.numeric(MCMC.covtensor$tensor.summary[(n+2),3:dim(MCMC.covtensor$tensor.summary)[2]]))

#Function to do projection
proj<- function(G, b) t(b) %*% G %*% (b)

#genetic variance along E1 e1, 1 for each MCMC sample of each replicate line
e11.proj <- apply(Garray, 3:4, proj, b = e11)
#genetic variance along E2 e1, 1 for each MCMC sample of each replicate line
e21.proj <- apply(Garray, 3:4, proj, b = e21)
#genetic variance along E2 e2, 1 for each MCMC sample of each replicate line
e22.proj <- apply(Garray, 3:4, proj, b = e22)
#These projections can be used to visualize the posterior mean and 95% HPD interval for the genetic variance along the direction of e11 and e21 for each replicate line (Figure 7 in the tutorial).
str(e11.proj)

#For e11
flx1.mean.e11 <- mean(e11.proj[1,])
flx1.CI.e11 <- HPDinterval(as.mcmc(e11.proj[1,]))
cfm1.mean.e11 <- mean(e11.proj[2,])
cfm1.CI.e11 <- HPDinterval(as.mcmc(e11.proj[2,]))
cwt1.mean.e11 <- mean(e11.proj[3,])
cwt1.CI.e11 <- HPDinterval(as.mcmc(e11.proj[3,]))
flx2.mean.e11 <- mean(e11.proj[4,])
flx2.CI.e11 <- HPDinterval(as.mcmc(e11.proj[4,]))
cfm2.mean.e11 <- mean(e11.proj[5,])
cfm2.CI.e11 <- HPDinterval(as.mcmc(e11.proj[5,]))
cwt2.mean.e11 <- mean(e11.proj[6,])
cwt2.CI.e11 <- HPDinterval(as.mcmc(e11.proj[6,]))
#Plot
bp <- barplot(c(flx1.mean.e11, cfm1.mean.e11, cwt1.mean.e11, flx2.mean.e11, cfm2.mean.e11, cwt2.mean.e11), ylab="Lambda", xlab="Treatment", names=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), ylim=c(0,2))
arrows(bp, c(flx1.CI.e11[1], cfm1.CI.e11[1], cwt1.CI.e11[1], flx2.CI.e11[1], cfm2.CI.e11[1], cwt2.CI.e11[1]), bp, c(flx1.CI.e11[2], cfm1.CI.e11[2], cwt1.CI.e11[2], flx2.CI.e11[2], cfm2.CI.e11[2], cwt2.CI.e11[2]), code=3, angle=90)
#The main difference along e1 is between CWT and FLX in rep 1, and between CFM and FLX in rep 2. Again, this dimension seems to represent variation in the amount of additive genetic variance in male and female locomotory activity - most in CWT and least in FLX for both sexes in rep 1. 

#For e21
flx1.mean.e21 <- mean(e21.proj[1,])
flx1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[1,]))
cfm1.mean.e21 <- mean(e21.proj[2,])
cfm1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[2,]))
cwt1.mean.e21 <- mean(e21.proj[3,])
cwt1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[3,]))
flx2.mean.e21 <- mean(e21.proj[4,])
flx2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[4,]))
cfm2.mean.e21 <- mean(e21.proj[5,])
cfm2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[5,]))
cwt2.mean.e21 <- mean(e21.proj[6,])
cwt2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[6,]))
#Plot
bp <- barplot(c(flx1.mean.e21, cfm1.mean.e21, cwt1.mean.e21, flx2.mean.e21, cfm2.mean.e21, cwt2.mean.e21), ylab="Lambda", xlab="Treatment", names=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), ylim=c(0,1))
arrows(bp, c(flx1.CI.e21[1], cfm1.CI.e21[1], cwt1.CI.e21[1], flx2.CI.e21[1], cfm2.CI.e21[1], cwt2.CI.e21[1]), bp, c(flx1.CI.e21[2], cfm1.CI.e21[2], cwt1.CI.e21[2], flx2.CI.e21[2], cfm2.CI.e21[2], cwt2.CI.e21[2]), code=3, angle=90)
#The main difference along E2 is between FLX and CFM for both replicates. This dimension mainly captures variance in male activity in e21. Values are also lower overall for rep 2.

#For e22
flx1.mean.e22 <- mean(e22.proj[1,])
flx1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[1,]))
cfm1.mean.e22 <- mean(e22.proj[2,])
cfm1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[2,]))
cwt1.mean.e22 <- mean(e22.proj[3,])
cwt1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[3,]))
flx2.mean.e22 <- mean(e22.proj[4,])
flx2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[4,]))
cfm2.mean.e22 <- mean(e22.proj[5,])
cfm2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[5,]))
cwt2.mean.e22 <- mean(e22.proj[6,])
cwt2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[6,]))
#Plot
bp <- barplot(c(flx1.mean.e22, cfm1.mean.e22, cwt1.mean.e22, flx2.mean.e22, cfm2.mean.e22, cwt2.mean.e22), ylab="Lambda", xlab="Treatment", names=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), ylim=c(0,1.3))
arrows(bp, c(flx1.CI.e22[1], cfm1.CI.e22[1], cwt1.CI.e22[1], flx2.CI.e22[1], cfm2.CI.e22[1], cwt2.CI.e22[1]), bp, c(flx1.CI.e22[2], cfm1.CI.e22[2], cwt1.CI.e22[2], flx2.CI.e22[2], cfm2.CI.e22[2], cwt2.CI.e22[2]), code=3, angle=90)
#The main difference along e22 is between CWT and the other treatments in rep 1, and CFM and the other treatmentst in rep 2. It is mainly explained by differences in female activity, partly in combination with male dessication resistance. This dimension seems to represent changes in the additive genetic variance in female locomotory activity, as well as the cross-sex genetic correlation between desiccation resistance and locomotory activity.

#Since e21 and e22 explain similar amounts of variation, I think we really need to try to interpret them in combination.
#Plot of the projections of e21 and e22 against each other:
plot(c(flx1.mean.e21, cfm1.mean.e21, cwt1.mean.e21, flx2.mean.e21, cfm2.mean.e21, cwt2.mean.e21),c(flx1.mean.e22, cfm1.mean.e22, cwt1.mean.e22, flx2.mean.e22, cfm2.mean.e22, cwt2.mean.e22), xlab="e21", ylab="e22", xlim=c(0,1), ylim=c(0,1.2), pch=c(16,17,18), col=rep(c("black","blue"), each=3), cex=2)
arrows(flx1.CI.e21[1], flx1.mean.e22, flx1.CI.e21[2], flx1.mean.e22, code=3, angle=90) #FLX1 in e21
arrows(flx1.mean.e21, flx1.CI.e22[1],  flx1.mean.e21, flx1.CI.e22[2], code=3, angle=90) #FLX1 in e22
arrows(cfm1.CI.e21[1], cfm1.mean.e22, cfm1.CI.e21[2], cfm1.mean.e22, code=3, angle=90) #CFM1 in e21
arrows(cfm1.mean.e21, cfm1.CI.e22[1],  cfm1.mean.e21, cfm1.CI.e22[2], code=3, angle=90) #CFM1 in e22
arrows(cwt1.CI.e21[1], cwt1.mean.e22, cwt1.CI.e21[2], cwt1.mean.e22, code=3, angle=90) #CWT1 in e21
arrows(cwt1.mean.e21, cwt1.CI.e22[1],  cwt1.mean.e21, cwt1.CI.e22[2], code=3, angle=90) #CWT1 in e22
arrows(flx2.CI.e21[1], flx2.mean.e22, flx2.CI.e21[2], flx2.mean.e22, code=3, angle=90, col="blue") #FLX2 in e21
arrows(flx2.mean.e21, flx2.CI.e22[1],  flx2.mean.e21, flx2.CI.e22[2], code=3, angle=90, col="blue") #FLX2 in e22
arrows(cfm2.CI.e21[1], cfm2.mean.e22, cfm2.CI.e21[2], cfm2.mean.e22, code=3, angle=90, col="blue") #CFM2 in e21
arrows(cfm2.mean.e21, cfm2.CI.e22[1],  cfm2.mean.e21, cfm2.CI.e22[2], code=3, angle=90, col="blue") #CFM2 in e22
arrows(cwt2.CI.e21[1], cwt2.mean.e22, cwt2.CI.e21[2], cwt2.mean.e22, code=3, angle=90, col="blue") #CWT2 in e21
arrows(cwt2.mean.e21, cwt2.CI.e22[1],  cwt2.mean.e21, cwt2.CI.e22[2], code=3, angle=90, col="blue") #CWT2 in e22
legend("topright", legend=c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2"), pch=c(16,17,18), col=rep(c("black","blue"), each=3))

###############################################################################

#Applying Walters et al 2018 randomization procedure
#Randomize phenotypic data
#Run animal model
#Plot difference between true mean and posterior mean for all samples up to that point, check where 95% confidence intervals seem to flatten out
#Use that value to determine number of iterations necessary for each randomization.
#Do covariance tensor analysis on all 1000 randomizations.
#Extract eigenvalues for all non-zero eigentensors (E1-E6) from the randomized analyses.
#Use the distribution of the eigenvalues as the null expectation.

#Use sample with replace=FALSE to randomize phenotypic data
#According to Morrissey et al 2019, Walters et al 2018 randomized phenotypic data relative to the pedigree, i.e. measurements from the same individual are kept together but they are randomized among families. We decided to do this and randomize within rounds and sexes in order to keep the differences between the rounds and sexes intact.

#Do one full randomization per population in order to check that results are similar.

#FLX1 full randomization
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

data <- read.csv("data.csv")
str(data)
data<-as.data.frame(data) # 6858 obs. of  13 variables
names(data)[1]<-"animal"
for(x in c(1:5,12:13))data[,x]<-as.factor(data[,x])
data$mLA.active<-ifelse(data$sex=="male",data$LA_active,NA)
meanActive.m<-mean(data$mLA.activ, na.rm=TRUE)
data$meanActive.m<-data$mLA.active/meanActive.m #mean standardization

data$mLA.passive<-ifelse(data$sex=="male",data$LA_passive,NA)
data$fLA.active<-ifelse(data$sex=="female",data$LA_active,NA)
meanActive.f<-mean(data$fLA.activ, na.rm=TRUE)
data$meanActive.f<-data$fLA.active/meanActive.f

data$fLA.passive<-ifelse(data$sex=="female",data$LA_passive,NA)
data$mWS<-ifelse(data$sex=="male",data$WS,NA)
meanWs.m<-mean(data$mWS,na.rm=TRUE)
data$meanWs.m<-data$mWS/meanWs.m

data$fWS<-ifelse(data$sex=="female",data$WS,NA)
meanWS.f<-mean(data$fWS,na.rm=TRUE)
data$meanWs.f<-data$fWS/meanWS.f

data$mDR<-ifelse(data$sex=="male",data$DR ,NA)
meanDR.m<-mean(data$mDR,na.rm=TRUE)
data$meanDr.m<-data$mDR/meanDR.m

data$fDR<-ifelse(data$sex=="female",data$DR,NA)
meanDR.f<-mean(data$fDR,na.rm=TRUE)
data$meanDr.f<-data$fDR/meanDR.f

rep1<-data[data$replicat=="rep1",]
FLX1<-rep1[rep1$treatment=="FLX",]
CFM1<-rep1[rep1$treatment=="CFM",]
CWT1<-rep1[rep1$treatment=="CWT",]

rep2<-data[data$replicat=="rep2",]
FLX2<-rep2[rep2$treatment=="FLX",]
CFM2<-rep2[rep2$treatment=="CFM",]
CWT2<-rep2[rep2$treatment=="CWT",]

rep1FLX_p <- read.csv("rep1FLX.p.csv")
ped1<-rep1FLX_p
head(ped1)
str(ped1)
ped1<-as.data.frame(ped1)

rep1CFM_p <- read.csv("rep1CFM.p.csv")
ped2<-rep1CFM_p
head(ped2)
str(ped2)
ped2<-as.data.frame(ped2)

rep1CWT_p <- read.csv("rep1CWT.p.csv")
ped3<-rep1CWT_p
head(ped3)
str(ped3)
ped3<-as.data.frame(ped3)

rep2FLX_p <- read.csv("rep2FLX.p.csv")
ped4<-rep2FLX_p
head(ped4)
str(ped4)
ped1<-as.data.frame(ped4)

rep2CFM_p <- read.csv("rep2CFM.p.csv")
ped5<-rep2CFM_p
head(ped5)
str(ped5)
ped5<-as.data.frame(ped5)

rep2CWT_p <- read.csv("rep2CWT.p.csv")
ped6<-rep2CWT_p
head(ped6)
str(ped6)
ped6<-as.data.frame(ped6)

library(MCMCglmm)

# Rep1 FLX
#Randomization
#4 rounds, 2 sexes = 8 subgroups. The randomization below isn't very elegant, but it gets the job done.
#Round 1 females
r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
FLX1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(FLX1.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1flxl.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                   random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=FLX1.rand,prior=p2,verbose=TRUE,
                   #family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1flxl.rand,file = "rep1flxl.rand")

# Rep1 CFM
#Randomization
#Round 1 females
r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
CFM1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(CFM1.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1cfm.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                        random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped2,data=CFM1.rand,prior=p2,verbose=TRUE,
                        family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cfm.rand,file = "rep1cfm.rand")

# Rep1 CWT
#Randomization
#Round 1 females
r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
CWT1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(CWT1.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep1cwt.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                       random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped3,data=CWT1.rand,prior=p2,verbose=TRUE,
                       family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep1cwt.rand,file = "rep1cwt.rand")

# Rep2 FLX
#Randomization
#Round 1 females
r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
FLX2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(FLX2.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2flxl.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                        random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped4,data=FLX2.rand,prior=p2,verbose=TRUE,
                        family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2flxl.rand,file = "rep2flxl.rand")

# Rep2 CFM
#Randomization
#Round 1 females
r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
CFM2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(CFM2.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2cfm.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                       random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped5,data=CFM2.rand,prior=p2,verbose=TRUE,
                       family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cfm.rand,file = "rep2cfm.rand")

# Rep2 CWT
#Randomization
#Round 1 females
r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
str(r1F)
r1F.rand <- r1F[,c(1,5,13,14)]
rand.order <- sample(nrow(r1F))
r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
r1F.rand$meanWs.f <- r1F[rand.order,24]
r1F.rand$meanDr.f <- r1F[rand.order,28]
r1F.rand$meanActive.f <- r1F[rand.order,19]
#Round 2 females
r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
str(r2F)
r2F.rand <- r2F[,c(1,5,13,14)]
rand.order <- sample(nrow(r2F))
r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
r2F.rand$meanWs.f <- r2F[rand.order,24]
r2F.rand$meanDr.f <- r2F[rand.order,28]
r2F.rand$meanActive.f <- r2F[rand.order,19]
str(r2F.rand)
#Round 3 females
r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
str(r3F)
r3F.rand <- r3F[,c(1,5,13,14)]
rand.order <- sample(nrow(r3F))
r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
r3F.rand$meanWs.f <- r3F[rand.order,24]
r3F.rand$meanDr.f <- r3F[rand.order,28]
r3F.rand$meanActive.f <- r3F[rand.order,19]
str(r3F.rand)
#Round 4 females
r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
str(r4F)
r4F.rand <- r4F[,c(1,5,13,14)]
rand.order <- sample(nrow(r4F))
r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
r4F.rand$meanWs.f <- r4F[rand.order,24]
r4F.rand$meanDr.f <- r4F[rand.order,28]
r4F.rand$meanActive.f <- r4F[rand.order,19]
str(r4F.rand)
#Round 1 males
r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
str(r1M)
r1M.rand <- r1M[,c(1,5,13,14)]
rand.order <- sample(nrow(r1M))
r1M.rand$meanWs.m <- r1M[rand.order,22]
r1M.rand$meanDr.m <- r1M[rand.order,26]
r1M.rand$meanActive.m <- r1M[rand.order,16]
r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
#Round 2 males
r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
str(r2M)
r2M.rand <- r2M[,c(1,5,13,14)]
rand.order <- sample(nrow(r2M))
r2M.rand$meanWs.m <- r2M[rand.order,22]
r2M.rand$meanDr.m <- r2M[rand.order,26]
r2M.rand$meanActive.m <- r2M[rand.order,16]
r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
#Round 3 males
r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
str(r3M)
r3M.rand <- r3M[,c(1,5,13,14)]
rand.order <- sample(nrow(r3M))
r3M.rand$meanWs.m <- r3M[rand.order,22]
r3M.rand$meanDr.m <- r3M[rand.order,26]
r3M.rand$meanActive.m <- r3M[rand.order,16]
r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
#Round 4 males
r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
str(r4M)
r4M.rand <- r4M[,c(1,5,13,14)]
rand.order <- sample(nrow(r4M))
r4M.rand$meanWs.m <- r4M[rand.order,22]
r4M.rand$meanDr.m <- r4M[rand.order,26]
r4M.rand$meanActive.m <- r4M[rand.order,16]
r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
#Combine randomized subsets
CWT2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
str(CWT2.rand)

#Analysis
p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
rep2cwt.rand<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                       random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped6,data=CWT2.rand,prior=p2,verbose=TRUE,
                       family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=2150000,thin=2000,burnin=150000) 
#For shorter run time:
#family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=11000,thin=100,burnin=1000) 
saveRDS(rep2cwt.rand,file = "rep2cwt.rand")

#Check when error rates levels out
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\randomized multivariate matrices")
rep1flxl.rand<-readRDS("rep1flx.rand")
rep1cfm.rand<-readRDS("rep1cfm.rand")
rep1cwt.rand<-readRDS("rep1cwt.rand")
rep2flxl.rand<-readRDS("rep2flx.rand")
rep2cfm.rand<-readRDS("rep2cfm.rand")
rep2cwt.rand<-readRDS("rep2cwt.rand")

#FLX1
str(rep1flx.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep1flx.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep1flx.rand.SEtest$Ws.m[i] <- mean(rep1flx.rand$Sol[1:i,1] - 1)
  rep1flx.rand.SEtest$SE.Ws.m[i] <- (sd(rep1flx.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep1flx.rand.SEtest$Dr.m[i] <- mean(rep1flx.rand$Sol[1:i,2] - 1)
  rep1flx.rand.SEtest$SE.Dr.m[i] <- (sd(rep1flx.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep1flx.rand.SEtest$LA.m[i] <- mean(rep1flx.rand$Sol[1:i,3] - 1)
  rep1flx.rand.SEtest$SE.LA.m[i] <- (sd(rep1flx.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep1flx.rand.SEtest$Ws.f[i] <- mean(rep1flx.rand$Sol[1:i,1] - 1)
  rep1flx.rand.SEtest$SE.Ws.f[i] <- (sd(rep1flx.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep1flx.rand.SEtest$Dr.f[i] <- mean(rep1flx.rand$Sol[1:i,2] - 1)
  rep1flx.rand.SEtest$SE.Dr.f[i] <- (sd(rep1flx.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep1flx.rand.SEtest$LA.f[i] <- mean(rep1flx.rand$Sol[1:i,3] - 1)
  rep1flx.rand.SEtest$SE.LA.f[i] <- (sd(rep1flx.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep1flx.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep1flx.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep1flx.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep1flx.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep1flx.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep1flx.rand.SEtest, type="l")

#CFM1
str(rep1cfm.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep1cfm.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep1cfm.rand.SEtest$Ws.m[i] <- mean(rep1cfm.rand$Sol[1:i,1] - 1)
  rep1cfm.rand.SEtest$SE.Ws.m[i] <- (sd(rep1cfm.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep1cfm.rand.SEtest$Dr.m[i] <- mean(rep1cfm.rand$Sol[1:i,2] - 1)
  rep1cfm.rand.SEtest$SE.Dr.m[i] <- (sd(rep1cfm.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep1cfm.rand.SEtest$LA.m[i] <- mean(rep1cfm.rand$Sol[1:i,3] - 1)
  rep1cfm.rand.SEtest$SE.LA.m[i] <- (sd(rep1cfm.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep1cfm.rand.SEtest$Ws.f[i] <- mean(rep1cfm.rand$Sol[1:i,1] - 1)
  rep1cfm.rand.SEtest$SE.Ws.f[i] <- (sd(rep1cfm.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep1cfm.rand.SEtest$Dr.f[i] <- mean(rep1cfm.rand$Sol[1:i,2] - 1)
  rep1cfm.rand.SEtest$SE.Dr.f[i] <- (sd(rep1cfm.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep1cfm.rand.SEtest$LA.f[i] <- mean(rep1cfm.rand$Sol[1:i,3] - 1)
  rep1cfm.rand.SEtest$SE.LA.f[i] <- (sd(rep1cfm.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep1cfm.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep1cfm.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep1cfm.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep1cfm.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep1cfm.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep1cfm.rand.SEtest, type="l")

#CWT1
str(rep1cwt.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep1cwt.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep1cwt.rand.SEtest$Ws.m[i] <- mean(rep1cwt.rand$Sol[1:i,1] - 1)
  rep1cwt.rand.SEtest$SE.Ws.m[i] <- (sd(rep1cwt.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep1cwt.rand.SEtest$Dr.m[i] <- mean(rep1cwt.rand$Sol[1:i,2] - 1)
  rep1cwt.rand.SEtest$SE.Dr.m[i] <- (sd(rep1cwt.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep1cwt.rand.SEtest$LA.m[i] <- mean(rep1cwt.rand$Sol[1:i,3] - 1)
  rep1cwt.rand.SEtest$SE.LA.m[i] <- (sd(rep1cwt.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep1cwt.rand.SEtest$Ws.f[i] <- mean(rep1cwt.rand$Sol[1:i,1] - 1)
  rep1cwt.rand.SEtest$SE.Ws.f[i] <- (sd(rep1cwt.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep1cwt.rand.SEtest$Dr.f[i] <- mean(rep1cwt.rand$Sol[1:i,2] - 1)
  rep1cwt.rand.SEtest$SE.Dr.f[i] <- (sd(rep1cwt.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep1cwt.rand.SEtest$LA.f[i] <- mean(rep1cwt.rand$Sol[1:i,3] - 1)
  rep1cwt.rand.SEtest$SE.LA.f[i] <- (sd(rep1cwt.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep1cwt.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep1cwt.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep1cwt.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep1cwt.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep1cwt.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep1cwt.rand.SEtest, type="l")

#FLX2
str(rep2flx.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep2flx.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep2flx.rand.SEtest$Ws.m[i] <- mean(rep2flx.rand$Sol[1:i,1] - 1)
  rep2flx.rand.SEtest$SE.Ws.m[i] <- (sd(rep2flx.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep2flx.rand.SEtest$Dr.m[i] <- mean(rep2flx.rand$Sol[1:i,2] - 1)
  rep2flx.rand.SEtest$SE.Dr.m[i] <- (sd(rep2flx.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep2flx.rand.SEtest$LA.m[i] <- mean(rep2flx.rand$Sol[1:i,3] - 1)
  rep2flx.rand.SEtest$SE.LA.m[i] <- (sd(rep2flx.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep2flx.rand.SEtest$Ws.f[i] <- mean(rep2flx.rand$Sol[1:i,1] - 1)
  rep2flx.rand.SEtest$SE.Ws.f[i] <- (sd(rep2flx.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep2flx.rand.SEtest$Dr.f[i] <- mean(rep2flx.rand$Sol[1:i,2] - 1)
  rep2flx.rand.SEtest$SE.Dr.f[i] <- (sd(rep2flx.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep2flx.rand.SEtest$LA.f[i] <- mean(rep2flx.rand$Sol[1:i,3] - 1)
  rep2flx.rand.SEtest$SE.LA.f[i] <- (sd(rep2flx.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep2flx.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep2flx.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep2flx.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep2flx.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep2flx.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep2flx.rand.SEtest, type="l")

#CFM2
str(rep2cfm.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep2cfm.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep2cfm.rand.SEtest$Ws.m[i] <- mean(rep2cfm.rand$Sol[1:i,1] - 1)
  rep2cfm.rand.SEtest$SE.Ws.m[i] <- (sd(rep2cfm.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep2cfm.rand.SEtest$Dr.m[i] <- mean(rep2cfm.rand$Sol[1:i,2] - 1)
  rep2cfm.rand.SEtest$SE.Dr.m[i] <- (sd(rep2cfm.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep2cfm.rand.SEtest$LA.m[i] <- mean(rep2cfm.rand$Sol[1:i,3] - 1)
  rep2cfm.rand.SEtest$SE.LA.m[i] <- (sd(rep2cfm.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep2cfm.rand.SEtest$Ws.f[i] <- mean(rep2cfm.rand$Sol[1:i,1] - 1)
  rep2cfm.rand.SEtest$SE.Ws.f[i] <- (sd(rep2cfm.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep2cfm.rand.SEtest$Dr.f[i] <- mean(rep2cfm.rand$Sol[1:i,2] - 1)
  rep2cfm.rand.SEtest$SE.Dr.f[i] <- (sd(rep2cfm.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep2cfm.rand.SEtest$LA.f[i] <- mean(rep2cfm.rand$Sol[1:i,3] - 1)
  rep2cfm.rand.SEtest$SE.LA.f[i] <- (sd(rep2cfm.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep2cfm.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep2cfm.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep2cfm.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep2cfm.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep2cfm.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep2cfm.rand.SEtest, type="l")

#CWT1
str(rep2cwt.rand$Sol)
samples <- c(1:1000)
Ws.m <- rep(NA,1000)
SE.Ws.m <- rep(NA,1000)
Dr.m <- rep(NA,1000)
SE.Dr.m <- rep(NA,1000)
LA.m <- rep(NA,1000)
SE.LA.m <- rep(NA,1000)
Ws.f <- rep(NA,1000)
SE.Ws.f <- rep(NA,1000)
Dr.f <- rep(NA,1000)
SE.Dr.f <- rep(NA,1000)
LA.f <- rep(NA,1000)
SE.LA.f <- rep(NA,1000)

rep2cwt.rand.SEtest <- data.frame(samples,Ws.m,SE.Ws.m,Dr.m,SE.Dr.m,LA.m,SE.LA.m,Ws.f,SE.Ws.f,Dr.f,SE.Dr.f,LA.f,SE.LA.f)

#Creating plot using SE
par(mfrow=c(2,3))
for (i in 2:1000){
  rep2cwt.rand.SEtest$Ws.m[i] <- mean(rep2cwt.rand$Sol[1:i,1] - 1)
  rep2cwt.rand.SEtest$SE.Ws.m[i] <- (sd(rep2cwt.rand$Sol[1:i,1] - 1))/(sqrt(i))
  rep2cwt.rand.SEtest$Dr.m[i] <- mean(rep2cwt.rand$Sol[1:i,2] - 1)
  rep2cwt.rand.SEtest$SE.Dr.m[i] <- (sd(rep2cwt.rand$Sol[1:i,2] - 1))/(sqrt(i))
  rep2cwt.rand.SEtest$LA.m[i] <- mean(rep2cwt.rand$Sol[1:i,3] - 1)
  rep2cwt.rand.SEtest$SE.LA.m[i] <- (sd(rep2cwt.rand$Sol[1:i,3] - 1))/(sqrt(i))
  rep2cwt.rand.SEtest$Ws.f[i] <- mean(rep2cwt.rand$Sol[1:i,1] - 1)
  rep2cwt.rand.SEtest$SE.Ws.f[i] <- (sd(rep2cwt.rand$Sol[1:i,4] - 1))/(sqrt(i))
  rep2cwt.rand.SEtest$Dr.f[i] <- mean(rep2cwt.rand$Sol[1:i,2] - 1)
  rep2cwt.rand.SEtest$SE.Dr.f[i] <- (sd(rep2cwt.rand$Sol[1:i,5] - 1))/(sqrt(i))
  rep2cwt.rand.SEtest$LA.f[i] <- mean(rep2cwt.rand$Sol[1:i,3] - 1)
  rep2cwt.rand.SEtest$SE.LA.f[i] <- (sd(rep2cwt.rand$Sol[1:i,6] - 1))/(sqrt(i))
}
plot(SE.Ws.m~samples,data=rep2cwt.rand.SEtest, type="l")
plot(SE.Dr.m~samples,data=rep2cwt.rand.SEtest, type="l")
plot(SE.LA.m~samples,data=rep2cwt.rand.SEtest, type="l")
plot(SE.Ws.f~samples,data=rep2cwt.rand.SEtest, type="l")
plot(SE.Dr.f~samples,data=rep2cwt.rand.SEtest, type="l")
plot(SE.LA.f~samples,data=rep2cwt.rand.SEtest, type="l")

#Some variation among populations and traits, but 200 samples is enough to have passed the "elbow" on the SE plots.

##################################################################################

#Create randomized matrices
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

data <- read.csv("data.csv")
str(data)
data<-as.data.frame(data) # 6858 obs. of  14 variables
names(data)[1]<-"animal"
for(x in c(1:5,12:13))data[,x]<-as.factor(data[,x])
data$mLA.active<-ifelse(data$sex=="male",data$LA_active,NA)
meanActive.m<-mean(data$mLA.activ, na.rm=TRUE)
data$meanActive.m<-data$mLA.active/meanActive.m #mean standardization

data$mLA.passive<-ifelse(data$sex=="male",data$LA_passive,NA)
data$fLA.active<-ifelse(data$sex=="female",data$LA_active,NA)
meanActive.f<-mean(data$fLA.activ, na.rm=TRUE)
data$meanActive.f<-data$fLA.active/meanActive.f

data$fLA.passive<-ifelse(data$sex=="female",data$LA_passive,NA)
data$mWS<-ifelse(data$sex=="male",data$WS,NA)
meanWs.m<-mean(data$mWS,na.rm=TRUE)
data$meanWs.m<-data$mWS/meanWs.m

data$fWS<-ifelse(data$sex=="female",data$WS,NA)
meanWS.f<-mean(data$fWS,na.rm=TRUE)
data$meanWs.f<-data$fWS/meanWS.f

data$mDR<-ifelse(data$sex=="male",data$DR ,NA)
meanDR.m<-mean(data$mDR,na.rm=TRUE)
data$meanDr.m<-data$mDR/meanDR.m

data$fDR<-ifelse(data$sex=="female",data$DR,NA)
meanDR.f<-mean(data$fDR,na.rm=TRUE)
data$meanDr.f<-data$fDR/meanDR.f

rep1<-data[data$replicat=="rep1",]
FLX1<-rep1[rep1$treatment=="FLX",]
CFM1<-rep1[rep1$treatment=="CFM",]
CWT1<-rep1[rep1$treatment=="CWT",]

rep2<-data[data$replicat=="rep2",]
FLX2<-rep2[rep2$treatment=="FLX",]
CFM2<-rep2[rep2$treatment=="CFM",]
CWT2<-rep2[rep2$treatment=="CWT",]

rep1FLX_p <- read.csv("rep1FLX.p.csv")
ped1<-rep1FLX_p
head(ped1)
str(ped1)
ped1<-as.data.frame(ped1)

rep1CFM_p <- read.csv("rep1CFM.p.csv")
ped2<-rep1CFM_p
head(ped2)
str(ped2)
ped2<-as.data.frame(ped2)

rep1CWT_p <- read.csv("rep1CWT.p.csv")
ped3<-rep1CWT_p
head(ped3)
str(ped3)
ped3<-as.data.frame(ped3)

rep2FLX_p <- read.csv("rep2FLX.p.csv")
ped4<-rep2FLX_p
head(ped4)
str(ped4)
ped4<-as.data.frame(ped4)

rep2CFM_p <- read.csv("rep2CFM.p.csv")
ped5<-rep2CFM_p
head(ped5)
str(ped5)
ped5<-as.data.frame(ped5)

rep2CWT_p <- read.csv("rep2CWT.p.csv")
ped6<-rep2CWT_p
head(ped6)
str(ped6)
ped6<-as.data.frame(ped6)

library(MCMCglmm)

#Create random FLX1 matrices
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\randomized multivariate matrices")
filenames.flx1 <- c("rep1flx.rand.1","rep1flx.rand.2","rep1flx.rand.3","rep1flx.rand.4","rep1flx.rand.5","rep1.flx.rand.6")
for (i in 1:6) {
  r1F <- FLX1[FLX1$round=="1" & FLX1$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- FLX1[FLX1$round=="2" & FLX1$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- FLX1[FLX1$round=="3" & FLX1$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- FLX1[FLX1$round=="4" & FLX1$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- FLX1[FLX1$round=="1" & FLX1$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- FLX1[FLX1$round=="2" & FLX1$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- FLX1[FLX1$round=="3" & FLX1$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- FLX1[FLX1$round=="4" & FLX1$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  FLX1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep1flx.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                          random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=FLX1.rand,prior=p2,verbose=TRUE,
                          family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep1flx.rand.red,file = filenames.flx1[i])
}
#nitt=550000,thin=2000,burnin=150000
#Create random CFM1 matrices
filenames.cfm1 <- c("rep1cfm.rand.1","rep1cfm.rand.2","rep1cfm.rand.3","rep1cfm.rand.4","rep1cfm.rand.5","rep1.cfm.rand.6")
for (i in 1:6) {
  r1F <- CFM1[CFM1$round=="1" & CFM1$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- CFM1[CFM1$round=="2" & CFM1$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- CFM1[CFM1$round=="3" & CFM1$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- CFM1[CFM1$round=="4" & CFM1$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- CFM1[CFM1$round=="1" & CFM1$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- CFM1[CFM1$round=="2" & CFM1$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- CFM1[CFM1$round=="3" & CFM1$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- CFM1[CFM1$round=="4" & CFM1$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  CFM1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep1cfm.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                             random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=CFM1.rand,prior=p2,verbose=TRUE,
                             family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep1cfm.rand.red,file = filenames.cfm1[i])
}

#Create random CWT1 matrices
filenames.cwt1 <- c("rep1cwt.rand.1","rep1cwt.rand.2","rep1cwt.rand.3","rep1cwt.rand.4","rep1cwt.rand.5","rep1.cwt.rand.6")
for (i in 1:6) {
  r1F <- CWT1[CWT1$round=="1" & CWT1$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- CWT1[CWT1$round=="2" & CWT1$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- CWT1[CWT1$round=="3" & CWT1$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- CWT1[CWT1$round=="4" & CWT1$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- CWT1[CWT1$round=="1" & CWT1$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- CWT1[CWT1$round=="2" & CWT1$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- CWT1[CWT1$round=="3" & CWT1$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- CWT1[CWT1$round=="4" & CWT1$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  CWT1.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep1cwt.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                             random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=CWT1.rand,prior=p2,verbose=TRUE,
                             family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep1cwt.rand.red,file = filenames.cwt1[i])
}
#Create random FLX2 matrices
filenames.flx2 <- c("rep2flx.rand.1","rep2flx.rand.2","rep2flx.rand.3","rep2flx.rand.4","rep2flx.rand.5","rep2.flx.rand.6")
for (i in 1:6) {
  r1F <- FLX2[FLX2$round=="1" & FLX2$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- FLX2[FLX2$round=="2" & FLX2$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- FLX2[FLX2$round=="3" & FLX2$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- FLX2[FLX2$round=="4" & FLX2$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- FLX2[FLX2$round=="1" & FLX2$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- FLX2[FLX2$round=="2" & FLX2$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- FLX2[FLX2$round=="3" & FLX2$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- FLX2[FLX2$round=="4" & FLX2$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  FLX2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep2flx.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                             random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=FLX2.rand,prior=p2,verbose=TRUE,
                             family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep2flx.rand.red,file = filenames.flx2[i])
}

#Create random CFM2 matrices
filenames.cfm2 <- c("rep2cfm.rand.1","rep2cfm.rand.2","rep2cfm.rand.3","rep2cfm.rand.4","rep2cfm.rand.5","rep2.cfm.rand.6")
for (i in 1:6) {
  r1F <- CFM2[CFM2$round=="1" & CFM2$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- CFM2[CFM2$round=="2" & CFM2$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- CFM2[CFM2$round=="3" & CFM2$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- CFM2[CFM2$round=="4" & CFM2$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- CFM2[CFM2$round=="1" & CFM2$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- CFM2[CFM2$round=="2" & CFM2$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- CFM2[CFM2$round=="3" & CFM2$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- CFM2[CFM2$round=="4" & CFM2$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  CFM2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep2cfm.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                             random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=CFM2.rand,prior=p2,verbose=TRUE,
                             family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep2cfm.rand.red,file = filenames.cfm2[i])
}

#Create random CWT2 matrices
filenames.cwt2 <- c("rep2cwt.rand.1","rep2cwt.rand.2","rep2cwt.rand.3","rep2cwt.rand.4","rep2cwt.rand.5","rep2.cwt.rand.6")
for (i in 1:6) {
  r1F <- CWT2[CWT2$round=="1" & CWT2$sex=="female",]
  str(r1F)
  r1F.rand <- r1F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1F))
  r1F.rand$meanWs.m <- rep(NA, nrow(r1F))
  r1F.rand$meanDr.m <- rep(NA, nrow(r1F))
  r1F.rand$meanActive.m <- rep(NA, nrow(r1F))
  r1F.rand$meanWs.f <- r1F[rand.order,24]
  r1F.rand$meanDr.f <- r1F[rand.order,28]
  r1F.rand$meanActive.f <- r1F[rand.order,19]
  r2F <- CWT2[CWT2$round=="2" & CWT2$sex=="female",]
  str(r2F)
  r2F.rand <- r2F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2F))
  r2F.rand$meanWs.m <- rep(NA, nrow(r2F))
  r2F.rand$meanDr.m <- rep(NA, nrow(r2F))
  r2F.rand$meanActive.m <- rep(NA, nrow(r2F))
  r2F.rand$meanWs.f <- r2F[rand.order,24]
  r2F.rand$meanDr.f <- r2F[rand.order,28]
  r2F.rand$meanActive.f <- r2F[rand.order,19]
  str(r2F.rand)
  r3F <- CWT2[CWT2$round=="3" & CWT2$sex=="female",]
  str(r3F)
  r3F.rand <- r3F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3F))
  r3F.rand$meanWs.m <- rep(NA, nrow(r3F))
  r3F.rand$meanDr.m <- rep(NA, nrow(r3F))
  r3F.rand$meanActive.m <- rep(NA, nrow(r3F))
  r3F.rand$meanWs.f <- r3F[rand.order,24]
  r3F.rand$meanDr.f <- r3F[rand.order,28]
  r3F.rand$meanActive.f <- r3F[rand.order,19]
  str(r3F.rand)
  r4F <- CWT2[CWT2$round=="4" & CWT2$sex=="female",]
  str(r4F)
  r4F.rand <- r4F[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4F))
  r4F.rand$meanWs.m <- rep(NA, nrow(r4F))
  r4F.rand$meanDr.m <- rep(NA, nrow(r4F))
  r4F.rand$meanActive.m <- rep(NA, nrow(r4F))
  r4F.rand$meanWs.f <- r4F[rand.order,24]
  r4F.rand$meanDr.f <- r4F[rand.order,28]
  r4F.rand$meanActive.f <- r4F[rand.order,19]
  str(r4F.rand)
  r1M <- CWT2[CWT2$round=="1" & CWT2$sex=="male",]
  str(r1M)
  r1M.rand <- r1M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r1M))
  r1M.rand$meanWs.m <- r1M[rand.order,22]
  r1M.rand$meanDr.m <- r1M[rand.order,26]
  r1M.rand$meanActive.m <- r1M[rand.order,16]
  r1M.rand$meanWs.f <- rep(NA, nrow(r1M))
  r1M.rand$meanDr.f <- rep(NA, nrow(r1M))
  r1M.rand$meanActive.f <- rep(NA, nrow(r1M))
  r2M <- CWT2[CWT2$round=="2" & CWT2$sex=="male",]
  str(r2M)
  r2M.rand <- r2M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r2M))
  r2M.rand$meanWs.m <- r2M[rand.order,22]
  r2M.rand$meanDr.m <- r2M[rand.order,26]
  r2M.rand$meanActive.m <- r2M[rand.order,16]
  r2M.rand$meanWs.f <- rep(NA, nrow(r2M))
  r2M.rand$meanDr.f <- rep(NA, nrow(r2M))
  r2M.rand$meanActive.f <- rep(NA, nrow(r2M))
  r3M <- CWT2[CWT2$round=="3" & CWT2$sex=="male",]
  str(r3M)
  r3M.rand <- r3M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r3M))
  r3M.rand$meanWs.m <- r3M[rand.order,22]
  r3M.rand$meanDr.m <- r3M[rand.order,26]
  r3M.rand$meanActive.m <- r3M[rand.order,16]
  r3M.rand$meanWs.f <- rep(NA, nrow(r3M))
  r3M.rand$meanDr.f <- rep(NA, nrow(r3M))
  r3M.rand$meanActive.f <- rep(NA, nrow(r3M))
  r4M <- CWT2[CWT2$round=="4" & CWT2$sex=="male",]
  str(r4M)
  r4M.rand <- r4M[,c(1,5,13,14)]
  rand.order <- sample(nrow(r4M))
  r4M.rand$meanWs.m <- r4M[rand.order,22]
  r4M.rand$meanDr.m <- r4M[rand.order,26]
  r4M.rand$meanActive.m <- r4M[rand.order,16]
  r4M.rand$meanWs.f <- rep(NA, nrow(r4M))
  r4M.rand$meanDr.f <- rep(NA, nrow(r4M))
  r4M.rand$meanActive.f <- rep(NA, nrow(r4M))
  CWT2.rand <- rbind(r1F.rand,r2F.rand,r3F.rand,r4F.rand,r1M.rand,r2M.rand,r3M.rand,r4M.rand)
  p2<-list(G=list(G1=list(V=diag(6),nu=0.002),G2=list(V=diag(6),nu=0.002)),R=list(V=diag(6),nu=0.002))
  rep2cwt.rand.red<-MCMCglmm(cbind(meanWs.m,meanDr.m,meanActive.m,meanWs.f,meanDr.f,meanActive.f)~trait-1+trait:round,
                             random=~us(trait):animal+us(trait):round.sec,rcov=~us(trait):units,pedigree=ped1,data=CWT2.rand,prior=p2,verbose=TRUE,
                             family=c("gaussian","gaussian", "gaussian","gaussian","gaussian", "gaussian"),nitt=550000,thin=2000,burnin=150000) 
  saveRDS(rep2cwt.rand.red,file = filenames.cwt2[i])
}

#Randomized combinations of matrices
#Make sure to create covtensor function first (lines 122-200)

#Create random Garrays
set.seed(2112)

#number of MCMC samples
MCMCsamp <- 1000 
#number of traits 
n <- 6 
#number of matrices to compare
m <- 6
#number of random effects specified in the model. In our analyses these were animal, section and residual effects.
r <- 3
#trait names
traitnames <- c("Ws.m","Dr.m","Active.m","Ws.f","Dr.f","Active.f") 
#matrix labels 
Gnames <- c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2")

#Create empty data frame för results, only first 5 tensors worth looking at.
rand.tensors <- data.frame("E1"=rep(NA,1000),"E2"=rep(NA,1000),"E3"=rep(NA,1000),"E4"=rep(NA,1000),"E5"=rep(NA,1000))

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\randomized multivariate matrices")
filenames.flx1 <- c("rep1flx.rand.1","rep1flx.rand.2","rep1flx.rand.3","rep1flx.rand.4","rep1flx.rand.5","rep1.flx.rand.6")
filenames.cfm1 <- c("rep1cfm.rand.1","rep1cfm.rand.2","rep1cfm.rand.3","rep1cfm.rand.4","rep1cfm.rand.5","rep1.cfm.rand.6")
filenames.cwt1 <- c("rep1cwt.rand.1","rep1cwt.rand.2","rep1cwt.rand.3","rep1cwt.rand.4","rep1cwt.rand.5","rep1.cwt.rand.6")
filenames.flx2 <- c("rep2flx.rand.1","rep2flx.rand.2","rep2flx.rand.3","rep2flx.rand.4","rep2flx.rand.5","rep2.flx.rand.6")
filenames.cfm2 <- c("rep2cfm.rand.1","rep2cfm.rand.2","rep2cfm.rand.3","rep2cfm.rand.4","rep2cfm.rand.5","rep2.cfm.rand.6")
filenames.cwt2 <- c("rep2cwt.rand.1","rep2cwt.rand.2","rep2cwt.rand.3","rep2cwt.rand.4","rep2cwt.rand.5","rep2.cwt.rand.6")

#Start loop here
for (k in 1:1000) {
  #Select random combination of matrices
  model.flx1.rand<-readRDS(sample(filenames.flx1, 1))
  VC_rep1flx.rand <- model.flx1.rand$VCV
  model.cfm1.rand<-readRDS(sample(filenames.cfm1, 1))
  VC_rep1cfm.rand <- model.cfm1.rand$VCV
  model.cwt1.rand<-readRDS(sample(filenames.cwt1, 1))
  VC_rep1cwt.rand <- model.cwt1.rand$VCV
  model.flx2.rand<-readRDS(sample(filenames.flx2, 1))
  VC_rep2flx.rand <- model.flx2.rand$VCV
  model.cfm2.rand<-readRDS(sample(filenames.cfm2, 1))
  VC_rep2cfm.rand <- model.cfm2.rand$VCV
  model.cwt2.rand<-readRDS(sample(filenames.cwt2, 1))
  VC_rep2cwt.rand <- model.cwt2.rand$VCV
  #empty array 
  MCMCarray <- array(NA,c(MCMCsamp,(n^2)*r,m)) 
  #G1 stored as the 1st element of dim[3] 
  MCMCarray[,,1] <- as.matrix(VC_rep1flx.rand)
  #G2 stored as the 2nd element of dim[3]
  MCMCarray[,,2] <- as.matrix(VC_rep1cfm.rand)
  #G3 stored as the 3rd element of dim[3]
  MCMCarray[,,3] <- as.matrix(VC_rep1cwt.rand)
  #G4 stored as the 4th element of dim[3]
  MCMCarray[,,4] <- as.matrix(VC_rep2flx.rand) 
  #G5 stored as the 5th element of dim[3]
  MCMCarray[,,5] <- as.matrix(VC_rep2cfm.rand) 
  #G6 stored as the 6th element of dim[3]
  MCMCarray[,,6] <- as.matrix(VC_rep2cwt.rand) 
  #Reshaping the array and standardizing G
  rand.Garray <- array(NA,c(n,n,m,MCMCsamp))
  dimnames(rand.Garray) <- list(traitnames,traitnames,Gnames)
  rand.Parray <- array(NA,c(n,n,m,MCMCsamp))
  dimnames(rand.Parray) <- list(traitnames,traitnames,Gnames)
  for (i in 1:m){
    for (j in 1:MCMCsamp){
      G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
      CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
      R <- matrix(MCMCarray[j,(((n^2)*2)+1):((n^2)*3),i],ncol= n)
      rand.Garray[,,i,j] <- G
      rand.Parray[,,i,j] <- G + CE + R
    }
  }
  MCMC.covtensor.rand <- covtensor(rand.Garray)
  #Output posterior mode for each tensor
  rand.tensors[k,] <- posterior.mode(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[,1:nnonzero]))
}

HPD.eT.val <- cbind(mean.obs=posterior.mode(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero])), HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[,1:nnonzero]), prob=0.95), mean.rand=posterior.mode(as.mcmc(rand.tensors)), HPDinterval(as.mcmc(rand.tensors), prob=0.95))
round(HPD.eT.val, 3)

write.csv(round(HPD.eT.val, 3), "Randomization test output.csv")

HPD.eT.val <- read.csv("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\randomized multivariate matrices\\Randomization test output.csv")
str(HPD.eT.val)
HPD.eT.val$X <- factor(HPD.eT.val$X)

#Figure S4
par(mfrow=c(1,1))
bp<-barplot(t(as.matrix(HPD.eT.val[,c(2,5)])), legend=c("Observed","Randomized"), ylab="Alpha", ylim=c(0, 0.4), beside=TRUE, names.arg=(c("E1","E2","E3","E4","E5")))
arrows(bp, t(as.matrix(HPD.eT.val[,c(3,6)])), bp, t(as.matrix(HPD.eT.val[,c(4,7)])), code=3, angle=90)
#NB - The last two error bars for the randomized data have length zero and therefore produce a warning message.

# Figure 4
flx1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,1,]))
flx1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[1,1,]))
cfm1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,1,]))
cfm1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[2,1,]))
cwt1.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,1,]))
cwt1.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[3,1,]))
flx2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,1,]))
flx2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[4,1,]))
cfm2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,1,]))
cfm2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[5,1,]))
cwt2.mean.E1 <- posterior.mode(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,1,]))
cwt2.CI.E1 <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[6,1,]))

b <- barplot(c(flx1.mean.E1, cfm1.mean.E1, cwt1.mean.E1),
             ylab="", xlab="", names=c("FLX","CFM","CWT"), ylim=c(-2,0.25),yaxt = "n",cex.names = 1.5,
             col=c("#0DA7C9", "#FBD00A","#009966"))
arrows(b, c(flx1.CI.E1[1], cfm1.CI.E1[1], cwt1.CI.E1[1]),
       b, c(flx1.CI.E1[2], cfm1.CI.E1[2], cwt1.CI.E1[2]),code=3, angle=90,lwd =2,
       col=c("#6fddf5", "#d2ad03","#4dcf88"))
title(main ="Rep1",line=-0.5,adj = 0,cex.main=1.6,font.main=1)
title(ylab = "E1", line = 2.6,cex.lab=1.5)
axis(2,at=seq(-2,0,by=0.5),las=2, cex.axis=1.3)

b2 <- barplot(c(flx2.mean.E1, cfm2.mean.E1, cwt2.mean.E1),
              ylab="", xlab="", names=c("FLX","CFM","CWT"), ylim=c(-2,0.25),yaxt = "n",cex.names = 1.5,
              col=c("#0DA7C9", "#FBD00A","#009966"))
arrows(b2, c(flx2.CI.E1[1], cfm2.CI.E1[1], cwt2.CI.E1[1]),
       b2, c(flx2.CI.E1[2], cfm2.CI.E1[2], cwt2.CI.E1[2]),code=3, angle=90,lwd =2,
       col=c("#6fddf5", "#d2ad03","#4dcf88"))
title(main ="Rep2",line=-0.5,adj = 0,cex.main=1.6,font.main=1)
title(ylab = "E1", line = 2.6,cex.lab=1.5)
axis(2,at=seq(-2,0,by=0.5),las=2, cex.axis=1.3)


#Figure 5

#For e21
flx1.mean.e21 <- mean(e21.proj[1,])
flx1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[1,]))
cfm1.mean.e21 <- mean(e21.proj[2,])
cfm1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[2,]))
cwt1.mean.e21 <- mean(e21.proj[3,])
cwt1.CI.e21 <- HPDinterval(as.mcmc(e21.proj[3,]))
flx2.mean.e21 <- mean(e21.proj[4,])
flx2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[4,]))
cfm2.mean.e21 <- mean(e21.proj[5,])
cfm2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[5,]))
cwt2.mean.e21 <- mean(e21.proj[6,])
cwt2.CI.e21 <- HPDinterval(as.mcmc(e21.proj[6,]))
#For e22
flx1.mean.e22 <- mean(e22.proj[1,])
flx1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[1,]))
cfm1.mean.e22 <- mean(e22.proj[2,])
cfm1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[2,]))
cwt1.mean.e22 <- mean(e22.proj[3,])
cwt1.CI.e22 <- HPDinterval(as.mcmc(e22.proj[3,]))
flx2.mean.e22 <- mean(e22.proj[4,])
flx2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[4,]))
cfm2.mean.e22 <- mean(e22.proj[5,])
cfm2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[5,]))
cwt2.mean.e22 <- mean(e22.proj[6,])
cwt2.CI.e22 <- HPDinterval(as.mcmc(e22.proj[6,]))

#Rep1
par(mar = c(4, 4, 4, 2))
plot(c(flx1.mean.e21, cfm1.mean.e21, cwt1.mean.e21),
     c(flx1.mean.e22, cfm1.mean.e22, cwt1.mean.e22), 
     xlab="", ylab="", xlim=c(0,1), ylim=c(0,1.2), pch=c(15, 16, 17),
     col=c("#0DA7C9", "#FBD00A","#009966"),cex = 2,yaxt = "n",xaxt = "n")
title(main ="Rep1",line=0.5,adj = 0,cex.main=1.4,font.main=1)
title(ylab ="e22",line = 2.4,cex.lab=1.6)
title(xlab ="e21",line = 2.4,cex.lab=1.6)
axis(1,las=1, cex.axis=1.3)
axis(2,las=2, cex.axis=1.3)

arrows(flx1.CI.e21[1], flx1.mean.e22, flx1.CI.e21[2], flx1.mean.e22, code=3, angle=90, col ="#0DA7C9",lwd =2) #FLX1 in e21
arrows(flx1.mean.e21, flx1.CI.e22[1],  flx1.mean.e21, flx1.CI.e22[2], code=3, angle=90, col ="#0DA7C9",lwd =2) #FLX1 in e22
arrows(cfm1.CI.e21[1], cfm1.mean.e22, cfm1.CI.e21[2], cfm1.mean.e22, code=3, angle=90,col ="#FBD00A",lwd =2) #CFM1 in e21
arrows(cfm1.mean.e21, cfm1.CI.e22[1],  cfm1.mean.e21, cfm1.CI.e22[2], code=3, angle=90,col ="#FBD00A",lwd =2) #CFM1 in e22
arrows(cwt1.CI.e21[1], cwt1.mean.e22, cwt1.CI.e21[2], cwt1.mean.e22, code=3, angle=90,col ="#009966",lwd =2) #CWT1 in e21
arrows(cwt1.mean.e21, cwt1.CI.e22[1],  cwt1.mean.e21, cwt1.CI.e22[2], code=3, angle=90,col ="#009966",lwd =2) #CWT1 in e22

legend(0.75,1.2, legend=c("FLX","CFM","CWT"), 
       pch=c(15, 16, 17), col=c("#0DA7C9", "#FBD00A","#009966"),
       cex = 1.2,box.lty=0,x.intersp=0.7)

#Rep2
par(mar = c(4, 4, 4, 2))
plot(c(flx2.mean.e21, cfm2.mean.e21, cwt2.mean.e21),
     c(flx2.mean.e22, cfm2.mean.e22, cwt2.mean.e22), 
     xlab="", ylab="", xlim=c(0,1), ylim=c(0,1.2), pch=c(15, 16, 17),
     col=c("#0DA7C9", "#FBD00A","#009966"),cex = 2,yaxt = "n",xaxt = "n")
title(main ="Rep2",line=0.5,adj = 0,cex.main=1.4,font.main=1)
title(ylab ="e22",line = 2.4,cex.lab=1.6)
title(xlab ="e21",line = 2.4,cex.lab=1.6)
axis(1,las=1, cex.axis=1.3)
axis(2,las=2, cex.axis=1.3)

arrows(flx2.CI.e21[1], flx2.mean.e22, flx2.CI.e21[2], flx2.mean.e22, code=3, angle=90, col ="#0DA7C9",lwd =2) #FLX1 in e21
arrows(flx2.mean.e21, flx2.CI.e22[1],  flx2.mean.e21, flx2.CI.e22[2], code=3, angle=90, col ="#0DA7C9",lwd =2) #FLX1 in e22
arrows(cfm2.CI.e21[1], cfm2.mean.e22, cfm2.CI.e21[2], cfm2.mean.e22, code=3, angle=90,col ="#FBD00A",lwd =2) #CFM1 in e21
arrows(cfm2.mean.e21, cfm2.CI.e22[1],  cfm2.mean.e21, cfm2.CI.e22[2], code=3, angle=90,col ="#FBD00A",lwd =2) #CFM1 in e22
arrows(cwt2.CI.e21[1], cwt2.mean.e22, cwt2.CI.e21[2], cwt2.mean.e22, code=3, angle=90,col ="#009966",lwd =2) #CWT1 in e21
arrows(cwt2.mean.e21, cwt2.CI.e22[1],  cwt2.mean.e21, cwt1.CI.e22[2], code=3, angle=90,col ="#009966",lwd =2) #CWT1 in e22

legend(0.75,1.24, legend=c("FLX","CFM","CWT"), 
       pch=c(15, 16, 17), col=c("#0DA7C9", "#FBD00A","#009966"),
       cex = 1.2,box.lty=0,x.intersp=0.7,y.intersp = 1)

##############################################################

#Random skewers analysis

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")

model.flx1<-readRDS("rep1flx")
VC_rep1flx <- model.flx1$VCV
model.cfm1<-readRDS("rep1cfm")
VC_rep1cfm <- model.cfm1$VCV
model.cwt1<-readRDS("rep1cwt")
VC_rep1cwt <- model.cwt1$VCV
model.flx2<-readRDS("rep2flx")
VC_rep2flx <- model.flx2$VCV
model.cfm2<-readRDS("rep2cfm")
VC_rep2cfm <- model.cfm2$VCV
model.cwt2<-readRDS("rep2cwt")
VC_rep2cwt <- model.cwt2$VCV

#Set random seed for repeatability
set.seed(2112)

#number of MCMC samples
MCMCsamp <- 1000 
#number of traits 
n <- 6 
#number of matrices to compare
m <- 6
#number of random effects specified in the model. In our analyses these were animal, section and residual effects.
r <- 3
#trait names
traitnames <- c("Ws.m","Dr.m","Active.m","Ws.f","Dr.f","Active.f") 
#matrix labels 
Gnames <- c("FLX1","CFM1","CWT1","FLX2","CFM2","CWT2")

#empty array 
MCMCarray <- array(NA,c(MCMCsamp,(n^2)*r,m)) 
#G1 stored as the 1st element of dim[3] 
MCMCarray[,,1] <- as.matrix(VC_rep1flx)
#G2 stored as the 2nd element of dim[3]
MCMCarray[,,2] <- as.matrix(VC_rep1cfm)
#G3 stored as the 3rd element of dim[3]
MCMCarray[,,3] <- as.matrix(VC_rep1cwt)
#G4 stored as the 4th element of dim[3]
MCMCarray[,,4] <- as.matrix(VC_rep2flx) 
#G5 stored as the 5th element of dim[3]
MCMCarray[,,5] <- as.matrix(VC_rep2cfm) 
#G6 stored as the 6th element of dim[3]
MCMCarray[,,6] <- as.matrix(VC_rep2cwt) 

#Reshaping the array and standardizing G
Garray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(Garray) <- list(traitnames,traitnames,Gnames)
Parray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(Parray) <- list(traitnames,traitnames,Gnames)
for (i in 1:m){
  for (j in 1:MCMCsamp){
    G <- matrix(MCMCarray[j,1:(n^2),i],ncol= n)
    CE <- matrix(MCMCarray[j,((n^2)+1):((n^2)*2),i],ncol= n)
    R <- matrix(MCMCarray[j,(((n^2)*2)+1):((n^2)*3),i],ncol= n)
    Garray[,,i,j] <- G
    Parray[,,i,j] <- G + CE + R
  }
}

#Define standardization function c.f. Hansen and Houle (2008)
#Should be used if traits are not on the same scale. The HHGarray object will contain the standardized G array. 
#I assume that HHGarray should be used anywhere below where the Garray is called.
inv.rootP <- function (P){
  rootP <- matrix(0,n, n)  
  for (i in 1:n){
    val <- eigen(P)$values
    vec <- eigen(P)$vectors
    rootP <- rootP + (vec[,i] %*% t(vec[,i]))*sqrt(val[i])
  }
  solve(rootP)
}

HHGarray <- array(NA,c(n,n,m,MCMCsamp))
for (k in 1:MCMCsamp){
  for (j in 1:m){
    P <- inv.rootP(Parray[,,j,k])
    HHGarray[,,j,k] <- P %*% Garray[,,j,k] %*% P
  }
}

#Generating randomised G matrices for hypothesis tests. 
#Import pedigree data
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")

Ped.flx1 <- read.csv("rep1FLX.p.csv")
Ped.cfm1 <- read.csv("rep1CFM.p.csv")
Ped.cwt1 <- read.csv("rep1CWT.p.csv")
Ped.flx2 <- read.csv("rep2FLX.p.csv")
Ped.cfm2 <- read.csv("rep2CFM.p.csv")
Ped.cwt2 <- read.csv("rep2CWT.p.csv")

library(MCMCglmm)

#Simulate a set of G that have been sampled from the same population, so the only dissimilarity among them is random sampling error
rand.Garray <- array(NA,c(n,n,m,MCMCsamp))
dimnames(rand.Garray) <- list(traitnames,traitnames,Gnames)
for (i in 1:MCMCsamp){
  flx1.bv<-rbv(Ped.flx1,HHGarray[,,1,i])
  cfm1.bv<-rbv(Ped.cfm1,HHGarray[,,2,i])
  cwt1.bv<-rbv(Ped.cwt1,HHGarray[,,3,i])
  flx2.bv<-rbv(Ped.flx2,HHGarray[,,4,i])
  cfm2.bv<-rbv(Ped.cfm2,HHGarray[,,5,i])
  cwt2.bv<-rbv(Ped.cwt2,HHGarray[,,6,i])
  a.pop <- cumsum(c(dim(Ped.flx1)[1],dim(Ped.cfm1)[1],dim(Ped.cwt1)[1],dim(Ped.flx2)[1],dim(Ped.cfm2)[1],dim(Ped.cwt2)[1]))
  pop.bv <- rbind(flx1.bv,cfm1.bv,cwt1.bv,flx2.bv,cfm2.bv,cwt2.bv)
  rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1],replace=F),]
  rand.Garray[,,1,i] <- cov(rand.pop.bv[1:a.pop[1],])
  rand.Garray[,,2,i] <- cov(rand.pop.bv[(a.pop[1] + 1):a.pop[2],])
  rand.Garray[,,3,i] <- cov(rand.pop.bv[(a.pop[2] + 1):a.pop[3],])
  rand.Garray[,,4,i] <- cov(rand.pop.bv[(a.pop[3] + 1):a.pop[4],])
  rand.Garray[,,5,i] <- cov(rand.pop.bv[(a.pop[4] + 1):a.pop[5],])
  rand.Garray[,,6,i] <- cov(rand.pop.bv[(a.pop[5] + 1):a.pop[6],])
}

## Method 1. Random projections through G 

#START
R.proj <- function(Gs,p,vec){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
  for (i in 1:vec){
    b <- runif(n,-1,1)
    rand.vec[i,] <- b/(sqrt(sum(b^2)))
  }
  #generate unit length random vectors  
  proj<- function(G,b) t(b) %*% G %*% (b)
  #internal function to do projection
  G.proj <- array(,c(MCMCsamp, m, vec))
  colnames(G.proj) <- dimnames(Gs)[[3]]
  for (i in 1:vec){
    G.proj[,,i]<- t(apply(Gs, 3:4, proj, b = rand.vec[i,]))
  }
  #project each random vector through each MCMC sample of each G
  prs <- cbind(rep(1:m, each = m), 1:m) 
  prs.comp <- prs[prs[,1] < prs[,2], , drop = FALSE] 
  #setting up an index for HPD comparisons
  proj.score <-matrix(,vec,((m^2 - m)/2))
  for (k in 1:vec){
    HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
    proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0) 
  }
  #for a given random vector, examine if the HPD intervals of any pair of G matrices overlap
  vec.score <-cbind(rand.vec, proj.score)
  colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
  #collate the random vectors and the outcome of their projection on the G matrices
  sig.vec <- subset(vec.score, rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0) 
  #collate just the random vectors that resulted in significant differences in variance
  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
  else{
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }  
  #eigen analysis of the R matrix
  list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
}
#END

MCMC.R.proj <- R.proj(Garray, p = 0.95, vec = 1000)

# Significant differences in genetic variance among populations (significant differences are indexed with 1’s and those not finding significant difference are indexed with 0’s) 
MCMC.R.proj$vec.score[1:n,(n+1):(n+((m^2 - m)/2))]

#  The proportion of random vectors that found significant differences in genetic variance among populations.
table(rowSums(MCMC.R.proj$vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0 )
#468/1000 found differences

# To visualize the eigenstructure of the R matrix
lapply(MCMC.R.proj$eig.R, round, digits = 3)
#$values
#[1] 0.290 0.208 0.164 0.150 0.138 0.049

#$vectors
#         e1     e2     e3     e4     e5     e6
#Ws.m     0.013  0.280  0.029  0.959  0.026 -0.023
#Dr.m     0.005 -0.884 -0.133  0.250  0.349 -0.126
#Active.m 0.833  0.074 -0.027 -0.045  0.000 -0.545
#Ws.f     0.022  0.047  0.869 -0.054  0.489  0.000
#Dr.f     0.013 -0.348  0.475  0.107 -0.799 -0.059
#Active.f 0.552 -0.103 -0.004  0.043 -0.003  0.826

# Figure S5
par(mfrow=c(2,3),mar = c(5, 5, 2, 3))

HPD1 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,1]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,1]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD1))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames, cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,1]),1:m,HPD1[,1],length=0.1,angle = 90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,1]),1:m,HPD1[,2],length=0.1,angle=90)
mtext("A",side=3,at=0,font=2)

HPD2 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,2]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,2]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD2))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames, cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,2]),1:m,HPD2[,1],length=0.1,angle=90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,2]),1:m,HPD2[,2],length=0.1,angle=90)
mtext("B",side=3,at=0,font=2)


HPD3 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,3]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,3]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD3))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames,cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,3]),1:m,HPD3[,1],length=0.1,angle = 90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,3]),1:m,HPD3[,2],length=0.1,angle=90)
mtext("C",side=3,at=0,font=2)

HPD4 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,4]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,4]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD4))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames,cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,4]),1:m,HPD4[,1],length=0.1,angle = 90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,4]),1:m,HPD4[,2],length=0.1,angle=90)
mtext("D",side=3,at=0,font=2)

HPD5 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,5]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,5]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD5))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames,cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,5]),1:m,HPD5[,1],length=0.1,angle = 90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,5]),1:m,HPD5[,2],length=0.1,angle=90)
mtext("E",side=3,at=0,font=2)

HPD6 <- HPDinterval(as.mcmc(MCMC.R.proj$G.proj[,,6]))

plot(1:m,colMeans(MCMC.R.proj$G.proj[,,6]),yaxs = "r",
     ylim=c(0,ceiling(max(HPD6))),xlab="",ylab="Va",cex.lab=1.7,pch=16,cex=1.2,xaxt="n",
     frame.plot=F,cex.axis=1.5)
axis(1,at=1:m,labels=Gnames,cex.axis=1.5,las=2)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,6]),1:m,HPD6[,1],length=0.1,angle = 90)
arrows(1:m,colMeans(MCMC.R.proj$G.proj[,,6]),1:m,HPD6[,2],length=0.1,angle=90)
mtext("F",side=3,at=0,font=2)
