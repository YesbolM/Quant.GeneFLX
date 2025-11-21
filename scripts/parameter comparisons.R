#Import dataset
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript")
modest <- read.csv("all model estimates.csv", sep=";", dec=",")
str(modest)

#Correlation between univariate and bivariate models
#MCMCglmm
cor.test(modest$Uni_MCMC,modest$Bi_MCMC)
#p=4.228e-16, r=0.9275844
#WOMBAT
cor.test(modest$Uni_WOM,modest$Bi_WOM)
#P < 2.2e-16, r=0.9914726

#Correlation between methods of analysis for same type of model
#Univariate models, MCMCglmm vs WOMBAT
cor.test(modest$Uni_MCMC,modest$Uni_WOM)
#p-value < 2.2e-16, r=0.9549022
#Bivariate models, MCMCglmm vs WOMBAT
cor.test(modest$Bi_MCMC,modest$Bi_WOM)
#p-value < 2.2e-16. r=0.9432457
#Multivariate models, MCMCglmm vs WOMBAT
cor.test(modest$Multi_MCMC,modest$Multi_WOM)
#p-value = 1.384e-06, r=0.7077318

#Correlation between univariate and multivariate heritability estimates
#MCMCglmm
cor.test(modest$Uni_MCMC_h2,modest$Multi_MCMC_h2)
#p-value = 3.267e-07, r=0.7352507
#WOMBAT
cor.test(modest$Uni_WOM_h2,modest$Multi_WOM_h2)
#p-value = 2.429e-12, r=0.876617

#Extract range of estimates
aggregate(modest$mean_h2, by=list(modest$Trait), range)
#Group.1        x.1        x.2
#1     fDR 0.29000000 0.47166667
#2     fLA 0.06933333 0.21033333
#3     fWS 0.48500000 0.78333333
#4     mDR 0.24333333 0.68166667
#5     mLA 0.01583333 0.50786667
#6     mWS 0.53333333 0.69833333

#Extract mean of estimates
aggregate(modest$mean_h2, by=list(modest$Trait), mean)
# Group.1         x
#1     fDR 0.3722222
#2     fLA 0.1475667
#3     fWS 0.6366667
#4     mDR 0.4030556
#5     mLA 0.3408944
#6     mWS 0.6275000