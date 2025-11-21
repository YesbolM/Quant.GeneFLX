#Predicted response to selection in males when using CWT or FLX G-matrices

library(MCMCglmm)

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")

############# Rep1 ##############
# Predict using CWT G  #########

# Rep1CWt and Rep1FLX 
rep1cwt<-readRDS("rep1cwt")
rep1flx<-readRDS("rep1flx")

# Phenotypic data in females (CWT1)
pheno_CWT1 <- matrix(apply(rep1cwt$Sol[,4:(2*3)],2, mean),1,3)

# Phenotypic data in females (FLX1)
pheno_FLX1 <- matrix(apply(rep1flx$Sol[,4:(2*3)],2, mean),1,3)

# To be sure for the mean estimation
mean_pheno<-as.data.frame(rep1flx$Sol)
mean(mean_pheno$traitmeanWs.f) # 1.032637     
mean(mean_pheno$traitmeanDr.f) # 0.8601164    
mean(mean_pheno$traitmeanActive.f) # 1.648984

# 1) Get deltaZ for each of the three traits in females (Trait1FemaleFLX-Trait1FemaleCWT)
deltaZ<-pheno_FLX1-pheno_CWT1

#### 2) Get Beta selection gradients on females using only the female G matrix 
# (i.e. extract the 3x3 part of the full G that only contains female traits) from CWT

## female G matrix of rep1CWT
G.cwt1 <- matrix(apply(rep1cwt$VCV[,1:(6*6)],2, mean),6,6) #G.1
# select the G matrix in females
G.cwt1.f<-G.cwt1[4:6,4:6] # f.G

## 2) Get Beta selection gradients on females using only the female G matrix
# deltaZ<-G.cwt1.f*beta, beta=deltaZ/G.cwt1.f    <!!!> Technically, there is no such thing as matrix division.
                                       # Instead, we multiply the vector by the inverse of the matrix.
G.cwt1.f
inverse_G<-solve(G.cwt1.f) # to determine the inverse matrix 
beta=deltaZ%*%inverse_G

# 3) Add 0 for the males in the Beta vector
zero_0<-c(0,0,0)
beta_0<-c(zero_0,beta)

# 4) Predict DeltaZ in males by dividing DeltaZ by full G (i.e. DeltaZ/G) in a) CWT,   b) FLX
#   a)
G.cwt1 ## female G matrix in cwt1
deltaZ<-pheno_FLX1-pheno_CWT1 # Phenotypic differences btween Rep1FLX and Rep1CWT in females
Dz_CWT1<-beta_0%*%G.cwt1

############# Rep1 ##############
# Predict using FLX G  #########

# 1) Get deltaZ for each of the three traits in females (Trait1FemaleFLX-Trait1FemaleCWT)
deltaZ<-pheno_FLX1-pheno_CWT1

#### 2) Get Beta selection gradients on females using only the female G matrix 
# (i.e. extract the 3x3 part of the full G that only contains female traits) from FLX

## female G matrix of rep1FLX
G.flx1 <- matrix(apply(rep1flx$VCV[,1:(6*6)],2, mean),6,6) #G.1
# select the G matrix in females
G.flx1.f<-G.flx1[4:6,4:6] # f.G

## 2) Get Beta selection gradients on females using only the female G matrix
G.flx1.f
inverse_G<-solve(G.flx1.f) # to determine the inverse matrix 
beta=deltaZ%*%inverse_G

# 3) Add 0 for the males in the Beta vector
zero_0<-c(0,0,0)
beta_0<-c(zero_0,beta)

# 4) Predict DeltaZ in males by dividing DeltaZ by full G (i.e. DeltaZ/G) in b) FLX
#   a)
Dz_FLX1<-beta_0%*%G.flx1


############# Rep2 ##############
# Predict using CWT G  #########

# Rep2CWt and Rep2FLX 
rep2cwt<-readRDS("rep2cwt")
rep2flx<-readRDS("rep2flx")

# Phenotypic data in females (CWT2)
pheno_CWT2 <- matrix(apply(rep2cwt$Sol[,4:(2*3)],2, mean),1,3)

# Phenotypic data in females (FLX2)
pheno_FLX2 <- matrix(apply(rep2flx$Sol[,4:(2*3)],2, mean),1,3)

# 1) Get deltaZ for each of the three traits in females (Trait1FemaleFLX-Trait1FemaleCWT)
deltaZ<-pheno_FLX2-pheno_CWT2

#### 2) Get Beta selection gradients on females using only the female G matrix 
# (i.e. extract the 3x3 part of the full G that only contains female traits) from CWT

## female G matrix of rep2CWT
G.cwt2 <- matrix(apply(rep2cwt$VCV[,1:(6*6)],2, mean),6,6) 
# select the G matrix in females
G.cwt2.f<-G.cwt2[4:6,4:6] # f.G

## 2) Get Beta selection gradients on females using only the female G matrix
# deltaZ<-G.cwt1.f*beta, beta=deltaZ/G.cwt1.f    <!!!> Technically, there is no such thing as matrix division.
# Instead, we multiply the vector by the inverse of the matrix.
G.cwt2.f
inverse_G<-solve(G.cwt2.f) # to determine the inverse matrix 
beta=deltaZ%*%inverse_G

# 3) Add 0 for the males in the Beta vector
zero_0<-c(0,0,0)
beta_0<-c(zero_0,beta)

# 4) Predict DeltaZ in males by dividing DeltaZ by full G (i.e. DeltaZ/G) in a) CWT
#   a)
Dz_CWT2<-beta_0%*%G.cwt2


############# Rep2 ##############
# Predict using FLX G  #########

# 1) Get deltaZ for each of the three traits in females (Trait1FemaleFLX-Trait1FemaleCWT)
deltaZ<-pheno_FLX2-pheno_CWT2

#### 2) Get Beta selection gradients on females using only the female G matrix 
# (i.e. extract the 3x3 part of the full G that only contains female traits) from FLX

## female G matrix of rep2FLX
G.flx2 <- matrix(apply(rep2flx$VCV[,1:(6*6)],2, mean),6,6) #G.1
# select the G matrix in females
G.flx2.f<-G.flx2[4:6,4:6] 

## 2) Get Beta selection gradients on females using only the female G matrix
G.flx2.f
inverse_G<-solve(G.flx2.f) # to determine the inverse matrix 
beta=deltaZ%*%inverse_G

# 3) Add 0 for the males in the Beta vector
zero_0<-c(0,0,0)
beta_0<-c(zero_0,beta)

# 4) Predict DeltaZ in males by dividing DeltaZ by full G (i.e. DeltaZ/G) in b) FLX
#   a)
Dz_FLX2<-beta_0%*%G.flx2


###########################  Visualization # Figure 6 ###########################
library(ggplot2)
library(ggthemes)
library(cowplot)

# Phenotypic data in males (CWT1)
pheno_CWT1M <- matrix(apply(rep1cwt$Sol[,1:3],2, mean),1,3)

# Phenotypic data in males (FLX1)
pheno_FLX1M <- matrix(apply(rep1flx$Sol[,1:3],2, mean),1,3)

# Phenotypic data in males (CWT2)
pheno_CWT2M <- matrix(apply(rep2cwt$Sol[,1:3],2, mean),1,3)

# Phenotypic data in males (FLX2)
pheno_FLX2M <- matrix(apply(rep2flx$Sol[,1:3],2, mean),1,3)

#Wing size
x<-rep(c("PRED.CWT","PRED.FLX"),times=2)
Replicate<-rep(c("Rep1","Rep2"),each=2)
del<-c(Dz_CWT1[1],Dz_FLX1[1],Dz_CWT2[1],Dz_FLX2[1])
obs1<-pheno_FLX1M[1]-pheno_CWT1M[1]
obs2<-pheno_FLX2M[1]-pheno_CWT2M[1]
df<-data.frame(x,Replicate,del)
df$x<-factor(df$x)

wl<-ggplot(df, aes(x=x,
                       y=del,
                       group=Replicate,
                       color=Replicate))+
  geom_point(size=3)+
  geom_hline(yintercept=obs1,linetype="dashed", color = "#F8766D")+
  geom_hline(yintercept=obs2,linetype="dashed", color = "#00BFC4")+
  labs(title = "A",
       x="",
       y="DeltaZ WL") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "plain",color = "black",
                                   size = 12,angle = 30,vjust = 0.7),
        axis.text.y = element_text(face = "plain",color = "black",
                                   size = 12),
        text = element_text(size = 14))+
  theme(legend.position = "none")

wl


# DR
x<-rep(c("PRED.CWT","PRED.FLX"),times=2)
Replicate<-rep(c("Rep1","Rep2"),each=2)
del<-c(Dz_CWT1[2],Dz_FLX1[2],Dz_CWT2[2],Dz_FLX2[2])
obs1<-pheno_FLX1M[2]-pheno_CWT1M[2]
obs2<-pheno_FLX2M[2]-pheno_CWT2M[2]
df<-data.frame(x,Replicate,del)
df$x<-factor(df$x)

dr<-ggplot(df, aes(x=x,
                   y=del,
                   group=Replicate,
                   color=Replicate))+
  geom_point(size=3)+
  geom_hline(yintercept=obs1,linetype="dashed", color = "#F8766D")+
  geom_hline(yintercept=obs2,linetype="dashed", color = "#00BFC4")+
  labs(title = "B",
       x="",
       y="DeltaZ DR") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "plain",color = "black",
                                   size = 12,angle = 30,vjust = 0.7),
        axis.text.y = element_text(face = "plain",color = "black",
                                   size = 12),
        text = element_text(size = 14))+
  theme(legend.position = "none")

dr

# LA
x<-rep(c("PRED.CWT","PRED.FLX"),times=2)
Replicate<-rep(c("Rep1","Rep2"),each=2)
del<-c(Dz_CWT1[3],Dz_FLX1[3],Dz_CWT2[3],Dz_FLX2[3])
obs1<-pheno_FLX1M[3]-pheno_CWT1M[3]
obs2<-pheno_FLX2M[3]-pheno_CWT2M[3]
df<-data.frame(x,Replicate,del)
df$x<-factor(df$x)

la<-ggplot(df, aes(x=x,
                   y=del,
                   group=Replicate,
                   color=Replicate))+
  geom_point(size=3)+
  geom_hline(yintercept=obs1,linetype="dashed", color = "#F8766D")+
  geom_hline(yintercept=obs2,linetype="dashed", color = "#00BFC4")+
  labs(title = "C",
       x="",
       y="DeltaZ LA") +
  theme_classic() +
  theme(axis.text.x = element_text(face = "plain",color = "black",
                                   size = 12,angle = 30,vjust = 0.7),
        axis.text.y = element_text(face = "plain",color = "black",
                                   size = 12),
        text = element_text(size = 14))+
  theme(legend.position = "none")

la


ggdraw() +
  draw_plot(wl, x = 0.02, y = 0.3, width = 0.3, height = .4) +
  draw_plot(dr, x = 0.35, y = 0.3, width = 0.3, height = .4)+
  draw_plot(la, x = 0.68, y = 0.3, width = 0.3, height = .4)
