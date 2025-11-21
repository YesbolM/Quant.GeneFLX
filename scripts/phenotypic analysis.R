                             # Phenotypic analysis #

library(lme4)    # GLMMs
library(lmerTest)
library(blmeco)  # Overdispersion test. 
library(car)     # Anova
library(emmeans) # Post-hoc test.
library(DHARMa)

#Import dataset
setwd("G:\\My Drive\\Lund\\Flies\\Yesbol\\Quant gen manuscript")
setwd("E:\\My Drive\\Lund\\Flies\\Yesbol\\Quant gen manuscript")
data <- read.csv("data.csv")
str(data)
for(x in 1:5)data[,x]<-as.factor(data[,x])
data$round<-factor(data$round)
data$section<-factor(data$section)
data$round.sec<-factor(data$round.sec)

#Subset by sex
females <- subset(data, data$sex=="female")
males <- subset(data, data$sex=="male")

#Desiccation resistance
#Calculate family means
FDR <- aggregate(females$DR, by=list(mother=females$mother,replicat=females$replicat,treatment=females$treatment,round=females$round), FUN=mean, na.action=na.rm)
FDR
str(FDR)
View(FDR)

MDR <- aggregate(males$DR, by=list(mother=males$mother,replicat=males$replicat,treatment=males$treatment,round=males$round), FUN=mean, na.action=na.rm)
MDR
str(MDR)
View(MDR)

#Analysis
modelF.dr<- lmer(x~treatment+round+(1|treatment:replicat), data=FDR)
summary(modelF.dr)
anova(modelF.dr)
emmeans(modelF.dr,list(pairwise ~treatment),adjust="tukey")
plot(modelF.dr) 
qqnorm(residuals(modelF.dr)) 
qqline(residuals(modelF.dr)) 

modelM.dr<- lmer(x~treatment+round+(1|treatment:replicat), data=MDR)
summary(modelM.dr)
anova(modelM.dr)
emmeans(modelM.dr,list(pairwise ~treatment),adjust="tukey")
plot(modelM.dr) 
qqnorm(residuals(modelM.dr)) 
qqline(residuals(modelM.dr)) 

#Wing size
#Calculate family means
FWS <- aggregate(females$WS, by=list(mother=females$mother,replicat=females$replicat,treatment=females$treatment,round=females$round), FUN=mean, na.action=na.rm)
FWS
str(FWS)
View(FWS)

MWS <- aggregate(males$WS, by=list(mother=males$mother,replicat=males$replicat,treatment=males$treatment,round=males$round), FUN=mean, na.action=na.rm)
MWS
str(MWS)
View(MWS)

#Analysis
modelF.ws<- lmer(x~treatment+round+(1|treatment:replicat), data=FWS)
summary(modelF.ws)
anova(modelF.ws)
emmeans(modelF.ws,list(pairwise ~treatment),adjust="tukey")
plot(modelF.ws) 
qqnorm(residuals(modelF.ws)) 
qqline(residuals(modelF.ws)) 

modelM.ws<- lmer(x~treatment+round+(1|treatment:replicat), data=MWS)
summary(modelM.ws)
anova(modelM.ws)
emmeans(modelM.ws,list(pairwise ~treatment),adjust="tukey")
plot(modelM.ws) #Looks OK, a few outliers but variance mostly constant
qqnorm(residuals(modelM.ws)) #Normality better than previous models
qqline(residuals(modelM.ws)) 


#Locomotory activity
#Calculate family total counts
FLA <- aggregate(females$LA_active, by=list(mother=females$mother,replicat=females$replicat,treatment=females$treatment,round=females$round), FUN=sum)
FLA
str(FLA)
View(FLA)

MLA <- aggregate(males$LA_active, by=list(mother=males$mother,replicat=males$replicat,treatment=males$treatment,round=males$round), FUN=sum)
MLA
str(MLA)
View(MLA)

#Analysis
modelF.la<- glmer(x~treatment+round+(1|treatment:replicat), data=FLA, family=poisson)
summary(modelF.la)
Anova(modelF.la)
emmeans(modelF.la,list(pairwise ~treatment),adjust="tukey")
simulationOutput <- simulateResiduals(fittedModel=modelF.la, n=250)
plot(simulationOutput)
testDispersion(simulationOutput)

modelM.la<- glmer(x~treatment+round+(1|treatment:replicat), data=MLA, family=poisson)
summary(modelM.la)
Anova(modelM.la)
emmeans(modelM.la,list(pairwise ~treatment),adjust="tukey")
simulationOutput <- simulateResiduals(fittedModel=modelM.la, n=250)
plot(simulationOutput)
testDispersion(simulationOutput)

#Resampling to get unbiased p-values
set.seed(2112)
str(FLA)
nrow(FLA)
N <- 999
CStreat.f <- rep(NA, 999)
CSround.f <- rep(NA, 999)
for (i in 1:N) {
  newdataF <- cbind(FLA[,1:4],x=FLA[sample(1:nrow(FLA), nrow(FLA)),5])
  modelF.nd<- glmer(x~treatment+round+(1|treatment:replicat), family=poisson, data=newdataF)
  out <- Anova(modelF.nd, type=3)
  CStreat.f[i] <- out[2,1]
  CSround.f[i] <- out[3,1]
}
CStreat.f <- sort(CStreat.f)
hist(CStreat.f)
CSf.actual <- Anova(modelF.la, type=3) #0.9451 for treatment, 365.4120 for round
(1000-(sum(CStreat.f <= CSf.actual[2,1])))/1000 #p=0.800
CSround.f <- sort(CSround.f)
hist(CSround.f)
(1000-(sum(CSround.f <= CSf.actual[3,1])))/1000 #p<0.001

str(MLA)
nrow(MLA)
N <- 999
CStreat.m <- rep(NA, 999)
CSround.m <- rep(NA, 999)
for (i in 1:N) {
  newdataM <- cbind(MLA[,1:4],x=MLA[sample(1:nrow(MLA), nrow(MLA)),5])
  modelM.nd<- glmer(x~treatment+round+(1|treatment:replicat), family=poisson, data=newdataM)
  out <- Anova(modelM.nd, type=3)
  CStreat.m[i] <- out[2,1]
  CSround.m[i] <- out[3,1]
}
CStreat.m <- sort(CStreat.m)
hist(CStreat.m)
CSm.actual <- Anova(modelM.la, type=3) #0.8146 for treatment, 59.8034 for round
(1000-(sum(CStreat.m <= CSm.actual[2,1])))/1000 #p=0.831
CSround.m <- sort(CSround.m)
hist(CSround.m)
(1000-(sum(CSround.m <= CSm.actual[3,1])))/1000 #p=0.006

#Visualization
#Figure 2

library(dplyr) 
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(readxl)
library(cowplot)
library(ggplotify)
library(gridGraphics)
theme_set(theme_pubclean())

dev.off()

names(females)[3] <- "Replicate"
names(females)[4] <- "Treatment"
females$Treatment <-factor(females$Treatment, levels = c("FLX","CFM","CWT"))
names(males)[3] <- "Replicate"
names(males)[4] <- "Treatment"
males$Treatment <-factor(males$Treatment, levels = c("FLX","CFM","CWT"))

f.la <- females[,c(3:4,11)]
f.la <- na.omit(f.la)
m.la <- males[,c(3:4,11)]
m.la <- na.omit(m.la)

e <-ggplot(f.la,aes(x = Treatment, y = LA_proportion))
fLA<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Locomotory activity")+ 
  xlab("")+
  labs(title = "Females")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")

theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(fLA,"fLA")

e <-ggplot(m.la,aes(x = Treatment, y = LA_proportion))
mLA<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Locomotory activity")+ 
  xlab("")+
  labs(title = "Males")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")


theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(mLA,"mLA")

f.ws <- females[,c(3:4,7)]
f.ws <- na.omit(f.ws)
m.ws <- males[,c(3:4,7)]
m.ws <- na.omit(m.ws)

e <-ggplot(f.ws,aes(x = Treatment, y = WS))
fWL<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Wing length (mm)")+ 
  xlab("")+
  labs(title = "Females")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")


theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(fWL,"fWL")

e <-ggplot(m.ws,aes(x = Treatment, y = WS))
mWL<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Wing length (mm)")+ 
  xlab("")+
  labs(title = "Males")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")


theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(mWL,"mWL")

f.dr <- females[,c(3:4,8)]
f.dr <- na.omit(f.dr)
m.dr <- males[,c(3:4,8)]
m.dr <- na.omit(m.dr)

e <-ggplot(f.dr,aes(x = Treatment, y = DR))
fDR<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Desiccation resistance (hours)")+ 
  xlab("")+
  labs(title = "Females")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")


theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(fDR,"fDR")

e <-ggplot(m.dr,aes(x = Treatment, y = DR))
mDR<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = "mean",
                      geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = "mean_se",
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Desiccation resistance (hours)")+ 
  xlab("")+
  labs(title = "Males")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.position = "none")


theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  theme(legend.justification = c("right", "top"), legend.key.width=unit(1.5, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="black", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

saveRDS(mDR,"mDR")

fLA
mLA 
fWL
mWL 
fDR
mDR

fig2<-ggdraw() +
  draw_plot(fLA, x = 0.12, y = 0.64, width = 0.35, height = .33) +
  draw_plot(mLA, x = 0.6, y = 0.64, width = 0.35, height = .33)+
  draw_plot(fWL, x = 0.11, y = 0.33, width = 0.35, height = .33)+
  draw_plot(mWL, x = 0.6, y = 0.33, width = 0.35, height = .33)+
  draw_plot(fDR, x = 0.12, y = 0, width = 0.35, height = .33)+
  draw_plot(mDR, x = 0.6, y = 0, width = 0.35, height = .33)
fig2
saveRDS(fig2,"noLegFig2")


# To add the legend 
library("grid") 
library("gridExtra")

levels(m.la$Replicate) <- c("Rep.1", "Rep.2")

e <-ggplot(m.la,aes(x = Treatment, y = LA_proportion ))
ForL<-e + stat_summary(aes(color = Replicate,shape=Replicate),fun.y = mean,
                       geom = "point",size = 2,position = position_dodge(0.2)) +
  stat_summary(
    aes(color = Replicate),fun.data = mean_se,
    geom = "errorbar",linewidth = 0.6,width = 0.1,
    position = position_dodge(0.2)
  )+
  theme_classic()+ 
  ylab("Locomotory activity")+ 
  xlab("")+
  labs(title = "Males")+
  theme(plot.title = element_text(color = "black", size = 14,))+
  theme(axis.text.x = element_text(face="plain", color="black",size=14))+
  theme(axis.text.y = element_text(face="plain", color="black",size=12))+
  theme(axis.title = element_text(face="plain", color="black",size=14))+ 
  theme(legend.text = element_text(size = 13))+
  theme(legend.title = element_blank())+
  theme(legend.position = c(.83, .89), legend.key.width=unit(1.8, "lines"),
        legend.box.margin=margin(c(0.2,0.2,0.2,0.05)))+
  theme(legend.box.background = element_rect(color="white", size=0.5),
        legend.box.margin = margin(1, 1, 1, 0.8))

ggp_legend <- get_legend(ForL)

ggdraw() +
  draw_plot(fig2, x = 0.01, y = 0.01, width = 0.8, height = .95)
leg<-grid.draw(ggp_legend) 
