
# Figure 3. Visualization of the 95% confidence ellipses, ordered in par plotting order

library(MCMCglmm)
library(plotrix)

### Plotting confidence ellipses by hand following visiondummy #####
#https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

#Extract all full 6x6 matrices ####
setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
model<-readRDS("rep1flx")
pm <- posterior.mode(model$VCV)
G_flx1 <- matrix(pm[1:36], 6,6, byrow=TRUE)

model<-readRDS("rep1cfm")
pm <- posterior.mode(model$VCV)
G_cfm1 <- matrix(pm[1:36], 6,6, byrow=TRUE)

model<-readRDS("rep1cwt")
pm <- posterior.mode(model$VCV)
G_cwt1 <- matrix(pm[1:36], 6,6, byrow=TRUE)

model<-readRDS("rep2flx")
pm <- posterior.mode(model$VCV)
G_flx2 <- matrix(pm[1:36], 6,6, byrow=TRUE)

model<-readRDS("rep2cfm")
pm <- posterior.mode(model$VCV)
G_cfm2 <- matrix(pm[1:36], 6,6, byrow=TRUE)

model<-readRDS("rep2cwt")
pm <- posterior.mode(model$VCV)
G_cwt2 <- matrix(pm[1:36], 6,6, byrow=TRUE)

#Create a 6 x 6 plot using par, filled up columnwise ####
par(mfcol = c(6, 6), mar=c(1,1,1,1), oma=c(3,3,3,3), mgp=c(3, 0.5, 0))
#par(mfcol = c(1, 1)) to reset plotting multiple graphs

#plot 1.1 ####
plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "WL m", cex=2)

G_flx1_1 <- matrix(G_flx1[c(8,2,7,1)], 2,2, byrow=TRUE) 
G_cfm1_1 <- matrix(G_cfm1[c(8,2,7,1)], 2,2, byrow=TRUE) 
G_cwt1_1 <- matrix(G_cwt1[c(8,2,7,1)], 2,2, byrow=TRUE) 

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_1) 
ev2<-eigen(G_cfm1_1)
ev3<-eigen(G_cwt1_1)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5),c(-3.5,3.5), type="n", ylab="", xlab="",xaxt = "n",cex.axis=1.5)
  draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
               deg=F, border='red', lty=5, lwd=2) # '#0DA7C9'
  draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
               deg=F, border='#d2ad03', lty=3, lwd=2) 
  draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
               deg=F, border='#009966', lty=1, lwd=2)

G_flx1_2 <- matrix(G_flx1[c(15,3,13,1)], 2,2, byrow=TRUE) 
G_cfm1_2 <- matrix(G_cfm1[c(15,3,13,1)], 2,2, byrow=TRUE) 
G_cwt1_2 <- matrix(G_cwt1[c(15,3,13,1)], 2,2, byrow=TRUE) 

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_2) 
ev2<-eigen(G_cfm1_2)
ev3<-eigen(G_cwt1_2)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="",xaxt = "n",cex.axis=1.5)
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_3 <- matrix(G_flx1[c(22,4,19,1)], 2,2, byrow=TRUE) 
G_cfm1_3 <- matrix(G_cfm1[c(22,4,19,1)], 2,2, byrow=TRUE) 
G_cwt1_3 <- matrix(G_cwt1[c(22,4,19,1)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_3) 
ev2<-eigen(G_cfm1_3)
ev3<-eigen(G_cwt1_3)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="",xaxt = "n",cex.axis=1.5)
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)


G_flx1_4 <- matrix(G_flx1[c(29,5,25,1)], 2,2, byrow=TRUE) 
G_cfm1_4 <- matrix(G_cfm1[c(29,5,25,1)], 2,2, byrow=TRUE) 
G_cwt1_4 <- matrix(G_cwt1[c(29,5,25,1)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_4)
ev2<-eigen(G_cfm1_4)
ev3<-eigen(G_cwt1_4)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="",xaxt = "n",cex.axis=1.5) #main="1.5: ws.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_5 <- matrix(G_flx1[c(36,6,31,1)], 2,2, byrow=TRUE) 
G_cfm1_5 <- matrix(G_cfm1[c(36,6,31,1)], 2,2, byrow=TRUE) 
G_cwt1_5 <- matrix(G_cwt1[c(36,6,31,1)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_5)
ev2<-eigen(G_cfm1_5)
ev3<-eigen(G_cwt1_5)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="",cex.axis=1.5) #, main="1.6: ws.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_1 <- matrix(G_flx2[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_cfm1_1 <- matrix(G_cfm2[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_cwt1_1 <- matrix(G_cwt2[c(1,2,7,8)], 2,2, byrow=TRUE) 

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_1) 
ev2<-eigen(G_cfm1_1)
ev3<-eigen(G_cwt1_1)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="1.2: ws.m-dr.m"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)


#plot 2.2 ####
plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "DR m", cex=2)

G_flx1_6 <- matrix(G_flx1[c(15,9,14,8)], 2,2, byrow=TRUE) 
G_cfm1_6 <- matrix(G_cfm1[c(15,9,14,8)], 2,2, byrow=TRUE) 
G_cwt1_6 <- matrix(G_cwt1[c(15,9,14,8)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_6)
ev2<-eigen(G_cfm1_6)
ev3<-eigen(G_cwt1_6)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.3: dr.m-act.m"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_7 <- matrix(G_flx1[c(22,10,20,8)], 2,2, byrow=TRUE) 
G_cfm1_7 <- matrix(G_cfm1[c(22,10,20,8)], 2,2, byrow=TRUE) 
G_cwt1_7 <- matrix(G_cwt1[c(22,10,20,8)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_7)
ev2<-eigen(G_cfm1_7)
ev3<-eigen(G_cwt1_7)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.4: dr.m-ws.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_8 <- matrix(G_flx1[c(29,11,26,8)], 2,2, byrow=TRUE) 
G_cfm1_8 <- matrix(G_cfm1[c(29,11,26,8)], 2,2, byrow=TRUE) 
G_cwt1_8 <- matrix(G_cwt1[c(29,11,26,8)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_8)
ev2<-eigen(G_cfm1_8)
ev3<-eigen(G_cwt1_8)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.5: dr.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_9 <- matrix(G_flx1[c(36,12,32,8)], 2,2, byrow=TRUE) 
G_cfm1_9 <- matrix(G_cfm1[c(36,12,32,8)], 2,2, byrow=TRUE) 
G_cwt1_9 <- matrix(G_cwt1[c(36,12,32,8)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_9) 
ev2<-eigen(G_cfm1_9)
ev3<-eigen(G_cwt1_9)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", yaxt = "n",cex.axis=1.5) #, main="2.6: dr.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_2 <- matrix(G_flx2[c(1,3,13,15)], 2,2, byrow=TRUE) 
G_cfm1_2 <- matrix(G_cfm2[c(1,3,13,15)], 2,2, byrow=TRUE) 
G_cwt1_2 <- matrix(G_cwt2[c(1,3,13,15)], 2,2, byrow=TRUE) 

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_2) 
ev2<-eigen(G_cfm1_2)
ev3<-eigen(G_cwt1_2)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="1.3: ws.m-act.m"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_6 <- matrix(G_flx2[c(8,9,14,15)], 2,2, byrow=TRUE) 
G_cfm1_6 <- matrix(G_cfm2[c(8,9,14,15)], 2,2, byrow=TRUE) 
G_cwt1_6 <- matrix(G_cwt2[c(8,9,14,15)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_6) 
ev2<-eigen(G_cfm1_6)
ev3<-eigen(G_cwt1_6)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.3: dr.m-act.m"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

#plot 3.3 ####
plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "LA m", cex=2)

G_flx1_10 <- matrix(G_flx1[c(22,16,21,15)], 2,2, byrow=TRUE) 
G_cfm1_10 <- matrix(G_cfm1[c(22,16,21,15)], 2,2, byrow=TRUE) 
G_cwt1_10 <- matrix(G_cwt1[c(22,16,21,15)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_10) 
ev2<-eigen(G_cfm1_10)
ev3<-eigen(G_cwt1_10)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="3.4: act.m-ws.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_11 <- matrix(G_flx1[c(29,17,27,15)], 2,2, byrow=TRUE) 
G_cfm1_11 <- matrix(G_cfm1[c(29,17,27,15)], 2,2, byrow=TRUE) 
G_cwt1_11 <- matrix(G_cwt1[c(29,17,27,15)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_11) 
ev2<-eigen(G_cfm1_11)
ev3<-eigen(G_cwt1_11)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="3.5: act.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_12 <- matrix(G_flx1[c(36,18,33,15)], 2,2, byrow=TRUE) 
G_cfm1_12 <- matrix(G_cfm1[c(36,18,33,15)], 2,2, byrow=TRUE) 
G_cwt1_12 <- matrix(G_cwt1[c(36,18,33,15)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_12)
ev2<-eigen(G_cfm1_12)
ev3<-eigen(G_cwt1_12)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", yaxt = "n",cex.axis=1.5) #, main="3.6: act.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_3 <- matrix(G_flx2[c(1,4,19,22)], 2,2, byrow=TRUE) 
G_cfm1_3 <- matrix(G_cfm2[c(1,4,19,22)], 2,2, byrow=TRUE) 
G_cwt1_3 <- matrix(G_cwt2[c(1,4,19,22)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_3)
ev2<-eigen(G_cfm1_3)
ev3<-eigen(G_cwt1_3)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="1.4: ws.m-ws.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_7 <- matrix(G_flx2[c(8,10,20,22)], 2,2, byrow=TRUE) 
G_cfm1_7 <- matrix(G_cfm2[c(8,10,20,22)], 2,2, byrow=TRUE) 
G_cwt1_7 <- matrix(G_cwt2[c(8,10,20,22)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_7)
ev2<-eigen(G_cfm1_7)
ev3<-eigen(G_cwt1_7)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.4: dr.m-ws.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_10 <- matrix(G_flx2[c(15,16,21,22)], 2,2, byrow=TRUE) 
G_cfm1_10 <- matrix(G_cfm2[c(15,16,21,22)], 2,2, byrow=TRUE) 
G_cwt1_10 <- matrix(G_cwt2[c(15,16,21,22)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_10)
ev2<-eigen(G_cfm1_10)
ev3<-eigen(G_cwt1_10)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="3.4: act.m-ws.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

#plot 4.4 ####
plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "WL f", cex=2)

G_flx1_13 <- matrix(G_flx1[c(29,23,28,22)], 2,2, byrow=TRUE) 
G_cfm1_13 <- matrix(G_cfm1[c(29,23,28,22)], 2,2, byrow=TRUE) 
G_cwt1_13 <- matrix(G_cwt1[c(29,23,28,22)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_13)
ev2<-eigen(G_cfm1_13)
ev3<-eigen(G_cwt1_13)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="4.5: ws.f-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_14 <- matrix(G_flx1[c(36,24,34,22)], 2,2, byrow=TRUE) 
G_cfm1_14 <- matrix(G_cfm1[c(36,24,34,22)], 2,2, byrow=TRUE) 
G_cwt1_14 <- matrix(G_cwt1[c(36,24,34,22)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_14)
ev2<-eigen(G_cfm1_14)
ev3<-eigen(G_cwt1_14)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", yaxt = "n",cex.axis=1.5) #, main="4.6: ws.f-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_4 <- matrix(G_flx2[c(1,5,25,29)], 2,2, byrow=TRUE) 
G_cfm1_4 <- matrix(G_cfm2[c(1,5,25,29)], 2,2, byrow=TRUE) 
G_cwt1_4 <- matrix(G_cwt2[c(1,5,25,29)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_4) 
ev2<-eigen(G_cfm1_4)
ev3<-eigen(G_cwt1_4)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="1.5: ws.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_8 <- matrix(G_flx2[c(8,11,26,29)], 2,2, byrow=TRUE) 
G_cfm1_8 <- matrix(G_cfm2[c(8,11,26,29)], 2,2, byrow=TRUE) 
G_cwt1_8 <- matrix(G_cwt2[c(8,11,26,29)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_8) 
ev2<-eigen(G_cfm1_8)
ev3<-eigen(G_cwt1_8)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.5: dr.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_11 <- matrix(G_flx2[c(15,17,27,29)], 2,2, byrow=TRUE) 
G_cfm1_11 <- matrix(G_cfm2[c(15,17,27,29)], 2,2, byrow=TRUE) 
G_cwt1_11 <- matrix(G_cwt2[c(15,17,27,29)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_11) 
ev2<-eigen(G_cfm1_11)
ev3<-eigen(G_cwt1_11)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="3.5: act.m-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

G_flx1_13 <- matrix(G_flx2[c(22,23,28,29)], 2,2, byrow=TRUE) 
G_cfm1_13 <- matrix(G_cfm2[c(22,23,28,29)], 2,2, byrow=TRUE) 
G_cwt1_13 <- matrix(G_cwt2[c(22,23,28,29)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_13) 
ev2<-eigen(G_cfm1_13)
ev3<-eigen(G_cwt1_13)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="4.5: ws.f-dr.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)


#plot 5.5 ####
plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "DR f", cex=2)

G_flx1_15 <- matrix(G_flx1[c(36,30,35,29)], 2,2, byrow=TRUE) 
G_cfm1_15 <- matrix(G_cfm1[c(36,30,35,29)], 2,2, byrow=TRUE) 
G_cwt1_15 <- matrix(G_cwt1[c(36,30,35,29)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_15) 
ev2<-eigen(G_cfm1_15)
ev3<-eigen(G_cwt1_15)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", yaxt = "n",cex.axis=1.5) #, main="5.6: dr.f-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)


#Plot 6.1: ws.m-act.f rep2####
#G_flx1_1 <- matrix(pm[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_flx1_5 <- matrix(G_flx2[c(1,6,31,36)], 2,2, byrow=TRUE) 
G_cfm1_5 <- matrix(G_cfm2[c(1,6,31,36)], 2,2, byrow=TRUE) 
G_cwt1_5 <- matrix(G_cwt2[c(1,6,31,36)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_5) 
ev2<-eigen(G_cfm1_5)
ev3<-eigen(G_cwt1_5)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="1.6: ws.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)


#Plot 6.2: dr.m-act.f rep2####
#G_flx1_1 <- matrix(pm[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_flx1_9 <- matrix(G_flx2[c(8,12,32,36)], 2,2, byrow=TRUE) 
G_cfm1_9 <- matrix(G_cfm2[c(8,12,32,36)], 2,2, byrow=TRUE) 
G_cwt1_9 <- matrix(G_cwt2[c(8,12,32,36)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_9) 
ev2<-eigen(G_cfm1_9)
ev3<-eigen(G_cwt1_9)

#new covariance
cov<-(-0.92*sqrt(G_flx1_9[1,1]*G_flx1_9[2,2]))
G_flx1_9new<-matrix(c(G_flx2[8],cov,cov,G_flx2[36]),2,2, byrow=TRUE)
ev1<-eigen(G_flx1_9new)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="2.6: dr.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

#Plot 6.3: act.m-act.f rep2####
#G_flx1_1 <- matrix(pm[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_flx1_12 <- matrix(G_flx2[c(15,18,33,36)], 2,2, byrow=TRUE) 
G_cfm1_12 <- matrix(G_cfm2[c(15,18,33,36)], 2,2, byrow=TRUE) 
G_cwt1_12 <- matrix(G_cwt2[c(15,18,33,36)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_12) 
ev2<-eigen(G_cfm1_12)
ev3<-eigen(G_cwt1_12)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="3.6: act.m-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

#G_flx1_1 <- matrix(pm[c(1,2,7,8)], 2,2, byrow=TRUE) 
G_flx1_14 <- matrix(G_flx2[c(22,24,34,36)], 2,2, byrow=TRUE) 
G_cfm1_14 <- matrix(G_cfm2[c(22,24,34,36)], 2,2, byrow=TRUE) 
G_cwt1_14 <- matrix(G_cwt2[c(22,24,34,36)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_14) 
ev2<-eigen(G_cfm1_14)
ev3<-eigen(G_cwt1_14)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="4.6: ws.f-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

#Plot 6.5: dr.f-act.f rep2####
#G_flx1_1 <- matrix(pm[c(1,2,7,8)], 2,2, byrow=TRUE)
G_flx1_15 <- matrix(G_flx2[c(29,30,35,36)], 2,2, byrow=TRUE) 
G_cfm1_15 <- matrix(G_cfm2[c(29,30,35,36)], 2,2, byrow=TRUE) 
G_cwt1_15 <- matrix(G_cwt2[c(29,30,35,36)], 2,2, byrow=TRUE)

#The 95% confidence ellipse can be defined similarly to the axis-aligned case, with the major axis of length 2*sqrt(5.991*lambda_1) and the minor axis of length 2*sqrt(5.991*lambda_2), where lambda_1 and lambda_2 represent the eigenvalues of the covariance matrix.
#the angle of the largest eigenvector towards the x-axis: alfa = arctan(v1(y)/v1(x)), where where v1 is the eigenvector of the covariance matrix that corresponds to the largest eigenvalue

ev1<-eigen(G_flx1_15) 
ev2<-eigen(G_cfm1_15)
ev3<-eigen(G_cwt1_15)

#G_flx1_15[1,2]/sqrt(G_flx1_15[1,1]*G_flx1_15[2,2])
#rho=-1.122223, so the matrix is indefinit

#new covariance
cov<-(-0.96*sqrt(G_flx1_15[1,1]*G_flx1_15[2,2]))
G_flx1_15new<-matrix(c(G_flx2[29],cov,cov,G_flx2[36]),2,2, byrow=TRUE)
ev1<-eigen(G_flx1_15new)

majorlength1<-2*sqrt(5.991*ev1$values[1])
minorlength1<-2*sqrt(5.991*ev1$values[2])
angle1<-atan((ev1$vectors[1,1])/(ev1$vectors[2,1]))

majorlength2<-2*sqrt(5.991*ev2$values[1])
minorlength2<-2*sqrt(5.991*ev2$values[2])
angle2<-atan((ev2$vectors[1,1])/(ev2$vectors[2,1]))

majorlength3<-2*sqrt(5.991*ev3$values[1])
minorlength3<-2*sqrt(5.991*ev3$values[2])
angle3<-atan((ev3$vectors[1,1])/(ev3$vectors[2,1]))

plot(c(-3.5,3.5), c(-3.5,3.5), type="n",ylab="", xlab="", xaxt = "n", yaxt = "n") #, main="5.6: dr.f-act.f"
draw.ellipse(c(0), c(0), majorlength1, minorlength1, angle=angle1, 
             deg=F, border='red', lty=5, lwd=2) #Angle is in radians
draw.ellipse(c(0), c(0), majorlength2, minorlength2, angle=angle2,
             deg=F, border='#d2ad03', lty=3, lwd=2)  #
draw.ellipse(c(0), c(0), majorlength3, minorlength3, angle=angle3, 
             deg=F, border='#009966', lty=1, lwd=2)

plot(c(-3.5,3.5), c(-3.5,3.5), type="n", ylab="", xlab="", xaxt = "n", yaxt = "n")
text(x=0, y=0, "LA f", cex=2)

##############################################################################

#Figure S2, Va density plots for rep 1

library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)

par(mfrow = c(3, 2)) 

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
flx1<-readRDS("rep1flx")
cfm1<-readRDS("rep1cfm")
cwt1<-readRDS("rep1cwt")

# Va of WL in males
mWS.Va.flx1<-flx1$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.flx1<-as.data.frame(mWS.Va.flx1)
mWS.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

mWS.Va.cfm1<-cfm1$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.cfm1<-as.data.frame(mWS.Va.cfm1)
mWS.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

mWS.Va.cwt1<-cwt1$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.cwt1<-as.data.frame(mWS.Va.cwt1)
mWS.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

mWS<-rbind(mWS.Va.flx1,mWS.Va.cfm1,mWS.Va.cwt1)
mWS$Treatment<-as.factor(mWS$Treatment)
mWS$Treatment<-factor(mWS$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mWS %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Mwing1<-ggplot(mWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of DR in males
mDR.Va.flx1<-flx1$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.flx1<-as.data.frame(mDR.Va.flx1)
mDR.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

mDR.Va.cfm1<-cfm1$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.cfm1<-as.data.frame(mDR.Va.cfm1)
mDR.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

mDR.Va.cwt1<-cwt1$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.cwt1<-as.data.frame(mDR.Va.cwt1)
mDR.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

mDR<-rbind(mDR.Va.flx1,mDR.Va.cfm1,mDR.Va.cwt1)
mDR$Treatment<-as.factor(mDR$Treatment)
mDR$Treatment<-factor(mDR$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mDR %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

MdesR1<-ggplot(mDR, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of LA in males
mLA.Va.flx1<-flx1$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.flx1<-as.data.frame(mLA.Va.flx1)
mLA.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

mLA.Va.cfm1<-cfm1$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.cfm1<-as.data.frame(mLA.Va.cfm1)
mLA.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

mLA.Va.cwt1<-cwt1$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.cwt1<-as.data.frame(mLA.Va.cwt1)
mLA.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

mLA<-rbind(mLA.Va.flx1,mLA.Va.cfm1,mLA.Va.cwt1)
mLA$Treatment<-as.factor(mLA$Treatment)
mLA$Treatment<-factor(mLA$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mLA %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Mloco1<-ggplot(mLA, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of WL in females
fWS.Va.flx1<-flx1$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.flx1<-as.data.frame(fWS.Va.flx1)
fWS.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

fWS.Va.cfm1<-cfm1$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.cfm1<-as.data.frame(fWS.Va.cfm1)
fWS.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

fWS.Va.cwt1<-cwt1$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.cwt1<-as.data.frame(fWS.Va.cwt1)
fWS.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

fWS<-rbind(fWS.Va.flx1,fWS.Va.cfm1,fWS.Va.cwt1)
fWS$Treatment<-as.factor(fWS$Treatment)
fWS$Treatment<-factor(fWS$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fWS %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Fwing1<-ggplot(fWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of DR in females
fDR.Va.flx1<-flx1$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.flx1<-as.data.frame(fDR.Va.flx1)
fDR.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

fDR.Va.cfm1<-cfm1$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.cfm1<-as.data.frame(fDR.Va.cfm1)
fDR.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

fDR.Va.cwt1<-cwt1$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.cwt1<-as.data.frame(fDR.Va.cwt1)
fDR.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

fDR<-rbind(fDR.Va.flx1,fDR.Va.cfm1,fDR.Va.cwt1)
fDR$Treatment<-as.factor(fDR$Treatment)
fDR$Treatment<-factor(fDR$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fDR %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

FdesR1<-ggplot(fDR, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of LA in females
fLA.Va.flx1<-flx1$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.flx1<-as.data.frame(fLA.Va.flx1)
fLA.Va.flx1$Treatment <- matrix("FLX", nrow=1000)

fLA.Va.cfm1<-cfm1$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.cfm1<-as.data.frame(fLA.Va.cfm1)
fLA.Va.cfm1$Treatment <- matrix("CFM", nrow=1000)

fLA.Va.cwt1<-cwt1$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.cwt1<-as.data.frame(fLA.Va.cwt1)
fLA.Va.cwt1$Treatment <- matrix("CWT", nrow=1000)

fLA<-rbind(fLA.Va.flx1,fLA.Va.cfm1,fLA.Va.cwt1)
fLA$Treatment<-as.factor(fLA$Treatment)
fLA$Treatment<-factor(fLA$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fLA %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Floco1<-ggplot(fLA, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

Mwing1
MdesR1
Mloco1

Fwing1
FdesR1
Floco1

density1<-ggdraw() +
  draw_plot(Mwing1, x = 0.12, y = 0.66, width = 0.35, height = .33) +
  draw_plot(Fwing1, x = 0.6, y = 0.66, width = 0.35, height = .33)+
  draw_plot(MdesR1, x = 0.11, y = 0.33, width = 0.35, height = .33)+
  draw_plot(FdesR1, x = 0.6, y = 0.33, width = 0.35, height = .33)+
  draw_plot(Mloco1, x = 0.12, y = 0, width = 0.35, height = .33)+
  draw_plot(Floco1, x = 0.6, y = 0, width = 0.35, height = .33)
density1

#Add legend
#Create extra plot to extract legend information
FwingX<-ggplot(fWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = c(.83, .89))
gg_legend <- get_legend(FwingX)

ggdraw() +
  draw_plot(density1, x = 0.01, y = 0.01, width = 0.8, height = .95)
leg<-grid.draw(gg_legend)

##############################################################################

#Figure S3, Va density plots for rep 2

setwd("E:\\Lund\\Flies\\Yesbol\\Quant gen manuscript\\multivariate models")
flx2<-readRDS("rep2flx")
cfm2<-readRDS("rep2cfm")
cwt2<-readRDS("rep2cwt")

# Va of WL in males
mWS.Va.flx2<-flx2$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.flx2<-as.data.frame(mWS.Va.flx2)
mWS.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

mWS.Va.cfm2<-cfm2$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.cfm2<-as.data.frame(mWS.Va.cfm2)
mWS.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

mWS.Va.cwt2<-cwt2$VCV[,"traitmeanWs.m:traitmeanWs.m.animal"]
mWS.Va.cwt2<-as.data.frame(mWS.Va.cwt2)
mWS.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

mWS<-rbind(mWS.Va.flx2,mWS.Va.cfm2,mWS.Va.cwt2)
mWS$Treatment<-as.factor(mWS$Treatment)
mWS$Treatment<-factor(mWS$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mWS %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Mwing2<-ggplot(mWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of DR in males
mDR.Va.flx2<-flx2$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.flx2<-as.data.frame(mDR.Va.flx2)
mDR.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

mDR.Va.cfm2<-cfm2$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.cfm2<-as.data.frame(mDR.Va.cfm2)
mDR.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

mDR.Va.cwt2<-cwt2$VCV[,"traitmeanDr.m:traitmeanDr.m.animal"]
mDR.Va.cwt2<-as.data.frame(mDR.Va.cwt2)
mDR.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

mDR<-rbind(mDR.Va.flx2,mDR.Va.cfm2,mDR.Va.cwt2)
mDR$Treatment<-as.factor(mDR$Treatment)
mDR$Treatment<-factor(mDR$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mDR %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

MdesR2<-ggplot(mDR, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of LA in males
mLA.Va.flx2<-flx2$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.flx2<-as.data.frame(mLA.Va.flx2)
mLA.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

mLA.Va.cfm2<-cfm2$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.cfm2<-as.data.frame(mLA.Va.cfm2)
mLA.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

mLA.Va.cwt2<-cwt2$VCV[,"traitmeanActive.m:traitmeanActive.m.animal"]
mLA.Va.cwt2<-as.data.frame(mLA.Va.cwt2)
mLA.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

mLA<-rbind(mLA.Va.flx2,mLA.Va.cfm2,mLA.Va.cwt2)
mLA$Treatment<-as.factor(mLA$Treatment)
mLA$Treatment<-factor(mLA$Treatment,levels = c("FLX","CFM","CWT"))

mu <- mLA %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Mloco2<-ggplot(mLA, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of WL in females
fWS.Va.flx2<-flx2$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.flx2<-as.data.frame(fWS.Va.flx2)
fWS.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

fWS.Va.cfm2<-cfm2$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.cfm2<-as.data.frame(fWS.Va.cfm2)
fWS.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

fWS.Va.cwt2<-cwt2$VCV[,"traitmeanWs.f:traitmeanWs.f.animal"]
fWS.Va.cwt2<-as.data.frame(fWS.Va.cwt2)
fWS.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

fWS<-rbind(fWS.Va.flx2,fWS.Va.cfm2,fWS.Va.cwt2)
fWS$Treatment<-as.factor(fWS$Treatment)
fWS$Treatment<-factor(fWS$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fWS %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Fwing2<-ggplot(fWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of DR in females
fDR.Va.flx2<-flx2$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.flx2<-as.data.frame(fDR.Va.flx2)
fDR.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

fDR.Va.cfm2<-cfm2$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.cfm2<-as.data.frame(fDR.Va.cfm2)
fDR.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

fDR.Va.cwt2<-cwt2$VCV[,"traitmeanDr.f:traitmeanDr.f.animal"]
fDR.Va.cwt2<-as.data.frame(fDR.Va.cwt2)
fDR.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

fDR<-rbind(fDR.Va.flx2,fDR.Va.cfm2,fDR.Va.cwt2)
fDR$Treatment<-as.factor(fDR$Treatment)
fDR$Treatment<-factor(fDR$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fDR %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

FdesR2<-ggplot(fDR, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

# Va of LA in females
fLA.Va.flx2<-flx2$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.flx2<-as.data.frame(fLA.Va.flx2)
fLA.Va.flx2$Treatment <- matrix("FLX", nrow=1000)

fLA.Va.cfm2<-cfm2$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.cfm2<-as.data.frame(fLA.Va.cfm2)
fLA.Va.cfm2$Treatment <- matrix("CFM", nrow=1000)

fLA.Va.cwt2<-cwt2$VCV[,"traitmeanActive.f:traitmeanActive.f.animal"]
fLA.Va.cwt2<-as.data.frame(fLA.Va.cwt2)
fLA.Va.cwt2$Treatment <- matrix("CWT", nrow=1000)

fLA<-rbind(fLA.Va.flx2,fLA.Va.cfm2,fLA.Va.cwt2)
fLA$Treatment<-as.factor(fLA$Treatment)
fLA$Treatment<-factor(fLA$Treatment,levels = c("FLX","CFM","CWT"))

mu <- fLA %>%
  group_by(Treatment) %>%
  summarise(grp.mean = mean(var1))

Floco2<-ggplot(fLA, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = "none")

Mwing2
MdesR2
Mloco2

Fwing2
FdesR2
Floco2

density2<-ggdraw() +
  draw_plot(Mwing2, x = 0.12, y = 0.66, width = 0.35, height = .33) +
  draw_plot(Fwing2, x = 0.6, y = 0.66, width = 0.35, height = .33)+
  draw_plot(MdesR2, x = 0.11, y = 0.33, width = 0.35, height = .33)+
  draw_plot(FdesR2, x = 0.6, y = 0.33, width = 0.35, height = .33)+
  draw_plot(Mloco2, x = 0.12, y = 0, width = 0.35, height = .33)+
  draw_plot(Floco2, x = 0.6, y = 0, width = 0.35, height = .33)
density2

#Add legend
#Create extra plot to extract legend information
FwingX<-ggplot(fWS, aes(x = var1, colour = Treatment)) +
  geom_density(aes(fill = Treatment, color = Treatment), alpha = 0.3)+
  geom_vline(aes(xintercept = grp.mean, color = Treatment),linewidth=0.8,data = mu,linetype = 2) +
  labs(x="Additive genetic variance (Va)", y = "Density")+
  theme_classic()+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 14))+
  theme(legend.position = c(.83, .89))
gg_legend <- get_legend(FwingX)

ggdraw() +
  draw_plot(density2, x = 0.01, y = 0.01, width = 0.8, height = .95)
leg<-grid.draw(gg_legend)