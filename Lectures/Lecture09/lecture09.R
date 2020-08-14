## ----echo=FALSE,results='hide',message=FALSE----------------------------------
######################################################################
### Author: Michael Höhle <http://www.math.su.se/~hoehle>
### Course: The Mathematics and Statistics of Infectious Disease Outbreaks
###         (MT3002 - Summer 2020) at the Department of Mathematics,
###         Stockholm University.
###         (https://kurser.math.su.se/course/view.php?id=911)
###         given by Tom Britton and Michael Höhle
###    
### License (for the slides):
### This work is licensed under a <a rel="license"
### href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
### Attribution-ShareAlike 4.0 International License</a>.
### License (for the code):
### GNU General Public License v3 - https://www.gnu.org/licenses/gpl-3.0.en.html
###
### Description:
###  Rnw File for lecture 08
###
### History:
###  -- 04 Aug 2020 file created 
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)
library(lubridate)


## -----------------------------------------------------------------------------
#source("sweave-header.R")
options(width=90,prompt="R> ")
library("surveillance")
library("RColorBrewer")
library("MASS")
library("mgcv")
#Load EuroMOMO data
data("momo")
#Convert to data frame
momo.df <- as.data.frame(momo)



## -----------------------------------------------------------------------------
library("surveillance")
data(salmNewport)
#Load Salmonella Newport data
#load("../Data/sNewport.RData")
#Load extra plotting functionality
source("R/plotOne.R")
sNewport <- aggregate(salmNewport, by = "unit")
#Make a reduced time series for illustration
sNewportAll <- sNewport[epoch(sNewport) >= as.Date("2009-01-01") & epoch(sNewport) <= as.Date("2011-12-31"),]

#From where to start
start <- as.Date("2011-09-26")+7
isoWeekYear(start)



## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 0*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 1*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 2*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 3*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 4*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 5*7,doSurv=FALSE,ylim=c(0,20))


## -----------------------------------------------------------------------------
plotOne(sts=sNewportAll,start=start,now=start + 5*7,doSurv=TRUE,ylim=c(0,20))
legend(x="topleft",expression(q[0.95]* " of predictive distribution"),col="red",lwd=2L,lty=1)


## ----echo=FALSE---------------------------------------------------------------
plot(momo[year(momo)>=2000,],type=observed ~ time | unit,par.list=list(mar=c(4,4,1,1)),ylab="No. deaths")


## ----echo=FALSE---------------------------------------------------------------
momo2 <- momo
momo2@observed <- observed(momo)/population(momo) * 100000
plot(momo2[year(momo2)>=2000,],ylab="Deaths per 100.000",type=observed ~ time | unit,par.list=list(mar=c(4,4,1,1)))


## ----PREDINT------------------------------------------------------------------
library(gamlss)
set.seed(123)
n <- 5
mu <- 8
sigma <- 2
y <- rnorm(n=5,mean=mu,sd=sigma)
#Estimates
(ybar <- mean(y))
(sd <- sd(y))

#Calculate the density on a grid in the three cases
z.grid <- seq(mu-3*sigma,mu+3*sigma,length=1000)
d.true      <- dnorm(z.grid, mean=mu, sd=sigma)
d.plugin    <- dnorm(z.grid, mean=ybar, sd=sd)
d.sampleadj <- dTF(z.grid, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1)

#Compute two-sided 95% PI
alpha <- 0.05
pi.true <- mu + qnorm(1-alpha/2) * c(-1,1) * sigma
pi.plugin <- ybar + qnorm(1-alpha/2) * c(-1,1) * sd
pi.sampleadj <- ybar + qt(1-alpha/2,df=n-1) * c(-1,1) * sd * sqrt(1+1/n)
mu.ci <- ybar + c(-1,1) * qnorm(1-alpha) * sd / sqrt(n)

tab <- rbind(pi.true, pi.plugin, pi.sampleadj, mu.ci)
print(tab,digits=2)

#How much more length?
(diff(pi.sampleadj)/diff(pi.plugin)-1)*100

#Compare the last one
pi <- qTF(c(alpha/2,1-alpha/2), mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1)

# Check coverage
diff(pTF(pi,  mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))

#Coverage of different interval types
p_piplugin <- diff(pTF(pi.plugin, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))
p_muci <- diff(pTF(mu.ci, mu=ybar, sigma=sd*sqrt(1+1/n), nu=n-1))

## ----PLOTPREDINT--------------------------------------------------------------
layout(c(1,2),heights=c(5,1.5))
par(mar = c(2,5,1,1))

matplot(z.grid, cbind(d.plugin,d.sampleadj) ,type="l",xlab="z",ylab="density",col=c("royalblue","magenta"),lty=1,lwd=2)
rug(y,lwd=2)
legend(x="topleft",c("normal (plug-in)","non-standard t"),lwd=2,lty=1,col=c("royalblue","magenta"))

par(mar = c(2,5,1,1))

plot.pi <- function(pi,where=0,col="black") {
  lines( c(pi[1],pi[2]),rep(where,2),lty=1,lwd=2,col=col)
  lines( c(pi[1],pi[1]),where + c(-1,1)*0.25,lty=1,lwd=2,col=col)
  lines( c(pi[2],pi[2]),where + c(-1,1)*0.25,lty=1,lwd=2,col=col)
}

plot( c(pi.true[1],pi.true[2]),c(0,0),axes=FALSE,xlim=range(z.grid),xlab="",ylab="",ylim=c(0.75,2.25),col="orange",lty=1,lwd=2,type="n")

#plot.pi(pi.true,where=0,col="orange")
plot.pi(pi.plugin,where=1,col="royalblue")
plot.pi(pi.sampleadj,where=2,col="magenta")
#axis(1)
title(ylab="95% Pred.\n intervals",line=2)


## -----------------------------------------------------------------------------
#Setup monitoring and surveillance periods for the euromomo data
phase2 <- which(epoch(momo) >= "2007-10-01")
phase1 <- which(year(momo) == 2002 & epochInYear(momo) == 40):(phase2[1]-1)

######################################################################
# Illustrate how the Farrington algorithm works
######################################################################
par(mar=c(2, 4, 4, 2) + 0.1)
#Setup Farrington control object
alpha <- 0.005
cntrlFar <- list(range=phase2,alpha=alpha,b=5,w=4)#,powertrans="2/3")
#Adopt control object so exactly one computation is shown and
#setting the argument "plot" a picture is shown
one.t <- cntrlFar ; one.t$range <- min(one.t$range) ; one.t$plot <- TRUE
#Perform surveillance for one time point and save the graphics
onet <- farrington(momo[,"[75,85)"],control=one.t)


## ----echo=FALSE,results='hide', fig.keep="none"-------------------------------
#Load data and convert to an S4 object
data("hepatitisA")
hepatitisA <- disProg2sts(hepatitisA)

#Define parameters to use for Farrington algo
b <- 3
w <- 5
t0 <- 190

#Index of the reference values
t.ref <- 190 - (rep(1:b,each=2*w+1)*52 +  seq(-w,w,by=1))
#Make a data frame with the reference values
refvals <- data.frame(t=t.ref, y=observed(hepatitisA)[t.ref,])

#Show result
with(refvals,plot(t, y))

#Fit quasi-Poisson GLM
m <- glm(y ~ t, data=refvals, family=quasipoisson)
summary(m)

# Use the \texttt{predict} function to compute
#  $\hat{\mu}_{t_0}$ and $\sqrt{\Var(\hat{\mu}_{t_0})}$.
m.pred <- predict(m, newdata=data.frame(t=t0), se=TRUE,type="response")
phi <- summary(m)$dispersion

#Check computations
eta.hat <- coef(m) %*% c(1,t0)
mu.hat <- exp(eta.hat)
c(predict=m.pred$fit, manual=mu.hat)

var.log.mu.hat <- vcov(m)[1,1] + t0^2*vcov(m)[2,2] + 2*t0*vcov(m)[1,2]
c(predict=predict(m, newdata=data.frame(t=t0), se=TRUE,type="link")$se.fit,manual=sqrt(var.log.mu.hat))

var.mu.hat <- exp(eta.hat)^2 * var.log.mu.hat
c(predict.se=m.pred$se.fit,manual.se=sqrt(var.mu.hat))

#Aside: We actually know the true sd of the transformation,
#because result is a logN distribution
#sqrt(exp(var.log.mu.hat)-1)*exp(2*eta.hat + var.log.mu.hat)


#Compute the upper limit of a 99\% prediction interval for $y_{t_0}$
alpha <- 0.05
(pi <- m.pred$fit + c(-1,1)*qnorm(1-alpha/2)*sqrt( phi*m.pred$fit + m.pred$se.fit^2))

#Same result with surveillance package
control <- list(range=190,b=b,w=w,powertrans="none",alpha=alpha,trend=TRUE,reweight=FALSE,fitFun="algo.farrington.fitGLM")
surv <- farrington(hepatitisA, control)

upperbound(surv)


## -----------------------------------------------------------------------------
#Run the farrington algo and the improved version
cntrlFar$limit54 <- c(0,4)
cntrlFar2 <- modifyList(cntrlFar,
                        list(noPeriods=10,
                             populationOffset=FALSE,
                             fitFun="algo.farrington.fitGLM.flexible",
                             weightsThreshold=2.58,
                             pastWeeksNotIncluded=26,
                             pThresholdTrend=1,
                             thresholdMethod="nbPlugin"))

#Old call to farrington -- note that the prediction interval here is two sided
#and hence alpha has to be modified to be comparable
s.far <- farrington(momo[,"[75,85)"],control=modifyList(cntrlFar,list(alpha=alpha*2)))


#Show the results of the surveillance
plot(s.far,dx.upperbound=0,legend.opts=NULL,xlab="time (weeks)",alarm.symbol=list(pch=24,col=1,cex=1),col=c("darkgray",NA,1),lwd=c(1,1,2),lty=c(1,1,1),main="",ylim=c(-35,max(upperbound(s.far))))


#Add legend
legend(x=1,y=130,c("Alarm","NNBA"),pch=c(24,NA),lty=c(NA,1),horiz=TRUE,bg="white",lwd=c(2,3),col=c(1,1))


## ----eval=FALSE,echo=FALSE----------------------------------------------------
## alarmDates <- epoch(s.far[alarms(s.far) == 1,])
## formatDate(alarmDates,"%V")

