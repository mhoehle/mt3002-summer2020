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
###  Rnw File for lecture 06
###
### History:
###  -- 27 Jun 2020 file created 
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)
library(lubridate)


## -----------------------------------------------------------------------------
######################################################################
#Define values to use in example
######################################################################
d.max <- 50
logmu <- 2
logsd <- 0.6

#Sample an outbreak of a point source with n individuals
n <- 55
#Start of the outbreak
t0 <- 25
#Length (if used)
l <- 10
#Size of length outbreak
n.l <- 200


## -----------------------------------------------------------------------------
######################################################################
# PMF obtained by interval-censoring a log-normal
d.dlnorm <- function(d,meanlog=0,sdlog=1) {
  ifelse(d<1,0,plnorm(d,meanlog=meanlog,sdlog=sdlog) - plnorm(d-1,meanlog,sdlog))
}

#Create a grid and norm PMF
d.grid <- seq(0,d.max,by=1)
d.pmf <- d.dlnorm(d.grid,meanlog=logmu,sdlog=logsd)
d.pmf <- d.pmf/sum(d.pmf)
sum(d.pmf)

#Median length
mD <- d.grid[which(cumsum(d.pmf)>0.5)][1]

barplot(height=d.pmf,names.arg=d.grid,width=1,space=0,xlab="Delay D (in days)",ylab="Probability")


## -----------------------------------------------------------------------------
set.seed(123)
#Sample disease onsets from incubation time PMF
D <- sample(d.grid,size=n,prob=d.pmf, replace=TRUE)
t.grid <- seq(1,t0+max(d.grid)+25,by=1)

#Time series
y.ts <- table(factor(t0 + D, levels=t.grid))
plot(t.grid,as.numeric(y.ts),type="h",xlab="Time (days)",ylab="Cases")
lines(rep(t0,2),c(0,1e99),col="red",lty=2,lwd=2)


## -----------------------------------------------------------------------------
set.seed(1234)
#Sample disease onsets from incubation time PMF
D <- sample(d.grid,size=n.l,prob=d.pmf, replace=TRUE)
#Sample when infected in the interval where source is active
t.exposure <- sample(seq(t0,t0+l-1,by=1),size=n.l,replace=TRUE)

#Time series (add extra space to make sure ts is done)
t.grid <- seq(1,max(t.exposure+D+25),by=1)
y.l.ts <- table(factor(t.exposure + D, levels=t.grid))

plot(t.grid,as.numeric(y.l.ts),type="h",xlab="Time (days)",ylab="Cases")
lines(rep(t0,2),c(0,1e99),col="red",lty=2,lwd=2)
lines(rep(t0+l,2),c(0,1e99),col="red",lty=2,lwd=2)

#Exposure times as a time series
n.l.ts <- table(factor(t.exposure, levels=t.grid))
#lines(t.grid,n.l.ts,col="red",lty=2,lwd=2,type="l")


## ----echo=TRUE,results='markup'-----------------------------------------------
subtract.minmax <- function(y, d.pmf,eps=1e-3) {
  exposure.left <- head(which(y>eps),n=1) - ((0:d.max)[head(which(d.pmf>eps),n=1)])
  exposure.right <- tail(which(y>eps),n=1) - ((0:d.max)[tail(which(d.pmf>eps),n=1)])
  structure( c(exposure.left,exposure.right-exposure.left),names=c("t0","l"))
}
subtract.minmax(y.ts, d.pmf)
subtract.minmax(y.l.ts, d.pmf)


## ----echo=TRUE----------------------------------------------------------------
subtract.median <- function(y,d.pmf) {
  d.median <- (0:length(d.pmf)-1)[which(cumsum(d.pmf)>0.5)][1]
  structure(c(tail(y,n=-d.median),rep(0,d.median)),names=names(y))
}
subtract.median(y.ts,d.pmf)


## -----------------------------------------------------------------------------
matplot(cbind(1:length(y.ts),1:length(y.ts)+0.4),cbind(y.ts,subtract.median(y.ts,d.pmf)),type="h",xlab="Time (days)",ylab="Cases",col=c(1,"lightgray"),lty=1,lwd=3)
lines(rep(t0,2),c(0,1e99),col="red",lty=2,lwd=2)
legend(x="topright",c("Onset","Onset-median"),lty=1,lwd=3,col=c(1,"lightgray"))


## ----eval=FALSE---------------------------------------------------------------
## ######################################################################
## # Old figure with axis in other direction
## ######################################################################
## cex.global <- 7
## cex.text <- 1.3
## dx <- 0.1
## 
## ######################################################################
## # Function to draw a node annotated with its position
## ######################################################################
## node <- function(x,y,xstr=deparse(substitute(x)),ystr=deparse(substitute(y)),cex=cex.global,...) {
##   points(x=x,y=y,cex=cex,...)
##   text(x,y,parse(text=paste("n[\"",xstr,",",ystr,"\"]",sep="")),cex=cex.text)
## }
## 
## par(mar=c(1,3,1,1))
## plot(NA,xlim=c(0,6),ylim=c(-1,6),type="n",axes=FALSE,xlab="",ylab="")
## arrows(0,0,6,0,code=2)
## arrows(0.5,-0.5,0.5,6,code=2)
## 
## text(5.85,0.75,"Time t",cex=cex.text)
## text(1,6,"Delay d",cex=cex.text)
## 
## for (t in 1:5) {
##   for (d in (5-t):0) {
##     node(t,d,bg="white",pch=21,xstr=as.character(t),ystr=as.character(d))
##   }
## 
##   for (i in (-1):1) { text(t,5-t+1-i*dx,".")}
## }
## for (i in 1:6) {
##   lines(rep(i+0.5,2),c(-0.5,6),lty=2,col="lightgray")
## }
## for (i in 0:5) {
##   lines( c(0.5,6),rep(i+0.5,2),lty=2,col="lightgray")
## }
## for (i in 1:5) {
##   text(i,-1,substitute(n[i],list(i=i)),cex=cex.text)
## }

## -----------------------------------------------------------------------------
cex.global <- 7
cex.text <- 1.3
dx <- 0.1

######################################################################
# Function to draw a node annotated with its position
######################################################################
node <- function(x,y,xstr=deparse(substitute(x)),ystr=deparse(substitute(y)),cex=cex.global,...) {
  points(x=x,y=y,cex=cex,...)
  text(x,y,parse(text=paste("n[\"",xstr,",",ystr,"\"]",sep="")),cex=cex.text)
  #text(x,y,parse(text=paste("n[\"",ystr,",",xstr,"\"]",sep="")),cex=cex.text)
}

par(mar=c(1,1,1,1))
plot(NA,xlim=c(-1.5,7),ylim=c(0.5,6),type="n",axes=FALSE,xlab="",ylab="")
arrows(-0.5,0.5,6,0.5,code=2) #delay arrow
arrows(0,0,0,6,code=2) #time arrow

text(6,1,"Delay d",cex=cex.text)
text(0.75,6,"Time t",cex=cex.text)

for (t in 1:5) {
  for (d in (5-t):0) {
    node(d,t,bg="white",pch=21,xstr=as.character(t),ystr=as.character(d))
  }
  #3 dots
  text(t,5-t+1,expression(cdots))
}
#Vertical lines
for (i in -1:5) {
  lines(rep(i+0.5,2),c(0.5,5.5),lty=2,col="lightgray")
}
#Horiznal lines
for (i in 0:5) {
  lines( c(-0.5,5.5),rep(i+0.5,2),lty=2,col="lightgray")
}
for (i in 1:5) {
  text(-1,i,substitute(n[i],list(i=i)),cex=cex.text)
}


## ----warning=FALSE------------------------------------------------------------
require("surveillance")
require("RColorBrewer")
######################################################################
#Define values to use in example
######################################################################
d.max <- 50
logmu <- 2
logsd <- 0.6

#Sample an outbreak of a point source with n individuals
n <- 55
#Start of the outbreak
t0 <- 25
#Length (if used)
l <- 10
#Size of length outbreak
n.l <- 200

######################################################################
# PMF obtained by interval-censoring a log-normal
d.dlnorm <- function(d,meanlog=0,sdlog=1) {
  ifelse(d<1,0,plnorm(d,meanlog=meanlog,sdlog=sdlog) - plnorm(d-1,meanlog,sdlog))
}

#Create a grid and norm PMF
d.grid <- seq(0,d.max,by=1)
d.pmf <- d.dlnorm(d.grid,meanlog=logmu,sdlog=logsd)
d.pmf <- d.pmf/sum(d.pmf)
sum(d.pmf)

#Create wrapper functions for the PMF and CDF based on the vector
inc.pmf <- d.pmf
dincu <- function(x) {
  notInSupport <- x<0 | x>=length(inc.pmf)
    #Give index -1 to invalid queries
  x[notInSupport] <- -1
  return(c(0,inc.pmf)[x+2])
}

ylim.global <- NULL

plotConvolution <- function(sts, ...) {
  lambda <- upperbound(sts)
  y <- observed(sts)
  t.grid <- 1:length(y)

  #Compute masses for the individual delays
  M <- matrix(NA,ncol=length(t.grid),nrow=length(t.grid))
  for (t in 1:length(t.grid)) {
    for (d in seq_len(length(t.grid)-t)) {
      M[t,t+d] <- dincu(d)*lambda[t]
    }
  }
  #Plot the stacked convoluted mu_t's
  par(mar=c(4,4,1,1))
  pal <- brewer.pal(n=length(t.grid),"Set3")

  #Set ylim if no other value is given
  if (is.null(ylim.global)) {
    ylim <- c(0,max(apply(M,1,sum,na.rm=TRUE),lambda)+1)
  } else {
    ylim <- ylim.global
  }

  barplot(M[t0:tail(which(y>0),n=1),],col=pal,width=1,space=0,ylim=ylim,xlim=c(t0,tail(which(y>0),n=1)),xlab="Time (days)",ylab="HUS Cases")
  axis(1,at=seq(25,100,by=5)-0.5,label=seq(25,100,by=5))
  #n_t
  points(t.grid[lambda>0]-0.5,lambda[lambda>0],pch=16,col=1,cex=1.2)
  points(t.grid[lambda>0]-0.5,lambda[lambda>0],pch=20,col=pal)
  #y_t
  points(t.grid[y>0]-0.5,y[y>0],pch=4)
  legend(x="topright",c(expression(n[t],y[t])),pch=c(16,4))

  invisible()
}


#Fix seed value for generating the example
set.seed(1234)
#Sample disease incubation times from PMF
D <- sample(d.grid,size=n.l,prob=d.pmf, replace=TRUE)
#Example 2: Sample when infected in the interval where source is active
t.exposure <- sample(seq(t0,t0+l-1,by=1),size=n.l,replace=TRUE)

#Time series (add extra space to make sure ts is done)
t.grid <- seq(1,max(t.exposure+D+25),by=1)
y.l.ts <- table(factor(t.exposure + D, levels=t.grid))

#Exposure times as a time series
n.l.ts <- table(factor(t.exposure, levels=t.grid))

stsl <- new("sts",epoch=1:length(y.l.ts),observed=matrix(y.l.ts,ncol=1))
upperbound(stsl) <- matrix(n.l.ts,ncol=1)
#Make an an illustrative plot
plotConvolution(stsl,ylim=NULL)


## ----eval=FALSE---------------------------------------------------------------
## ######################################################################
## # Show animation
## ######################################################################
## ylim.global <- c(0,30)
## for (pos in 25:28) {
##   n.l.ts.copy <- n.l.ts
##   n.l.ts.copy[pos] <- 0
##   for (i in 0:30) {
##     n.l.ts.copy[pos] <- i
##     upperbound(stsl) <- matrix(n.l.ts.copy,ncol=1)
##    #Make an an illustrative plot
##     plotConvolution(stsl)
##     Sys.sleep(0.1)
##   }
## }
## 


## ----echo=TRUE,cache=TRUE,tidy=FALSE------------------------------------------
#Create vector with incubation time PMF values on (0,...,d_max)
incu.pmf <- c(0, (plnorm(1:d.max,logmu,logsd) - plnorm(0:(d.max-1),logmu,logsd))/plnorm(d.max,logmu,logsd))
#Create sts object
require("surveillance")
sts <- new("sts",epoch=1:length(y.ts),observed=matrix(y.ts,ncol=1))
#Backproject using the method by Becker et al. (1991)
bp.control <- list(k=0,eps=1e-3,iter.max=100,verbose=TRUE,eq3a.method="C")
sts.bp.k0 <- backprojNP(sts, incu.pmf=incu.pmf, control=bp.control)

## ----cache=TRUE---------------------------------------------------------------
sts.bp.k2 <- backprojNP(sts, incu.pmf=incu.pmf, modifyList(bp.control, list(k=2)))
sts.bp.k4 <- backprojNP(sts, incu.pmf=incu.pmf, modifyList(bp.control, list(k=4)))


## ----PLOT, eval=FALSE, echo=TRUE----------------------------------------------
## plot(sts.bp.k0,xaxis.labelFormat=NULL)


## ----PLOTIT3, echo=FALSE, fig.height=3----------------------------------------
par(mar=c(4.2,4,1,1))
plot(sts.bp.k0,xaxis.labelFormat=NULL,xlab="Time (days)",legend=NULL,dx.upperbound=0,,lty=c(1,1,1),lwd=c(1,1,3),main="")
legend(x="topleft",c(expression(Y[t]),expression(hat(lambda)[t])),col=c(1,"blue"),lwd=c(1,3),lty=c(1,1))


## -----------------------------------------------------------------------------
#Plot results for example 2
plotIt1 <- function(which) {
  plot(res[[which]],xaxis.labelFormat=NULL,xlab="Time (days)",legend=NULL,dx.upperbound=0,lty=c(1,1,1),lwd=c(1,1,3),main=paste("k=",res[[which]]@control$k,sep=""),ylim=c(0,28))
  #Add lines describing length of the outbreak
  lines(rep(t0,2),c(0,1e99),col="red",lty=2,lwd=2)
  #lines(rep(t0+l,2),c(0,1e99),col="red",lty=2,lwd=2)
  legend(x="topright",c(expression(Y[t]),expression(hat(lambda)[t]),"exposure point"),col=c("black","blue","red"),lty=c(1,1,2),lwd=c(2,2,1))
}

res <- list("4"=sts.bp.k4)
plotIt1("4")


## ----eval=FALSE---------------------------------------------------------------
## sts2 <- new("sts",epoch=1:length(y.l.ts),observed=matrix(y.l.ts,ncol=1))
## ylim.global <- c(0,30)
## backprojNP(sts2, incu.pmf.vec=incu.pmf, control=modifyList(bp.control,list(k=0,hookFun=plotConvolution)),ylim=c(0,25))
## 


## ----cache=TRUE---------------------------------------------------------------
sts2 <- new("sts",epoch=1:length(y.l.ts),observed=matrix(y.l.ts,ncol=1))
#Backproject using the method by Becker et al. (1991)
res <- list()
for (k in c(0,2,4)) {
  res[[paste(k)]] <- backprojNP(sts2, incu.pmf=incu.pmf, control=modifyList(bp.control,list(k=k)))
}

## -----------------------------------------------------------------------------
#Plot results for example 2
plotIt2 <- function(which) {
  plot(res[[which]],xaxis.labelFormat=NULL,xlab="Time (days)",legend=NULL,dx.upperbound=0,lty=c(1,1,1),lwd=c(1,1,3),main=paste("k=",res[[which]]@control$k,sep=""),ylim=c(0,28))
  #Add lines describing length of the outbreak
  lines(rep(t0,2),c(0,1e99),col="red",lty=2,lwd=2)
  lines(rep(t0+l,2),c(0,1e99),col="red",lty=2,lwd=2)
  legend(x="topright",c(expression(Y[t]),expression(hat(lambda)[t]),"exposure interval"),col=c("black","blue","red"),lty=c(1,1,2),lwd=c(2,2,1))
}


## -----------------------------------------------------------------------------
plotIt2("4")


## ----results='markup'---------------------------------------------------------
#Load data from Zeger, See and Diggle (1989), Stats in Medicine on the
#number of AIDS cases and their delay for north-east MSM
zeger <- read.csv("zeger_etal89-tab2.csv",sep="\t",row.names=1)
colnames(zeger) <- c(0:12)
zeger[rev(seq_len(nrow(zeger))),]


## ----echo=FALSE,tidy=FALSE----------------------------------------------------
#Function to convert reporting triangle matrix into data.frame
matrix2df <- function(zeger) {
  data.frame(n=as.numeric(as.matrix(zeger)),
             t=as.numeric(as.matrix(row(zeger)-1)),
             d=as.numeric(as.matrix(col(zeger)-1)))
}

#Convert to data.frame
zeger.df <- matrix2df(zeger)
#Fit log-linear model.
m <- glm( n ~ as.factor(t) + as.factor(d), data=zeger.df, subset=!is.na(n), family=poisson)

#Prediction m_{t,d} for ALL cells in the contingency table
mu.mle <- predict(m, newdata=zeger.df, type="response")


## ----echo=FALSE,cache=TRUE,tidy=FALSE-----------------------------------------
#Function to compute our target statistic
NtInf <- function(data) {
    as.numeric(with(data, tapply(n, t, sum, na.rm=TRUE)))
}

#Function to generate new data by parametric bootstrap
rntd <- function(data, mle) {
    #Indicator vector of what is observed
    observed <- !is.na(data$n)
    #Extra data copies (one to estimate, one to predict)
    data.estimate <- data.predict <- data

    #Make a new data matrix with observed values replaced
    data.estimate$n[observed] <- rpois(n=nrow(data),lambda=mle)[observed]

    #Fit Poisson GLM to the data to obtain estimates
    m.star <- glm( n ~ as.factor(t) + as.factor(d), data=data.estimate, subset=!is.na(n), family=poisson)

    #Add sampled values where missing
    data.predict$n[!observed] <- rpois(n=nrow(data), predict(m.star,newdata=data,type="response"))[!observed]
    #Done - return new data.frame
    return(data.predict)
}

set.seed(123)
b <- boot::boot(zeger.df, statistic=NtInf,sim="parametric",R=999, ran.gen=rntd,mle=mu.mle)

#Simple percentile intervals
predIntervals <- apply(rbind(b$t0,b$t),2,quantile,prob=c(0.025,0.975))


## -----------------------------------------------------------------------------
#Perform nowcasting only based on the mean
zeger.df2 <- zeger.df
zeger.df2$n[is.na(zeger.df2$n)] <- mu.mle[is.na(zeger.df2$n)]

plot(1:nrow(zeger),b$t0,type="l",ylab="No. cases",ylim=c(0,max(predIntervals)),axes=FALSE,xlab="",lwd=2,col="steelblue")
axis(2)
axis(1,at=1:nrow(zeger),label=rownames(zeger),las=2)
box()
lines(1:nrow(zeger),NtInf(zeger.df2),col="magenta",lwd=2)
matlines(1:nrow(zeger),t(predIntervals),lty=2,col="magenta",lwd=1)
legend(x="topleft",c("observed","delay adjusted"),lty=c(1,1),col=c("steelblue","magenta"),lwd=3)


## ----eval=FALSE---------------------------------------------------------------
## #Just make confidence intervals for mu_{t,d}.
## mu.hat <- predict(m,newdata=zeger.df,type="link",se=TRUE)
## #This ignores that errors are not independent! -> fix this
## waldci.log <- cbind(mu.hat$fit,mu.hat$fit) + t(outer(c(-1,1)*qnorm(1-0.05/2),mu.hat$se.fit))
## waldci <- exp(waldci.log)
## ci.low <- NtInf(within(zeger.df, n <- waldci[,1]))
## ci.up <- NtInf(within(zeger.df, n <- waldci[,2]))
## 
## matlines(1:nrow(zeger),cbind(ci.low,ci.up),col="green")

