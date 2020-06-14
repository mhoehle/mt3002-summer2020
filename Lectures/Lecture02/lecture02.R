## ----echo=FALSE,results='hide',message=FALSE----------------------------------
######################################################################
### Author: Michael Höhle <http://www.math.su.se/~hoehle>
### Course: The Mathematics and Statistics of Infectious Disease Outbreaks
###         (MT3002 - Summer 2020) at the Department of Mathematics,
###         Stockholm University.
###         (https://kurser.math.su.se/course/view.php?id=911)
###         given by Tom Britton and Michael Höhle
###    
### License:
### This work is licensed under a <a rel="license"
### href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
### Attribution-ShareAlike 4.0 International License</a>.
###
### Description:
###  Rnw File for lecture 2
###
### History:
###  -- 10 Jun 2020 file created 
######################################################################

options(width=100)
set.seed(123)
library(tidyverse)
library(RColorBrewer)


## -----------------------------------------------------------------------------
######################################################################
#Laplace transform of stoch. variable X=c=1
######################################################################
Lc <- function(theta) {
 return(exp(-theta*1))
}

#############################################################################
# Find final sizes of Reed Frost chain binomial model by sampling from
# it. 
#
# Parameters:
#  n - Initial number of susceptible
#  m - Initial infected
#  w - contact probability, c.f. Britton2000 p.4
#  samples - number of RF processes to sample.
#
# Note:
# The sampling function is a vector version sampling from all processes
# simultaneously. Ceased processes continue running, and their NA's
# are taken care of.
#############################################################################

fsize.RF <- function(n, m, w, samples) {
  #Initial susceptible
  xj <- matrix(data=n,nrow=samples,ncol=1) 
  #Initial infectives
  yj <- matrix(data=m,nrow=samples,ncol=1) 

  #Loop over all (samples) simulations until they all are ceased.
  while (sum(yj>0) & sum(xj>0)) {
    #Sample from all processes concurrently
    yj <- ifelse(xj > 0, rbinom(samples, xj, 1-(1-w)^yj), 0)
    #Update all xj
    xj <- xj - yj
  }
  #Done 
  return(n-xj) 
}

##################################################################  
# Calculate distribution of final size of standard SIR epidemic.
# C.f. Theorem 2.2.
#
# Parameters:
#   n - population size
#   m - initial number of infectious
#   lambda - transmission parameter
#   phi - Laplace transform of the stoch. var. I
#
# Note:  When n becomes large (n>35-40 in Fig 2.2 example the
#        solve command doesn't work anymore, coz of near singular
#        matrix. More stable implementation is necessary.
#
# Example:
#  Pfinal <- fsize.SIR(35,1,1.5,Lc)
##################################################################  

fsize.SIR <- function(n,m,lambda,phi) {
  #Declare matrices in linear eqn. system Ax=b.
  #Note: size+1 since we also include index 0.
  A <- matrix(data=0,nrow=n+1,ncol=n+1)
  b <- matrix(data=NA,nrow=n+1,ncol=1)

  #Fill them according to Theorem 2.2 in Britton 2000.
  for (l in 0:n) {
    for (k in 0:l) {
    A[l+1,k+1] <- choose(n-k,l-k) / (phi(lambda*(n-l)/n)^(k+m))
    }
    b[l+1] <- choose(n,l)
    #Process indicitor if slow: print(l)
  }
  #Solve the linear equation system. Note: A can be near singular if n
  #is large. Numerical improvement?
  #In R: Use  Pfinal <- forwardsolve(A,b), which gives a little more
  #stability
  Pfinal <- solve(A,b)
  # Done
  return(Pfinal)
}


fsize.hist <- function(n=50,m=1,lambda=1.5,samples=10000,ylim=NULL,plot=TRUE, title=TRUE) {
  #RF by simulation, note q=exp(-lambda/(n)) is the avoidance probability. w=1-q.
  if (n>30) {
    cat("Sample size too large. Using simulation.\n")
    f1 <- fsize.RF(n, m, 1 - exp(-lambda/(n)),samples)
    f2 <- table(factor(f1,levels=0:n))/length(f1)
  } else {
    #Use exact method
    f2 <- fsize.SIR(n,m,lambda,Lc)
  }

  
  #Plot the two for comparison. Note: We plot Z, the final
  #size and not Z+1 as in Fig 2.2
  if (plot) {
      if(title){
    barplot(t(f2),main=paste("w=",sprintf("%.3f",1-exp(-lambda/n)),sep=""),xlab="Final size",names.arg=0:n,ylab="Probability",ylim=ylim)
}
      else{
    barplot(t(f2),main="",xlab="Final size",names.arg=0:n,ylab="Probability",ylim=ylim)
}
      
}
  #Done
  invisible(as.numeric(f2))
}


## -----------------------------------------------------------------------------
#Figure for the presentation
n <- 20
w <- 0.02


## ----fig.width=7,fig.height=3.5-----------------------------------------------
lambdas <- -n*log(1-w)
#lambdas <- seq(0.6,4.32,length=8)
#w <- 1-exp(-lambdas/n) ; lambdas <- -log(1-0.25)*n

maxpmf <- max(sapply(lambdas, function(lambda) {
  max(fsize.hist(n=n,m=1,lambda=lambda,samples=10000,plot=FALSE))
}))

#par(mfrow=c(2,4))
invisible(sapply(lambdas, function(lambda) {
  fsize.hist(n=n,m=1,lambda=lambda,samples=10000,ylim=c(0,maxpmf),plot=TRUE, title=FALSE)
}))



## -----------------------------------------------------------------------------
n <- 100
w <- seq(0.01, 0.045, length=8)
lambdas <- -n*log(1-w)
#lambdas <- seq(1.0,4.0,length=8)
maxpmf <- max(sapply(lambdas, function(lambda) {
  max(fsize.hist(n=n,m=1,lambda=lambda,samples=10000,plot=FALSE))
}))

par(mfrow=c(2,4))
invisible(sapply(lambdas, function(lambda) {
  fsize.hist(n=n,m=1,lambda=lambda,samples=10000,ylim=c(0,maxpmf),plot=TRUE)
}))



## ----echo=TRUE, code = capture.output(dump('fsize.RF', ''))-------------------


## ----echo=TRUE,results='markup',tidy=FALSE------------------------------------
######################################################################
# Likelihood function for the Reed-Frost model
#
# Parameters:
#  w.logit - logit(w) to have unrestricted parameter space
#  x       - vector containing the number of susceptibles at each time
#  y       - vector containing the number of infectious   at each time
#
######################################################################

l <- function(w.logit,x,y) {
  if (length(x) != length(y)) { stop("x and y need to be the same length") }
  K <- length(x)
  w <- plogis(w.logit)
  p <- 1 - (1-w)^y
  return(sum(dbinom( y[-1], size=x[-K], prob=p[-K],log=TRUE)))
}

# Epidemic D in Table 4.1 of Daley and Gani (1999), assuming all susceptibles got infected
y <- c(1, 4, 14, 10, 1, 0)
x <- numeric(length(y))
x[1] <- sum(y[-1])
x[2:length(x)] <- x[1]-cumsum(y[-length(y)])

mle <- optim(par=0,fn=l,method="BFGS",x=x,y=y,control=list(fnscale=-1),hessian=TRUE)
# Maximum likelihood estimator
(w.hat <- plogis(mle$par))

## ----echo=FALSE, results="hide"-----------------------------------------------
# 95% confidence interval
(w.95ci <- plogis( mle$par + c(-1,1)*qnorm(0.975)* sqrt(-1/as.vector(mle$hess))))


## ----fig.align="center"-------------------------------------------------------
x0 <- seq(sum(y[-1]), sum(y[-1])+20, 1)
w.hat <- numeric(length(x0))
log.lik <- numeric(length(x0))
for(i in 1:length(x0)){
    x[1] <- x0[i]
    x[2:length(x)] <- x[1]-cumsum(y[-length(y)])
    
    mle <- optim(par=0,fn=l,method="BFGS",x=x,y=y,control=list(fnscale=-1),hessian=TRUE)
    w.hat[i] <- plogis(mle$par)
    log.lik[i] <- mle$value
}

par(mfrow=c(1,1))
plot(x0, log.lik, type="b", ylab="profile log likelihood", axes=TRUE)


## ----fig.width=7,fig.height=4.5-----------------------------------------------

csfv <- read.table("Data/csfv.txt",col.names=c("t","I"))
pal <- brewer.pal(3,"Set1")
plot(csfv$t,csfv$I,lwd=3,type="l",lty=1,ylab="No. infectious herds",xlab="time (weeks after first infection)",col=pal[1])



## ----echo=FALSE, results="hide"-----------------------------------------------
# Load pkg.
suppressPackageStartupMessages(library(deSolve))


## ----SIRDIFF------------------------------------------------------------------
##############################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list containing dS(t)/dt and dI(t)/dt
##############################################################################

sir <- function(t,y, parms) {
  beta <- parms[1]
  gamma <- parms[2]
  S <- y[1]
  I <- y[2]
  return(list(c(S=-beta*S*I,I=beta*S*I-gamma*I)))
}

## -----------------------------------------------------------------------------
gamma <- 0.3
## beta.grid <- c(0.001,0.00005, 0.00003)
beta.grid <- c(1.5e-04, 4.5e-05, 3.0e-05)
N <- 20000
times <- seq(0,50,length=1000)

#R0 for the 3 examples
(R0 <- beta.grid*N/gamma)


## ----fig.width=7,fig.height=4.5-----------------------------------------------
I <- sapply(beta.grid, function(beta) {
  lsoda(y=c(N-1,1), times=times, func=sir,parms=c(beta,gamma))[,3]
})
pal <- brewer.pal(length(beta.grid),"Set1")
matplot(times,I,type="l",lwd=3,lty=1,col=pal,xlab="time",ylab="I(t)", ylim=c(0, N))
#legend(x="topright",paste("beta=",sprintf("%.5f",beta.grid),sep=""),lty=1,lwd=3,col=pal)
#Expressions in the legend
beta.grid.str <- format(beta.grid, scientific=TRUE)
leg <- c(as.expression(substitute(beta==val,list(val=beta.grid.str[1]))),as.expression(substitute(beta==val,list(val=beta.grid.str[2]))),as.expression(substitute(beta==val,list(val=beta.grid.str[3]))))
legend(x="topright",legend=leg,lty=1,lwd=3,col=pal)


## ----echo=TRUE----------------------------------------------------------------
##############################################################################
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list containing dS(t)/dt and dI(t)/dt
##############################################################################

sir <- function(t,y, parms) {
  beta <- parms[1]
  gamma <- parms[2]
  S <- y[1]
  I <- y[2]
  return(list(c(S=-beta*S*I,I=beta*S*I-gamma*I)))
}


## ----echo=TRUE, results="show"------------------------------------------------
sim <- lsoda(y=c(N-1,1), times=times, func=sir,parms=c(beta.grid[1],gamma))
head(sim, n=3)


## ----echo=TRUE----------------------------------------------------------------
# Step width of the Euler method
h <- 0.1
y <- matrix(NA, nrow=ceiling(20/h), ncol=3, dimnames=list(NULL, c("t","S","I")))
# Initial value
y[1,] <- c(0,N-1,1)
# Loop
for (i in 2:nrow(y)){
  y[i,] <- c(y[i-1,"t"]+h, y[i-1,c("S","I")] +
         h * sir(y[i-1,"t"], y[i-1,c("S","I")], parms=c(beta.grid[1],gamma))[[1]])
}


## -----------------------------------------------------------------------------
# Show Euler solve
plot(y[,"t"], y[,"I"], type="l", xlab="time", ylab="I(t)")
# Add lsoda (which uses a more advanced method)
lines(times, lsoda(y=c(N-1,1), times=times, func=sir,parms=c(beta.grid[1],gamma))[,3], col="#377EB8")
legend(x="topright", c("Euler-method (h=0.1)", "lsoda"), col=c(1,"#377EB8"), lty=1)
title(substitute(beta == a, list(a=beta.grid[1])))


## ----fig.width=7,fig.height=4.5, fig.align="center"---------------------------
S <- sapply(beta.grid, function(beta) {
  lsoda(y=c(N-1,1), times=times, func=sir,parms=c(beta,gamma))[,2]
})
pal <- brewer.pal(length(beta.grid),"Set1")
matplot(times,S,type="l",lwd=3,lty=1,col=pal,xlab="time",ylab="S(t)")
legend(x="topright",legend=leg,lty=1,lwd=3,col=pal)


## ----echo=FALSE---------------------------------------------------------------

N <- 21500
#sumy <- sum(csfv$I)
sumy <- 429
f <- sumy/N
R0 <- -log(1-f)/f



## ----echo=TRUE----------------------------------------------------------------
######################################################################
#Least-squares fit
######################################################################

ll.gauss <- function(theta, take.sqrt=FALSE) {
  #Solve ODE using the parameter vector theta
  res <- lsoda(y=c(N-1,1), times=csfv$t, func=sir, parms=exp(theta))
  #Squared difference?
  if (take.sqrt==FALSE) {
    return(sum(dnorm(csfv$I,mean=res[,3],sd=1,log=TRUE)))
  } else { 
    return(sum(dnorm(sqrt(csfv$I),mean=sqrt(abs(res[,3])),sd=1,log=TRUE)))
  }
}


## ----echo=TRUE----------------------------------------------------------------
#Determine MLE
N <- 21500
mle <- optim(log(c(0.00002,3)), fn=ll.gauss,control=list(fnscale=-1))

#Show estimates and resulting R0 estimate
beta.hat <- exp(mle$par)[1]
gamma.hat <- exp(mle$par)[2]
R0.hat <- beta.hat*N/gamma.hat


## ----echo=TRUE, results="show"------------------------------------------------
mu <- lsoda(y=c(N-1,1), times=csfv$t, func=sir,parms=exp(mle$par))
head(mu, n=3)


## ----echo=FALSE---------------------------------------------------------------
mle2 <- optim(log(c(0.00002,3)), fn=ll.gauss, take.sqrt=TRUE, control=list(fnscale=-1))
beta.hat2 <- exp(mle2$par)[1]
gamma.hat2 <- exp(mle2$par)[2]
R0.hat2 <- beta.hat2*N/gamma.hat2
mu2 <- lsoda(y=c(N,1), times=csfv$t, func=sir,parms=exp(mle2$par))


## ----fig.width=7,fig.height=4.0, fig.align="center"---------------------------
pal <- brewer.pal(4,"Set1")
matplot(mu[,1],cbind(csfv$I,mu[,3],mu2[,3]),type="l",lwd=3,lty=1,ylab="No. infectious herds",xlab="time (weeks after first infection)",col=pal)
legend(x="topright",c("CSFV outbreak","LS fit", "LS-sqrt fit"), lty=1,col=pal,lwd=3)


## ----echo=FALSE, warning=FALSE------------------------------------------------
######################################################################
# Poisson likelihood
######################################################################

ll.pois <- function(theta) {
  #Solve ODE using the parameter vector theta
  res <- lsoda(y=c(N-1,1), times=csfv$t, func=sir, parms=exp(theta))
  #Poisson likelihood
  return(sum(dpois(round(csfv$I),lambda=res[,3],log=TRUE)))
}

#Determine MLE
mle <- optim(mle$par, fn=ll.pois,control=list(fnscale=-1))

#Show results
(beta.hat <- exp(mle$par)[1])
(gamma.hat <- exp(mle$par)[2])
(R0.hat <- beta.hat*N/gamma.hat)
muP <- lsoda(y=c(N-1,1), times=csfv$t, func=sir,parms=exp(mle$par))

## ----fig.width=7,fig.height=4.0-----------------------------------------------
matplot(mu[,1],cbind(csfv$I,mu[,3],mu2[,3], muP[,3]),type="l",lwd=3,lty=1,ylab="No. infectious herds",xlab="time (weeks after first infection)",col=pal)
legend(x="topright",c("CSFV outbreak","LS fit", "LS-sqrt fit", "Poisson fit"),lty=1,col=pal,lwd=3)


## ----echo=FALSE,eval=FALSE----------------------------------------------------
## #Estimate R0 from final size data
## f <- 429/N
## (R0.hat.finalsize <- -log(1-f)/f)


## ----cache=TRUE---------------------------------------------------------------
######################################################################
# Simulate simple SIR model using the Gillespie-Algorithm
#
# Params:
#  T - max time
#  lambda - transmission rate
#  gamma - recovery rate (mu(X) = d*X) 
#  n - initial number of susceptibles. Always initially one infective
#  m - initial number of infectives
######################################################################

rSIR <- function(T, beta, gamma, n, m) {
  #Initialize (x= number of susceptibles)
  t <- 0
  x <- n
  y <- m

  #Possible events
  eventLevels <- c("S->I","I->R")
  #Initialize result
  events <- data.frame(t=t,x=x,y=y,event=NA)
  #Loop until we are past time T
  while (t < T & (y>0)) {
      #Draw the waiting type for each possible event
      wait <-  rexp(2,c("S->I"=beta*x*y,"I->R"=gamma*y))
      #Which event occurs first
      i <- which.min(wait)
      #Advance Time
      t <- t+wait[i]
      #Update population according to the eventy type
      if (eventLevels[i]=="S->I") { x <- x-1 ; y <- y+1}
      if (eventLevels[i]=="I->R") { y <- y-1 }
      #Store result
      events <- rbind(events,c(t,x,y,i))
  }
  #Recode event type and return
  events$event <- factor(eventLevels[events$event], levels=eventLevels)
  return(events)
}

##Select parameters as in the example
beta <- 0.01
gamma <- 0.5
S0 <- 100
(R0 <- beta*S0/gamma)
T <- 500
nSim <- 10
#Do a few simulations
set.seed(124)
trajs <- lapply( 1:nSim, function (x) {
    rSIR(T, beta=beta, gamma=gamma, n=S0, m=1)
})


## -----------------------------------------------------------------------------
for (i in 1:length(trajs)) {
    plotf <- if (i==1) plot else lines
    plotf(trajs[[i]]$t,trajs[[i]]$x,type="s",ylim=c(0,S0),xlab="Time",ylab="Susceptibles")
    #lines(c(trajs[[i]]$t[nrow(trajs[[i]])],1e99),rep(trajs[[1]]$x[nrow(trajs[[i]])],2))
    points(trajs[[i]]$t[nrow(trajs[[i]])],trajs[[i]]$x[nrow(trajs[[i]])],cex=1,pch=20)
}


## ----echo=TRUE, code = capture.output(dump('rSIR', ''))-----------------------


## ----COMPARE------------------------------------------------------------------
compare <- function(S0) {
  # Simulate 250
  trajs <- lapply(1:250, function(x) rSIR(T, beta=beta, gamma=gamma, n=S0, m=1))
  max_T <- max(unlist(lapply(trajs, function(traj) max(traj$t))))
  max_Y <- max(unlist(lapply(trajs, function(traj) max(traj$y))))
  for (i in 1:length(trajs)) {
    plotf <- if (i==1) plot else lines
    plotf(trajs[[i]]$t,trajs[[i]]$y,type="s",ylim=c(0,max_Y),xlab="Time",ylab="Susceptibles",xlim=c(0,max_T),col=rgb(0,0,0,0.2))
  }
  # Solve ODE
  sim <- lsoda(y=c(S0,1), times=times, func=sir,parms=c(beta,gamma))
  lines(sim[,"time"], sim[,3], col="#377EB8", lwd=2)
  invisible()
}
compare(S0=S0)

## ----eval=FALSE---------------------------------------------------------------
## # Compare for different S0
## par(mfcol=c(2,1))
## compare(S0=100)
## title("S(0) = 100, m=1")
## compare(S0=500)
## title("S(0) = 500, m=1")

