######################################################################
# Helper function to make a plot showing year and week of a surveillance
# time series in interpretable fashion.
######################################################################

plotOne <- function(sts,start,now,doSurv=TRUE,extraYlab="",ylim=c(0,max(observed(sts))),dy.alarm=-0.5) {
  #Remove superfluous margin
  par(mar=c(5,4,0.2,1))
  #Put those value not up to now to zero
  observed(sts)[epoch(sts) > now] <- 0
  
  plot(sts,dx.upperbound=0,legend.opts=NULL,ylab=iconv(paste("No. reported cases",extraYlab,sep=""),"utf-8","latin1"),main="",axes=FALSE,xlab="Year/Reporting Week",xlim=c(0,nrow(sts)),ylim=ylim)

  week <- isoWeekYear(epoch(sts))$ISOWeek
  year <- isoWeekYear(epoch(sts))$ISOYear
  #Where the by 10 divisible weeks
  is.subset1 <- week %% 10 == 0 | week == 1
  #Where the by 5 divisible weeks
  is.subset2 <- week %% 5 == 0 
  #Where the last week of the year
  is.year <- week %% 52 == 0
  #Show axis with this
  axis(1,line=1,cex.axis=0.8,labels=NA,lwd.ticks=0)
  axis(1,(1:length(week))[is.subset1],label=week[is.subset1],line=1,cex.axis=0.8)
  axis(1,(1:length(week))[is.subset2],label=rep("",sum(is.subset2)),line=1,tck=-0.01)
  axis(1,(1:length(week))[is.year]+0.5,label=rep("",sum(is.year)),line=1,tcl=1.3,lwd.ticks=1)
  #Where to plot the years
  where <- tapply(1:length(year),factor(year),mean)
  axis(1,where,label=unique(format(epoch(sts),"%Y")),line=-1.4,lwd.ticks=NA,tick=FALSE)
  axis(2)

  #Show today
  idx <- which(epoch(sts) >= start & epoch(sts) <= now)
  #Upper bound
  points(max(idx),0,pch=20,cex=1.5,col="blue")

  #Do surveillance for this time point
  if (doSurv) {
    s.sts <- farrington(sts, control=list(range=idx,alpha=0.001,b=2,w=3,trend=FALSE,correct=FALSE,limit54=c(0,1)))
    for (i in 1:length(idx)) {
      lines( c(idx[i]-1,idx[i])+0.5, rep(upperbound(s.sts)[i,],2),lwd=3,col="red")
      if (i<length(idx)) {
        lines( rep(idx[i]+0.5,2), upperbound(s.sts)[c(i,i+1),],lwd=1,col="red")
      }
    }
    #Show alarms
    par(xpd=NA)
    points(idx[alarms(s.sts)>0],rep(dy.alarm,sum(alarms(s.sts)>0)),pch=24,cex=1,col="red",lwd=1)
  }
  
  invisible()
}
