# average plots
setwd(paste(directory,"/Results GM",sep=""))

# long term fitness estimation #################################################

nperiods <- floor(totalT/mu_T) # number of cycles/periods

amu_period <- amu[c(1:nperiods)*mu_T/delta_t] # average mu in each period

Lambda <- amu_period[nperiods]

nu <- mu_f

# long term fitness in detail ##################################################

png(filename = paste("GM Model ",modelname," average mu periods (",totalT,"h).png",sep=""), units = "in",width=5,height=5, res=600)

# mean amplitude of mu oscillation
amplitude <- mean(abs(mu_opt - muOGS))

plot(c(1:nperiods)*mu_T, amu_period, xlab="Number of periods",ylab="Average fitness Bar(mu)",
     pch=19,xlim=c(0,totalT),ylim=c(varphi0, muOGS + 0.1*(muOGS - varphi0) ) ) #,ylim=c(varphi0, 1.001*max(amu_period))
abline(a=muOGS,b=0, col=2)
abline(a=varphi0,b=0, col=4)
legend("bottomright",c(paste("mu*       =",muOGS),paste("Lambda =",Lambda),paste("varphi0  =",varphi0), paste("amplitude =",amplitude)  ),col=c(2,1,4,1), cex = 0.8 , pch=16)

dev.off()

if (totalT >= 100) {
  
  #pdf(paste("GM Model ",modelname," averages (T=",totalT,").pdf",sep=""),title="Optimization results",
  #    width=5,height=8)
  
  png(paste("GM Model ",modelname," averages (T=",totalT,").png",sep=""),res=600, units = "in", width=5,height=8)
  
  par(mfrow=c(1,1))
  
  # calculate average state (t)
  
  aqt <- qt
  for (t in 2:nt) aqt[t,] <- colSums(qt[1:t,])/t
  
  # last 1% of nt
  nt1 <- round(nt/100)
  
  # estimation of converging q
  q_inf <- colSums(aqt[(nt-nt1):nt,])/(nt1 + 1)
  
  # mu of q_inf
  mu_inf <- mean(amu[(nt-nt1):nt])
  
  # plot limits
  
  xrange <- c(signif(0.9999*q0_OGS[r-1],4),signif(1.0001*q0_OGS[r-1],4))
  
  yrange <- c(signif(0.9999*q0_OGS[r],4),signif(1.0001*q0_OGS[r],4))
  
  # difference from qOGS to q_inf, proportinal to fdifference to initial q0 (percentage)
  pdq <- signif( 100*sqrt( sum( (q0_OGS[-1] - q_inf[-1])^2 ) ) /  sqrt( sum( (q0_OGS[-1] - qt[1,-1])^2 ) ) , digits= 2)
  
  # estimating Lambda by restrictinf time T to a multiple of the period
  number_timepoints <- floor(totalT/mu_f)
  
  #Lambda <<- mean(mu_opt)
  
  # difference in growth rate (percentage)
  pdmu <- signif( 100*( Lambda - muOGS )/muOGS , digits= 2)
  
  ################################################################################
  
  png(filename = paste("GM Model ",modelname," average q (",totalT,"h).png",sep=""), units = "in",width=5,height=5, res=600)
  
  par(mfrow=c(1,1))
  
  nt_final <- 0.98*nt
  
  plot(aqt[(nt-nt_final):nt,r-1],aqt[(nt-nt_final):nt,r], type = "l",xlim=xrange,ylim=yrange,xlab = bquote(q[r-1]),ylab = bquote(q[r]))
  par(new=TRUE)
  plot(q0_OGS[r-1],q0_OGS[r],col=4,pch=19,xlim=xrange,ylim=yrange, axes = FALSE,xlab = bquote(q[r-1]),ylab = bquote(q[r]))
  par(new=TRUE)
  plot(q_inf[r-1],q_inf[r],col=2,pch=19,xlim=xrange,ylim=yrange, axes = FALSE,xlab = bquote(q[r-1]),ylab = bquote(q[r]))
  legend("bottomright",c(paste("d(q*,q_T)/d(q*,q0) =",pdq,"%"),paste("(Lambda - mu*)/mu* =",pdmu,"%" )),col=1:2, cex = 0.8 )
  
  dev.off()
  
  ################################################################################
  
  png(filename = paste("GM Model ",modelname," average mu (",totalT,"h).png",sep=""), units = "in",width=10,height=5, res=600)
  
  plot(time,amu,xlim=c(0,totalT),ylim = c(0.99999*muOGS,1.00001*muOGS),xlab=bquote("Time" ~ (h)), ylab=bquote("Average" ~ mu ~ (h^-1)))
  par(new=TRUE)
  abline(a=muOGS,b=0, col=2)
  
  dev.off()
  
}

setwd(directory)
