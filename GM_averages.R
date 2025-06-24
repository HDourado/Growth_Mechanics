# average plots
setwd(paste(directory,"/Results GM",sep=""))

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
q_inf <- c( mean(aqt[(nt-nt1):nt,r-1]) , mean(aqt[(nt-nt1):nt,r]) )

# estimation of converging q
q_inf <- colSums(aqt[(nt-nt1):nt,-1])/(nt1 + 1)

mu_inf <- mean(amu[(nt-nt1):nt])

# plot limits

xrange <- c(signif(0.9999*q0_OGS[r-1],4),signif(1.0001*q0_OGS[r-1],4))

yrange <- c(signif(0.9999*q0_OGS[r],4),signif(1.0001*q0_OGS[r],4))

# difference from qOGS to q_inf, proportinal to fdifference to initial q0 (percentage)
pdq <- signif( 100*sqrt( sum( (q0_OGS[-1] - q_inf)^2 ) ) /  sqrt( sum( (q0_OGS[-1] - qt[1,-1])^2 ) ) , digits= 2)

# difference in growth rate (percentage)
pdmu <- signif( 100*( mu_inf - muOGS )/muOGS , digits= 2)

png(filename = paste("GM Model ",modelname," average q (",totalT,"h).png",sep=""), units = "in",width=5,height=5, res=600)

par(mfrow=c(1,1))

plot(aqt[(nt-98000):nt,2],aqt[(nt-98000):nt,3], type = "l",xlim=xrange,ylim=yrange,xlab = bquote(q[1]),ylab = bquote(q[2]))
par(new=TRUE)
plot(q0_OGS[2],q0_OGS[3],col=4,pch=19,xlim=xrange,ylim=yrange, axes = FALSE,xlab = bquote(q[1]),ylab = bquote(q[2]))
par(new=TRUE)
plot(q_inf[1],q_inf[2],col=2,pch=19,xlim=xrange,ylim=yrange, axes = FALSE,xlab = bquote(q[1]),ylab = bquote(q[2]))
legend("bottomright",c(paste("d(q*,q_T)/d(q*,q0) =",pdq,"%"),paste("(mu_T - mu*)/mu* =",pdmu,"%" )),col=1:2, cex = 0.8 )

dev.off()

setwd(directory)
