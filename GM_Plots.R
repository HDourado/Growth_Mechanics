# plot results #######################################################################################################

setwd(paste(directory,"/Results GM",sep=""))

#pdf(paste("GM Model ",modelname," (T=",totalT,").pdf",sep=""),title="Optimization results", width=10,height=5)

png(paste("GM Model ",modelname," (T=",totalT,").png",sep=""), width=10,height=5,res=600, units = "in")
par(mfrow=c(1,1))

maxmu <- max(mu_opt)*1.1

# graphic
coloraxis <- "darkgrey"
colorlab <- "black"

# Growth rate vs. log10(first external concentration) ################################################################ 

mu_opt   <- opt_state[,"mu"]
mainx <- opt_state[,4]

plot(time,mainx,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 1,
     ylim=c(0,1.1*max(mainx)),xlim=c(0,totalT), ylab=paste("Concentration",reactant[1]," (g/L)"), 
     xlab=bquote("Time" ~ (h)), yaxs="i", xaxs="i")

# Growth rate vs. time ###############################################################################################

mu0 <- signif(mu_opt[1],digits = 4)

muT <- signif(mu_opt[length(mu_opt)],digits=4)

plot(time,mu_opt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 1,
     ylim=ylimits(mu_opt),xlim=c(0,totalT), ylab=bquote("Specific growth rate " ~ mu ~~ (h^-1)),
     xlab=bquote("Time" ~ (h)), yaxs="i", xaxs="i",type = "l",main=paste("mu(0)=",
     mu0,", mu(T)=",muT, ", mean mu = ",signif(mean(mu_opt),digits = 7), ", muOGS =",muOGS))
#axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
if (const_env==1) abline(a=muOGS,b=0,col=2,lty=2)
#abline(a=mean(mu_opt),b=0,col=1,lty=2)

# Dynamic trajectory fitness vs. Balanced growth fitness #######################

fitness_ratio <- 0
for (t in 1:nt) fitness_ratio[t] <- mean(mu_opt[1:t])/muOGS

plot(time,fitness_ratio,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= 2,
     ylab=bquote("Fitness(OGT)/Fitness(OGS)" ), ylim=c(0.95*min(fitness_ratio),1.05*max(fitness_ratio)),
     xlab=bquote("Time" ~ (h)), yaxs="i", xaxs="i",type = "l")
#axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, at=seq(minx1,maxx1,1), mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
abline(a=1,b=0,col=1,lty=2)

# long term fitness in detail #######################################################################################
nperiods <- floor(totalT/mu_T)

amu_period <- amu[c(1:nperiods)*mu_T/delta_t]

plot(1:nperiods, amu_period, xlab="Number of periods",ylab="Average fitness Bar(mu)") 

# proportions
plot(1:nperiods, amu_period/muOGS, xlab="Number of periods",ylab="Bar(mu)/mu*",
     ylim=c(mu_E0/muOGS,max(amu_period/muOGS)), main="mu* (blue), average mu (black), mu_E (red)" ) 
abline(a=mu_E0/muOGS,b=0,col=2,lty=2)
abline(a=1,b=0,col=4,lty=2)
abline(a=last(amu_period)/muOGS,b=0,col=1,lty=2)

# Internal reactant concentrations vs.mu #############################################################################

colori <- rainbow(p)

if (p > 8) colori <- c(1:8, rainbow(p-8))

maxci <- 1.2*max(cit)

matplot(time,cit,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colori,
        ylim=c(0,maxci),xlim=c(0,totalT),xlab=bquote("Time" ~ (h)), 
        ylab=bquote("Internal concentrations " ~ ~ c^i ~ " (g/L)"), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",i_reactant,col=colori,pch=16,cex=0.7, horiz = T)

# Metabolite concentrations vs. mu ###################################################################################

maxcm <- 1.2*(max(cmt))

matplot(time,cmt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab,xlim=c(0,totalT), frame.plot=F, 
        col= colori,ylim=c(0,maxcm),  xlab=bquote("Time" ~ (h)),
        ylab=bquote("Metabolite concentrations " ~ c^m ~ " (g/L)"), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",i_reactant[-p],col=colori,pch=16,cex=0.7, horiz = T)

# Now the averages #############################################################

# average cm up to each point in time
acmt <- cmt
for (t in 2:nt) acmt[t,] <- colSums(cmt[1:t,]/t)

c_OGS <- rho*b(q0_OGS)

matplot(time,cmt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab,xlim=c(0,totalT), frame.plot=F,
        col= colori,ylim=c(0,maxcm),  xlab=bquote("Time" ~ (h)),
        ylab=bquote("Metabolite concentrations " ~ c^m ~ " (g/L)"), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",i_reactant[-p],col=colori,pch=16,cex=0.7, horiz = T)
for(i in 1:(p-1)) abline(a = c_OGS[i],b=0,col=colori[i],lty=2)

for (i in 1:(p-1)) {

  plot(1:nperiods,acmt[1:nperiods ,i],ylab=i_reactant[i] )
  abline(a = c_OGS[i],b=0,col=colori[i],lty=2)
  legend("left",paste("average c",i))

}

# Protein concentrations #############################################################################################

colorj <- 1:r

if (r > 8) colorj <- c(1:8, rainbow(r-8))

pt <- opt_state[,paste("p",reaction)]

maxp <- ceiling(max(pt)*1.05)

minp <- floor(min(pt)*0.95)

p_OGS <- prot(q0_OGS,0*q0,0)

matplot(time,pt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        col= colorj,ylim=c(minp,maxp), xlab=bquote("Time" ~ (h)),
        ylab=bquote("Protein concentrations " ~  p ~ " (g/L)"), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# Now the averages #############################################################

# average p up to each point in time
apt <- pt
for (t in 2:nt) apt[t,] <- colSums(pt[1:t,]/t)

matplot(time,apt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        col= colorj,ylim=c(minp,maxp), xlab=bquote("Time" ~ (h)),
        ylab=bquote("Protein concentrations " ~  p ~ " (g/L)"), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)
for(j in 1:r) abline(a = p_OGS[j],b=0,col=colorj[j],lty=2)

for (j in 1:r) {
  
  plot(1:nperiods,apt[1:nperiods ,j], main=paste("average p", j) )
  abline(a = p_OGS[j],b=0,col=colorj[j],lty=2)

}

# chis #########################################################################

maxchit <- max(chit)*1.1

matplot(time,chit,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        col= colorj,ylim=c(0,maxchit), xlab=bquote("Time" ~ (h)), 
        ylab=bquote("Ribosome allocation " ~  chi ), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# Protein fractions ##################################################################################################

maxp <- ceiling(max(pt)*1.1)

matplot(time,phit,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,
        ylim=c(0,1.1*max(phit)),xlab=bquote("Time" ~ (h)), ylab=bquote("Proteome fractions " ~  phi ), 
        yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# Fluxes vs. mu ######################################################################################################

maxv <- ceiling(max(vt)*1.1)

matplot(time,vt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,ylim=c(0,maxv), 
        xlab=bquote("Time" ~ (h)), ylab=bquote("Fluxes " ~  v ~ " (g/L/h)"), yaxs="i", xaxs="i",type = "l")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topleft",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# q vs. mu ##############################################################################################

qt <-  opt_state[,paste("q",reaction)]

matplot(time,qt,pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,
        ylim=c(0,1.1*max(qt)),xlim=c(0,totalT),xlab=bquote("Time" ~ (h)), ylab=bquote("Generalized coordinate" ~  q), 
        yaxs="i", xaxs="i",type = "l")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# dq vs. mu ##############################################################################################

dqt <-  opt_state[,paste("dq",reaction)]

matplot(time,dqt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,
        xlab=bquote("Time" ~ (h)), ylab=bquote("Generalized velocity" ~  dq), 
        yaxs="i", xaxs="i",type = "l")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)


# Turnover times ####################################################################################################

tau_opt <- opt_state[,paste("tau",reaction)]

maxtau <- 1.2*(max(tau_opt))

matplot(time,tau_opt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,ylim=c(0,maxtau), 
        xlab=bquote("Time" ~ (h)), ylab=bquote("Turnover times " ~  tau ~ (h)), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topright",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# Turnover frequencies ##############################################################################################

k_opt <- 1/tau_opt

maxf <- 1.2*(max(k_opt))

matplot(time,k_opt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj,ylim=c(0,maxf), 
        xlab=bquote("Time" ~ (h)), ylab=bquote("Apparent turnover numbers " ~  k[app] ~ (h^-1)), yaxs="i", xaxs="i")
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
#axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
legend("topleft",reaction,col=colorj,pch=16,cex=0.7, horiz = T)

# Singular Growth Modes ############################################################################################

kappat <- t(t(V)%*%t(qt))

sgms <- paste("x",c(1:r))

matplot(time,kappat,xlim=c(0,totalT),ylim=c(min(kappat) - 0.1,1.1*max(kappat)),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, col= colorj, 
        xlab=bquote("Time" ~ (h)), ylab=bquote("x"), yaxs="i", xaxs="i",main="Singular Growth Modes")
legend("topright",sgms,col=colorj,pch=16,cex=0.7)


# Growth laws ##################################################################

lowmu  <-  max(mu_opt)/4 # min(mu_opt) 

highmu <- max(mu_opt) #  4*max(mu_opt)/5

mediummu <- lowmu < mu_opt & mu_opt < highmu

phi0 <- 0
slope <- 0
for (j in 1:r) {
  
  maxphi <- 1.1*(max(phit[,j]))
  
  philm <- lm(phit[mediummu,j] ~ mu_opt[mediummu] )
  
  phi0[j] <- signif(philm$coefficients[1],3)
  
  slope[j] <- signif(philm$coefficients[2],3)
  
  # plot #######################################################################
  
  maxmu <- max(mu_opt)*1.1
  
  maxphi <- 1.1*max(phit[,j])
  
  # main=paste(reaction[j],": offset =",phi0[j],", slope = ",slope[j])
  matplot(mu_opt,phit[,j],pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
          main=paste(reaction[j], "1/slope:",signif(1/slope[j],digits=3),"kcatf:",
                     signif(kcatf[j],digits = 3)),
          col= colorj[j],xlim=c(0,maxmu),ylim=c(0,maxphi), xaxt = "n", yaxt = "n", 
          xlab=bquote("Growth rate"~ mu ~ (h^-1)), ylab=bquote("Proteome fraction " ~  phi ), 
          yaxs="i", xaxs="i")
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
  axis(1, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=0.02,las=1)
  axis(2, mgp=c(3, 0.4, 0), col=coloraxis, col.ticks=coloraxis, col.axis=coloraxis, tck=-0.02,las=1)
  #legend("top",reaction,col=colorj,pch=16,cex=0.7)
  if (!is.na(slope[j])) abline(philm)
  #abline(v=lowmu,lty=2)
  #text(0.5*lowmu,0.95*maxphi,labels=bquote(0.2*mu[max]))
  #abline(v=highmu,lty=2)
  #text(1.13*highmu,0.95*maxphi,labels=bquote(0.8*mu[max]))
  #abline(a=ephi0[j],b=0,lty=2,col=3)
  #text(1.13*highmu,1.1*ephi0[j],labels=ephi0[j])
  
}

# Phase space #######################################################################################################

qt <- opt_state[,paste("q",reaction)]

dqt <- opt_state[,paste("dq",reaction)]

matplot(qt,dqt,type = "l", xlab=bquote(q(t)),ylab=TeX(r"($\dot{q}(t)$)"))
legend("topleft",reaction,col=colorj,pch=16,cex=0.7)

#3D plot #######################################################################

#if (r == 4) scatter3D(qt[,3], qt[,2], qt[,4], colvar=time)

if (r == 4) {
  
  lines3D(qt[,3], qt[,2], qt[,4], colvar=time, xlab="q3",ylab="q2",zlab="q4")


  lines3D(qt[,1], qt[,2], qt[,3], colvar=time, xlab="q1",ylab="q2",zlab="q3")
  
  lines3D(qt[,2], qt[,1], qt[,4], colvar=time, xlab="q2",ylab="q1",zlab="q4")
  
  lines3D(qt[,4], qt[,2], qt[,3], colvar=time, xlab="q4",ylab="q2",zlab="q3")
  
  
}

if (r > 2) lines3D(qt[,r-2], qt[,r-1], qt[,r], colvar=time, xlab=paste("q_",r-2,sep=""),
                   ylab=paste("q_",r-1,sep=""),zlab=paste("q_",r,sep=""))


 lines3D(time, qt[,r], mu_opt, colvar=time, 
         xlab=c("time (h)"), ylab="q_r",zlab="mu (1/h)")
 
 lines3D(time, qt[,r-1], qt[,r], colvar=time, 
         xlab=c("time (h)"), ylab=paste("q_",r-1,sep=""),zlab=paste("q_",r,sep=""))
 
 lines3D(qt[,r], qt[,r-1],mu_opt, colvar=time, xlab=paste("q_",r-1,sep=""),ylab=paste("q_",r,sep=""),
         zlab="mu (1/h)")
 
 lines3D(qt[,r], qt[,r-1],amu, colvar=time, xlab=paste("q_",r-1,sep=""),ylab=paste("q_",r,sep=""),
         zlab="amu (1/h)")
 
 lines3D(A_it[,r-2], A_it[,r-1],mu_opt, colvar=time, xlab=paste("p_",r-2,sep=""),ylab=paste("p_",r-1,sep=""),
         zlab="mu (1/h)")
 
 lines3D(dqt[,r], dqt[,r-1],mu_opt, colvar=time, xlab=paste("dq_",r-1,sep=""),ylab=paste("dq_",r,sep=""),
         zlab="mu (1/h)")
 
 lines3D(phit[,r-1], phit[,r],mu_opt, colvar=time, xlab=paste("phi_",r-1,sep=""),ylab=paste("phi_",r,sep=""),
         zlab="mu (1/h)")
 
 lines3D(vt[,r-1], vt[,r],mu_opt, colvar=time, xlab=paste("v_",r-1,sep=""),ylab=paste("v_",r,sep=""),
         zlab="mu (1/h)")
 
 lines3D(tau_opt[,r-1], tau_opt[,r],mu_opt, colvar=time, xlab=paste("tau_",r-1,sep=""),ylab=paste("tau_",r,sep=""),
         zlab="mu (1/h)")
 
 lines3D(kappat[,r-1], kappat[,r],mu_opt, colvar=time, xlab=paste("z_",r-1,sep=""),ylab=paste("z_",r,sep=""),
         zlab="mu (1/h)")
 
 ##### 
 
 if (r == 3) {
 
 # q2,q3 range
 
 q1 <- seq(min(qt[,2]*0.99),max(qt[,2]*1.01),by=q0_OGS[r-1]/1000) 
 
 q2 <- seq(min(qt[,3]*0.99),max(qt[,3]*1.01),by=q0_OGS[r]/1000) 
 
 Eq <- q1%*%t(q2)*0
 
 for (i in 1:dim(Eq)[1]) for (j in 1:dim(Eq)[2]) Eq[i,j] <- mu_E(c(1, q1[i], q2[j] ),nt)
 
 contourplot(Eq,cuts=50, xlab="q_1",ylab="q_2",row.values =q1, column.values =q2, xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)), labels = F  )
 
 plot(qt[,2],qt[,3], xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2)) ,asp=1, axes = FALSE, col= 2)
 
 levelplot(Eq,cuts=10, xlab="q_1",ylab="q_2",row.values =q1, column.values =q2, xlim=c(min(q1),max(q1)), ylim=c(min(q2),max(q2))   )
 
}
 
# Plots Model ##################################################################

par(mfrow=c(1,1))

size <- 5/r

textplot(M,cex=size)
title("M")

textplot(K,cex=size)
title("K")

textplot(KA,cex=size)
title("KA")

textplot(rbind(kcatf,kcatb),cex=size)
title("kcat")

# Calculates the corresponding equilibrium constants Keq
Keq <- 0
for (j in 1:r) Keq[j] <- kcatf[j]*prod(KP[KP[,j] < Inf,j])/(kcatb[j]*prod(KS[KS[,j] < Inf,j]))

textplot(Keq,cex=size)
title("Keq")

# Plots Parameters #############################################################

textplot(q0,cex=size)
title("q0")

textplot(dq0,cex=size)
title("dq0")

textplot(delta_t,cex=1)
title("delta_t")

textplot(rtol,cex=1)
title("rtol")

dev.off()

setwd(directory)

