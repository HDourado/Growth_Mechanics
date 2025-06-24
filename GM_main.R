# plot results #######################################################################################################

setwd(paste(directory,"/Results GM",sep=""))

#pdf(paste("GM Model ",modelname," main (",totalT,"h).pdf",sep=""),title="Optimization results",
#    width=5,height=8)

png(filename = paste("GM Model ",modelname," main (",totalT,"h).png",sep=""), units = "in",width=5,height=8, res=600)

par(mfrow=c(1,1))

# plots ########################################################################
par(mfcol = c(6, 1), mar = numeric(4), oma = c(4, 4, .5, .5), 
    mgp = c(2, .6, 0))

colorj <- 1:r
if (r > 8) colorj <- c(1:8, rainbow(r-8))

if (r == 3) colorj <- c(3,4,2)

colorm <- (r+1):(r+p)

coloraxis <- "darkgrey"
colorlab <- "black"

# function to calculate y-axis limits in the plot automatically
ylimits <- function(something) {
  
  yrange <- max(something) - min(something)
  
  if (yrange == 0) yrange <- 1
  
  ymax <- max(something) + 0.1*yrange
  
  ymin <- max( min(something) - 0.1*yrange, 0)
  
  return(c(ymin,ymax))
}

# function to calculate y-axis limits in the plot automaticaly
ylimitslab <- function(something) {
  
  yrange <- max(something) - min(something)
  
  if (yrange == 0) yrange <- 1
  
  ymax <- max(something) + 0.5*yrange
  
  ymin <- max( min(something) - 0.1*yrange, 0)
  
  return(c(ymin,ymax))
}

# average mu
amu <- mu_opt
for (t in 1:nt) amu[t] <- sum(mu_opt[1:t])/t

cols <- colorRampPalette(c( "green","darkgreen"))

cole <- colorRampPalette(c( "deepskyblue","navyblue"))

colorj <- c(cols(n_tr),cole(ne),2)

colm <- colorRampPalette(c( "lightgrey","grey50"))

colorm <- colm(p-1)

if (p > 4) colorm <- rainbow(p-1)

# moving average function
moving_av <- function(data_in,window_hour) {
  
  extratime <- window_hour/delta_t
  
  data_new <- c(rep(data_in[1],extratime),data_in,rep(last(data_in),extratime))
  
  moving <- 0
  for (t in 1:nt) moving[t] <- mean( data_new[ t:(t + extratime) ] )
  
  return(moving)
}

if (sum(testx == 0) == nt) const_env <<- 1

# plots
plot(time,a1,xlim=c(0,totalT), ylim=ylimits(a1),
     ylab=bquote("External concentration " ~ x^1 ~~ (a.u.)), yaxs="i", xaxs="i",
     type = "l", axes = FALSE)
axis(2L)
box()
# growth rate
plot(time,mu_opt, xlim=c(0,totalT), ylim=ylimits(mu_opt),
     yaxs="i", xaxs="i", ylab=bquote("Specific growth rate " ~ mu ~~ (h^-1)),type = "l", axes = FALSE)
lines(time,amu,lty=2)
if (const_env==1) abline(a=muOGS,b=0,col="grey")
axis(2L)
box()
matplot(time,phit[,r],xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorj[r],ylim=ylimits(phit[,r]), 
        ylab=bquote("Proteome fractions " ~  phi ), yaxs="i", xaxs="i", axes = FALSE)
axis(2L)
box()
matplot(time,phit[,-r],xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorj[-r],ylim=ylimitslab(phit[,-r]), 
        ylab=bquote("Proteome fractions " ~  phi ), yaxs="i", xaxs="i", axes = FALSE)
legend("topright",reaction[-r],col=colorj[-r],pch=16,cex=1,horiz=T)
axis(2L)
box()
matplot(time,bmt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorm,ylim=ylimitslab(bmt), 
        xlab=bquote("Time" ~ (h)), yaxs="i", xaxs="i", axes = FALSE)
axis(2L)
box()
legend("topright",metabolite,col=colorm,pch=16,cex=1,horiz=T)
matplot(time,chit,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        col= colorj,ylim= ylimitslab(chit), xlab=bquote("Time" ~ (h)), 
        ylab=bquote("Ribosome allocation " ~  chi ), yaxs="i", xaxs="i")
box()
legend("topright",reaction,col=colorj,pch=16, horiz = T)
mtext(bquote("Time" ~ (h)), side = 1, outer = TRUE, line = 2.2)
mtext(bquote(a[C] ~~ (g/L)), side = 2, outer = TRUE, line = 2.2, adj = 9.5/10)
mtext(bquote(mu ~~ (h^-1)), side = 2, outer = TRUE, line = 2.2, adj = 7.8/10)
mtext(bquote(phi[r])       , side = 2, outer = TRUE, line = 2.2, adj = 5.9/10)
mtext(bquote(phi)         , side = 2, outer = TRUE, line = 2.2, adj = 4.25/10)
mtext(bquote(b           ), side = 2, outer = TRUE, line = 2.2, adj = 2.4/10)
mtext(bquote(chi         ), side = 2, outer = TRUE, line = 2.2, adj = 0.85/10)

dev.off()

################################################################################

#Fig.4C ########################################################################

if ("ATP" %in% i_reactant) {
  
  png(paste("GM Model ",modelname," Fig 4C (",totalT,"h).png",sep=""),width=4,height=4, units = "in", res=600)
  par(mfrow=c(1,1))
  
  b_ATP <- bt[,i_reactant=="ATP"]
  
  plot(mu_opt,b_ATP, mgp=c(2, 0.5, 0), col.lab = colorlab, ylim=ylimits(b_ATP), type="l",
       xlim=ylimits(mu_opt), xlab=bquote("Instantaneous growth rate" ~ mu(t) ~ (h^-1)),
       ylab=bquote("ATP biomass fraction" ~~ b[ATP](t)), xaxs="i",yaxs="i")
  
  dev.off()  
  
}

setwd(directory)
