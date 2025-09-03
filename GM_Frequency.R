# Frequency and dynamical study

setwd(paste(directory,"/Results GM",sep=""))

# Calculates frequencies

# extrema of mu
extmu <- extrema(mu_opt)

if (length(extmu$maxima[,1]) > 3) {
  
  mumaxt <- extmu$maxima[,1] 
  
  mumaxt <- mumaxt[-c(1,length(mumaxt))]
  
  mu_T <<- mean(diff(mumaxt)*delta_t)
  
  mu_f <<- 1/mu_T
  
  sd_mu_f <- sd(1/diff(mumaxt)*delta_t)
  
  # phir
  extphir <- extrema(phit[,r])
  
  phirmaxt <- extphir$maxima[-c(1:round(length(extphir[[1]][,1])/2)),1]
  
  phirmaxt <- phirmaxt[-length(phirmaxt)]
  
  phir_T <- mean(diff(phirmaxt)*delta_t)
  
  phir_f <- 1/phir_T
  
  sd_phir_f <- sd(1/diff(phirmaxt)*delta_t)
  
} else  {
  
  mu_T <<- abs( extmu$minima[2,1] - extmu$maxima[2,1] )*delta_t*2
  
  mu_f <<- 1/mu_T
  
  sd_mu_f <- 0
}

# mu FFT #######################################################################

mufft <- fft(mu_opt- mean(mu_opt))

magnitudes <- abs(mufft) 

max_index_sine <- which.max(magnitudes[1:round(nt/2)]) - 1

maxf <- (max_index_sine)/nt/delta_t

frequencies <- (1:nt -1)/nt/delta_t

singsdiffmag <- sign(diff(magnitudes))

# finds maxima in magnitudes
maxifreq <- 0
for (i in 1:(max_index_sine*20)) maxifreq[i] <- frequencies[i+1]*(singsdiffmag[i] - singsdiffmag[i+1])/2

maximag    <- magnitudes[2:(max_index_sine*20 +1)]

tablefm <- cbind(maxifreq,maximag)

tablefm <- tablefm[maxifreq > 0,]

if (length(tablefm) > 2) tablefm <- tablefm[order(tablefm[,2],decreasing=T),]

if (length(tablefm) > 2) tablefm <- tablefm[tablefm[,2] > max(magnitudes)/5 | (tablefm[,1] < maxf & tablefm[,2] > max(magnitudes)/25), ]

if (length(dim(tablefm))==0 ) {
  
  nlines <-  1
  
  freqs <- tablefm[1]
  
} else {
  
  nlines <-  dim(tablefm)[1]
  
  freqs <- tablefm[1:2,1]
  
}

scaledmag <- magnitudes[-1]/max(magnitudes[-1]) 

freq_legend <- paste(signif(freqs,digits=3),"/h",sep="")

if (length(freq_legend) > 6) freq_legend <- freq_legend[1:6]

################################################################################

#pdf(paste("GM Model ",modelname," frequency (",totalT,"h).pdf",sep=""),title="Frequency study",
#    width=5,height=3)

png(paste("GM Model ",modelname," frequency (",totalT,"h).png",sep=""),res=600, units = "in", width=5,height=3)

par(mfrow=c(1,1))

plot(frequencies[-1],scaledmag, type = "l",xlim=c(0,2*max(freqs)), mgp=c(2, 0.5, 0), col.lab = colorlab,
     xlab="Frequency (1/h)",ylab="Amplitude [a.u.]", xaxs="i")
legend("topright",freq_legend,col=1:nlines +1,pch=16,cex=0.8)
for (i in 1:nlines) abline(v=freqs[i],col=i+1,lty=2)

dev.off()

################################################################################

png(paste("GM Model ",modelname," dynamics q (",totalT,"h).png",sep=""),res=600, units = "in", width=5,height=5)
par(mfrow=c(1,1))

qrm1 <- paste("q_",r-1,sep="")

lines3D(qt[,r], qt[,r-1],mu_opt, colvar=time, xlab=paste("q_",r-1,sep=""),ylab=paste("q_",r-2,sep=""),
        zlab="mu (1/h)")

dev.off()

################################################################################

# all fft #########################################

xfft <- function(something) {
  
mufft <- fft(something- mean(something))

magnitudes <- abs(mufft) 

max_index_sine <- which.max(magnitudes[1:round(nt/2)]) - 1

maxf <- (max_index_sine)/nt/delta_t

frequencies <- (1:nt -1)/nt/delta_t

singsdiffmag <- sign(diff(magnitudes))

# finds maxima in magnitudes
maxifreq <- 0
for (i in 1:(max_index_sine*20)) maxifreq[i] <- frequencies[i+1]*(singsdiffmag[i] - singsdiffmag[i+1])/2

maximag    <- magnitudes[2:(max_index_sine*20 +1)]

tablefm <- cbind(maxifreq,maximag)

tablefm <- tablefm[maxifreq > 0,]

if (length(tablefm) > 2) tablefm <- tablefm[order(tablefm[,2],decreasing=T),]

if (length(tablefm) > 2) tablefm <- tablefm[tablefm[,2] > max(magnitudes)/5 | (tablefm[,1] < maxf & tablefm[,2] > max(magnitudes)/25), ]

if (length(dim(tablefm))==0 ) {
  
  nlines <-  1
  
  freqs <- tablefm[1]
  
} else {
  
  nlines <-  dim(tablefm)[1]
  
  freqs <- tablefm[,1]
  
}

scaledmag <- magnitudes[-1]/max(magnitudes[-1]) 

# takes only first 5 first frequencies in legend

freq_legend <- paste(signif(freqs,digits=3),"/h",sep="")

if (length(freq_legend) > 6) freq_legend <- freq_legend[1:6]

}

# take the main frequency as determined by a FFT of mu(t), if n >=5
if (r > 4) {
  
  mu_f <<- freqs

}

if (mu_f > 10*varphi0) mu_f <- freqs[1]

mu_T <-1/mu_f

# Make pie/clock chart with order which biomass allocation peaks
period <- round(mu_T/delta_t)

bt  <- cit/rho

bpeaks  <- matrix(rep(0,p*2),ncol=2)

phipeaks <- matrix(rep(0,r*2),ncol=2)

pbpeaks <- phipeaks

tn0 <- grep(max(pt[1:period,1]), pt[1:period,1]) #- 1

# define one cycle based on first transporter
tinterval <- c( tn0:(tn0 + period) )   # c( extrema(pt[,1])$maxima[3,1]:(extrema(pt[,1])$maxima[4,1] - 1) )

# it may not oscillate first, depending on the medium, so we need to check this

if (period < nt) {

# defines one cycle based on mu_T

for (i in 1:p) bpeaks[i,]   <- c( grep(max(bt[tinterval,i]),bt[tinterval,i])[1] , max(bt[tinterval,i]) )

for (j in 1:r) phipeaks[j,] <- c( grep(max(phit[tinterval,j]),phit[tinterval,j])[1] , max(phit[tinterval,j]) ) 

for (j in 1:r) pbpeaks[j,]  <- c( grep(max(pt[tinterval,j]),pt[tinterval,j])[1], max(pt[tinterval,j])/rho ) 

peaks <- rbind(bpeaks[-p,],pbpeaks)

peaks[,1] <- peaks[,1]*delta_t # converts to hour

#pdf(paste("GM Model ",modelname," peak times (",totalT,"h).pdf",sep=""),title="Peak allocation times",
#    width=7,height=3)

png(paste("GM Model ",modelname," peak times (",totalT,"h).png",sep=""),res=600, units = "in", width=7,height=3)
par(mfrow=c(1,1))

plot(peaks,col=c(rep(8,p-1),colorj ) , pch=c(rep(19,p-1), rep(15,r)), ylim=c(0,1.4*max(peaks[,2])),
     ylab="biomass fraction", xlab="peak allocation time (h)",xlim=c(0,mu_T))
text(peaks,labels=c(metabolite,reaction), pos=3 , cex=0.6)

dev.off()

}

setwd(directory)
  
  