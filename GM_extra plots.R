setwd(paste(directory,"/Meta",sep=""))

Afreq <-  read.csv(paste("frequencies ",modelname,".csv",sep=""))

Aphir <-  read.csv(paste("phir ",modelname,".csv",sep=""))

# natural frequencies

a_all     <- Afreq[,2]

Lambda_all    <- Afreq[,3]

Amplitude <- Afreq[,4]

varphi_all   <- Afreq[,5]

bp        <- Afreq[,6]

if (dim(Afreq)[2] == 7)  {
  
  nu <- Afreq[,7] 
  
} else {

freqs     <- Afreq[,7:dim(Afreq)[2]]

nu <- 0
for (i in 1:length(Lambda_all)) nu[i]  <- min(freqs[i,])

}

#Period <- 1/nu

# graphic
coloraxis <- "darkgrey"
colorlab <- "black"

nu0     <- lm(nu ~ Lambda_all)$coefficients[1]

nuslope <- lm(nu ~ Lambda_all)$coefficients[2]

phir    <- Aphir[,4]

muphir  <- Aphir[,3]

colphir <- Aphir[,2]

# mean muphir, phir
mean_muphir <- 0
mean_phir   <- 0
nphir       <- max(colphir)
for (i in 1:nphir) {
  
  mean_muphir[i] <- mean(muphir[colphir==i])
  
  mean_phir[i]   <- mean(phir[colphir==i])
  
}

phi0 <- lm(mean_phir ~ mean_muphir)[[1]][1]

#####################################################################

#pdf(paste("Fig4B ",modelname,".pdf",sep=""),title="Fig4B",width=4,height=4)

png(paste("Fig4B ",modelname,".png",sep=""),width=4,height=4,res=600, units = "in")
par(mfrow=c(1,1))

plot(muphir,phir, mgp=c(2, 0.5, 0),col=c(colphir),pch=20,
     col.lab = colorlab, ylim=c(0,1.02*max(phir)), xlim=c(0,1.05*max(muphir)),
     xlab=bquote("Instantaneous growth rate" ~ mu(t) ~ (h^-1)),ylab=bquote("Ribosome proteome fraction" ~~ phi[r](t)), xaxs="i",yaxs="i")
abline(lm(mean_phir ~ mean_muphir), col=1 , lty=2)
arrows(0.3, phi0, 0, phi0, length = 0.25, angle = 30, code = 1, col = 2, lty = NULL, xpd = FALSE)
text(0.77, phi0, bquote("Ribosome offset" ~~ phi[r0]), col=2)

print( paste("phir0 =",phi0) )

print( paste("phir0 r2 =",cor(mean_phir,mean_muphir)^2 ) )

dev.off()

#####################################################################

#pdf(paste("Fig4A ",modelname,".pdf",sep=""),title="Fig4A",width=4,height=4)

png(paste("Fig4A ",modelname,".png",sep=""),width=4,height=4,res=600, units = "in")
par(mfrow=c(1,1))

nu0 <- lm(nu ~ Lambda_all)[[1]][1]

plot(Lambda_all,nu, mgp=c(2, 0.5, 0),col=1:nphir, col.lab = colorlab, ylim=c(0,1.02*max(nu)), pch=19,
     xlim=c(0,1.05*max(Lambda_all)), xlab=bquote("Long term growth rate" ~ Lambda ~ (h^-1)),
     ylab=bquote("Frequency" ~~ nu ~~ (h^-1)), xaxs="i",yaxs="i")
abline(lm(nu ~ Lambda_all), col=1 , lty=2)
arrows(0.3, nu0, 0, nu0, length = 0.25, angle = 30, code = 1, col = 2, lty = NULL, xpd = FALSE)
text(0.75, nu0, bquote("Frequency offset" ~ nu[0]), col=2)

print(paste("nu0 =",nu0) )

print(paste("nu r2 =",cor(nu,Lambda_all)^2) )

dev.off()  

#####################################################################

#pdf(paste("FigS21 ",modelname,".pdf",sep=""),title="FigS21",width=4,height=4)

png(paste("FigS21 ",modelname,".png",sep=""),width=4,height=4,res=600, units = "in")
par(mfrow=c(1,1))

bp0 <- lm(bp  ~ Lambda_all)[[1]][1]

plot(Lambda_all,bp , mgp=c(2, 0.5, 0),col=1:nphir, col.lab = colorlab, ylim=c(0,1.1*max(bp)), pch=19,
     xlim=c(0,1.05*max(Lambda_all)), xlab=bquote("Long term growth rate" ~ Lambda ~ (h^-1)),
     ylab=bquote("Proteome biomass fraction" ~~ b[p]), xaxs="i",yaxs="i")
abline(lm(bp ~ Lambda_all), col=1 , lty=2)

print(paste("bp0 =",bp0) )

print(paste("bp r2 =",cor(bp,Lambda_all)^2) )

dev.off()  

#####################################################################

png(paste("Fig Lambda x varphi ",modelname,".png",sep=""),width=4,height=4,res=600, units = "in")
par(mfrow=c(1,1))

#bp0 <- lm(Lambda_all  ~ Lambda_all)[[1]][1]

plot(varphi_all , Lambda_all,  mgp=c(2, 0.5, 0),col=1:nphir, col.lab = colorlab, xlim=c(0,1.1*max(varphi_all)), pch=19,
     ylim=c(0,1.05*max(Lambda_all)), ylab=bquote("Long term growth rate" ~~ Lambda ~ (h^-1)),
     xlab=bquote("Scalar potential" ~~ varphi ~ (h^-1)), xaxs="i",yaxs="i")
abline(lm(varphi_all ~ Lambda_all), col=1 , lty=2)

#print(paste("bp0 =",bp0) )

print(paste("r2 =",cor(,varphi_all, Lambda_all)^2) )

dev.off()  

# h
h_all <- (1 - varphi_all/Lambda_all)*log(2)

plot(Lambda_all, h_all,  mgp=c(2, 0.5, 0),col=1:nphir, col.lab = colorlab, xlim=c(0,1.1*max(Lambda_all)), pch=19,
     ylim=c(0,1.05*max(h_all)), ylab=bquote("h"),
     xlab=bquote("Long term growth rate" ~~ Lambda ~ (h^-1)), xaxs="i",yaxs="i")
abline(lm(varphi_all ~ e_all), col=1 , lty=2)

setwd(directory)

