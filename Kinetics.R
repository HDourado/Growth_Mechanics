# Turnover time function for different kinetic rate laws, and their partial
# derivatives with respect to b and a 

# Reversible MM  ######################################################################################

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

KS[KS == 0] <- Inf

KP[KP == 0] <- Inf

KI[KI == 0] <- Inf

# Global kinetic function ######################################################

tau <- function(t,cint) {
  
  ca <- c(at(t),cint)
  
  tauj <- rep(0,r)
  
  for (j in 1:r) {
    
    subr <- ca/KS[,j]
    
    subr <- subr[subr != 0]
    
    prodr <- ca/KP[,j]
    
    prodr <- prodr[prodr != 0]
    
    tauj[j] <- as.numeric( prod(1 + KA[,j]/ca)*prod(1 + ca/KI[,j])*( prod(1 + subr) + prod(1 + prodr) - 1 )
      
      /( kcatf[j]*prod(subr) - kcatb[j]*prod(prodr)  ) 
      
    )
    
  }
  
  tauj
  
}

# Global kinetic derivatives function ##########################################

dtau <- function(t,cint) {
  
  ca <- c(at(t),cint)
  
  tauj <- rep(0,r)
  
  dtauj <- matrix(rep(0,r*p),ncol=p)
  
  for (j in 1:r) {
    
    subr <- ca/KS[,j]
    
    subr <- subr[subr != 0]
    
    prodr <- ca/KP[,j]
    
    prodr <- prodr[prodr != 0]
    
    # tau = tauA*tauB/tauC
    
    # tauA
    taujA <- prod(1 + KA[,j]/ca)*prod(1 + ca/KI[,j])
    
    taujB <- prod(1 + subr) + prod(1 + prodr) - 1
    
    taujC <- kcatf[j]*prod(subr) - kcatb[j]*prod(prodr)
    
    # calculates two terms that depend only on a
    dtaujA <- rep(0,p)
    dtaujB <- rep(0,p)
    dtaujC <- rep(0,p)
    
    for (i in 1:p) {
      
      # y = i position counting for environment reactants n
      
      y <- i + n_a
      
      subr <- ca[-y]/KS[-y,j]
      
      subr0 <- subr[subr != 0]
      
      subr1 <- subr0 + 1
      
      prodr <- ca[-y]/KP[-y,j]
      
      prodr0 <- prodr[prodr != 0]
      
      prodr1 <- prodr0 + 1
      
      # dtau A
     
      dtaujA[i] <-  -(KA[y,j]/(ca[y]^2))*prod(1 + KA[-y,j]/ca[-y])*prod(1 + ca/KI[,j]) + 
        
                     prod(1 + KA[,j]/ca)*(1/KI[y,j])*prod(1 + ca[-y]/KI[-y,j])
      
      
      dtaujB[i] <-  prod(subr1)/KS[y,j] + prod(prodr1)/KP[y,j]
      
      
      dtaujC[i] <-  (kcatf[j]/KS[y,j])*prod(subr0) - (kcatb[j]/KP[y,j])*prod(prodr0) 
      
    }
    
    dtauj[j,] <- as.numeric(( dtaujA*taujB + taujA*dtaujB - taujA*taujB*dtaujC/taujC  )/taujC)
    
  }
  
  dtauj
  
}

# Derivatives wrt medium concentrations a ######################################

dtauda  <- function(t,cint) {
  
  ca <- c(at(t),cint)
  
  tauj <- rep(0,r)
  
  dtauj <- matrix(rep(0,r*n_a),ncol=n_a)
  
  for (j in 1:r) {
    
    subr <- ca/KS[,j]
    
    subr <- subr[subr != 0]
    
    prodr <- ca/KP[,j]
    
    prodr <- prodr[prodr != 0]
    
    # calculates tauj first
    
    tauj <- as.numeric( ( prod(1 + subr) + prod(1 + prodr) - 1 )
                        
                        /( kcatf[j]*prod(subr) - kcatb[j]*prod(prodr)  ) )
    
    
    # calculates two terms that depend only on a
    term1 <- rep(0,n_a)
    term2 <- rep(0,n_a)
    
    for (n in 1:n_a) {
      
      subr <- ca[-n]/KS[-n,j]
      
      subr0 <- subr[subr != 0]
      
      subr1 <- subr0 + 1
      
      prodr <- ca[-n]/KP[-n,j]
      
      prodr0 <- prodr[prodr != 0]
      
      prodr1 <- prodr0 + 1
      
      
      term1[n] <-  prod(subr1)/KS[n,j] + prod(prodr1)/KP[n,j]
      
      term2[n] <-  (kcatf[j]/KS[n,j])*prod(subr0) - (kcatb[j]/KP[n,j])*prod(prodr0)                                              
      
    }
    
    subr  <- ca/KS[,j]
    
    subr0 <- subr[subr != 0]
    
    subr1 <- subr0 + 1
    
    prodr <- ca/KP[,j]
    
    prodr0 <- prodr[prodr != 0]
    
    prodr1 <- prodr0 + 1
    
    
    dtauj[j,] <- tauj*as.numeric( term1/(prod(subr1) + prod(prodr1) - 1)
                                  
                                  - term2/ (kcatf[j]*prod(subr0) - kcatb[j]*prod(prodr0) ) )
    
  }
  
  dtauj
  
}

