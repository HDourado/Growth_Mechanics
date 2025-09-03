# Interface

# Clear variables
rm(list=ls(all=TRUE))

require('rstudioapi') 
require('readODS')
require('nloptr')
require('Matrix')
require('MASS')
require('lpSolve')
require("latex2exp")
require("lattice")
require("deSolve")
require("rgl")
require("plot3D")
require('gplots')
require('Rlibeemd')
require('sfsmisc')

options(digits=10)

# Sets working directory to source file location ###############################

directory <- dirname(getActiveDocumentContext()$path)

setwd(directory) 

# Growth Balance Analysis (GBA) function for balanced growth simulations #######

GBA <- function(modelname,predict.parameters,is.reversible) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible")))
  
  modelname <<- modelname
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  suppressMessages(source("GBA.R"))
  
}

# Example:
# GBA("L3",3,0)

# Growth Mechanics (GM) function for dynamic simulations #######################

# modelname: name of a document in the Models file in quotes, e.g. "L3"

# predict.parameters: 0 = does not predict, take from file, 1,2,3,... = predicts using this parameter as K_ratio, as described in "methods".

# is.reversible (only relevant if what to predict kinetic parameters): 0 = not , 1 = yes.

# delta = interval of time in hours

# totalT: total time in hours

# Example:
# GM("L3",3,0,0.001,3)

# for accurate average mu calculations, use smaller dt, e.g. dt = 0.00001
# GM("L3",3,0,0.00001,3)

GM <- function(modelname,predict.parameters,is.reversible,delta_t,totalT) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                  "delta_t","totalT")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters if > 0, using cm/Km = predict.parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  rtol <<- 1e-8
  
  options(digits=10)
  
  suppressMessages(source("GM.R"))
  
}

# Example:
# GM("L3",3,0,0.001,3)

# For more accurate (and slower) average mu calculations, use smaller dt, e.g. dt = 0.00001
# GM("L3",3,0,0.00001,3)

# Function to do many GM simulations for a model in different static media
GMmeta <- function(modelname,predict.parameters,is.reversible,delta_t,totalT,nmeta) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                          "delta_t","totalT","nmeta")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  ameta1 <<- c(10, 0.25, 0.1, 0.05, 0.03, 0.018, 0.01)
  
  #ameta1 <<- c(10, 9, 0.25, 0.245, 0.1, 0.099, 0.05, 0.049,0.025, 0.024, 0.015, 0.0145)
  
  for (itest in 1:length(ameta1)) {
    
    # number of external medium concentrations
    nmc <- length(at(0))
    
    if ( nmc == 1 )  {
      
      at <<- function(t) ameta1[itest]
      
    } else at <<- function(t) c(ameta1[itest], rep(10, nmc - 1 ) )
    
    rtol <<- 1e-8
    
    print("----------------------")
    
    print(paste("Optimization",itest,"of",length(ameta1)))
    
    source("GM.R")
    
    # # do it again if negative chi
    # while (min(chit[-c(1:5),]) < 0) { # ignores first 5 points, due to possible numerical artifacts
    # 
    #   dq0[r] <<- dq0[r]*0.9
    # 
    #   source("GM.R")
    # 
    # }
    
    dq0[r] <<- dq0[r]*0.9
    
    # results to save (x,freq, mu, amplitude)
    if(itest == 1) results_freq <- c(ameta1[itest], Lambda,max(mu_opt) - min(mu_opt), varphi0, mean(bt[,p]), freqs )
    
    if(itest > 1)  results_freq <- rbind(results_freq, c(ameta1[itest], Lambda,max(mu_opt) - min(mu_opt), varphi0, mean(bt[,p]), freqs) )
    
    if(itest == 1) results_phir <- cbind(itest,mu_opt,phit[,r])
    
    if(itest > 1) results_phir  <- rbind(results_phir, cbind(itest,mu_opt,phit[,r]))
    
  }
    
    setwd(paste(directory,"/Meta",sep=""))
    
    # export results
    write.csv(results_freq, file = paste("frequencies ",modelname,".csv",sep=""))
    
    write.csv(results_phir, file = paste("phir ",modelname,".csv",sep=""))
    
    setwd(directory)
    
    print("----- END ------")
    
    source("GM_extra plots.R")
  
}

# Example:
# GMmeta("L3",3,0,0.001,10)

# Function to do many GM simulations for a model in a static media with different initial conditions
GMinitial <- function(modelname,predict.parameters,is.reversible,delta_t,totalT,n_initial) {
  
  # first deletes all variables from previous calculations 
  rm(list=setdiff(ls(), c("modelname","predict.parameters", "is.reversible",
                          "delta_t","totalT","n_initial")))
  
  modelname <<- modelname
  
  # Forces all reactions to be irreversible if desired
  is.reversible <<- is.reversible
  
  # Predicts kinetic parameters
  predict.parameters <<- predict.parameters
  
  # delta t
  delta_t <<- delta_t
  
  # Total simulation time
  totalT <<- totalT
  
  # reads model
  suppressMessages(source("Readmodelods.R"))
  
  delta_t <<- 0.00001
  
  dqr0 <<- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
  
  for (itest in 1:length(dqr0)) {
    
    rtol <<- 1e-8
    
    print("----------------------")
    
    print(paste("Optimization",itest,"of",length(dqr0)))
    
    dq0[r] <<- dqr0[itest] 
    
    source("GM.R")
    
    # results to save (at,freq, mu, amplitude)
    if(itest == 1) results_freq <- c(dqr0[itest], Lambda , max(mu_opt) - min(mu_opt), varphi0, mean(bt[,p]), nu )
    
    if(itest > 1)  results_freq <- rbind(results_freq, c(dqr0[itest], Lambda, max(mu_opt) - min(mu_opt), varphi0, mean(bt[,p]), nu) )
    
  }
  
  setwd(paste(directory,"/Meta",sep=""))
  
  # export results
  write.csv(results_freq, file = paste("initial ",modelname,".csv",sep=""))
  
  png(paste("Initial ",modelname,".png",sep=""),res=600, units = "in", width=5,height=5)
  
  par(mfrow=c(1,1))
  
  plot(amplitude_i,100*gap,pch=16, xlab=bquote("Amplitude of" ~  mu(t) ), ylab=bquote("Difference between " ~  Lambda ~ "and" ~ mu ~ "* (%)" ) )
  
  dev.off()
  
  setwd(directory)
  
  print("----- END ------")
  
}

# Example:
# GMinitial("L3",3,0,0.00001,3)

dev.off()

