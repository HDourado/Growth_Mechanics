# exporting file with the optimal states ci, tau, mu, v, p 

setwd(paste(directory,"/Results GM",sep=""))

simulation_time <- signif(st[[1]], digits= 4)

opt_state <- matrix(0,nt,(n_a+3+p+r+r+r+r+r))
for (cond in 1:nt) {
  
  t <- delta_t*(cond - 1)
  
  a <- at(t)
  
  da <- dadt(t)
  
  qj <- qt[cond,]
  
  dqj <- dqt[cond,]
  
  opt_state[cond,] <- c(conv[cond],t,mu_opt[cond],a,b(qj),tau(t,rho*b(qj)),v(qj,dqj,t),prot(qj,dqj,t),qj,dqj)
}

colnames(opt_state) <- c("convergence","time","mu",reactant,paste("tau",reaction),
                         paste("v",reaction),paste("p",reaction),paste("q",reaction),paste("dq",reaction))

# export results
write.csv(opt_state, file = paste("GM Model ",modelname,
                            " (",simulation_time,"s).csv",sep=""))

setwd(directory)