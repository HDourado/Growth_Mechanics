# exporting file with the optimal states ci, tau, mu, v, p 

setwd(paste(directory,"/Results GBA",sep=""))

opt_state <- matrix(rep(0,(n_a+2+p+r+r+r+r)*n_conditions),nrow = n_conditions)
for (cond in 1:n_conditions) {
  
  t <- cond - 1
  
  a  <- at(t)
  
  q <- q_opt[cond,]
  
  opt_state[cond,] <- c(conv[cond],mu(q),a,ci(q),tau(t,ci(q)),v(q),prot(q),q)
}

colnames(opt_state) <- c("convergence","mu",reactant,paste("tau",reaction),
                         paste("v",reaction),paste("p",reaction),paste("q",reaction))

# export results
write.csv(opt_state, file = paste("GBA Model ",modelname,
                            ", mean time (",mean_time,"s) results.csv",sep=""))


setwd(directory)