# Functions ####################################################################

# mu ############## (growth rate)
mu <- function(q) as.numeric(M[p,r]*q[r]/(tau(t,rho*b(q))%*%q)) 

# negative mu (for minimization)
negative_mu <- function(q) -as.numeric(M[p,r]*q[r]/(tau(t,rho*b(q))%*%q)) 

# fluxes
v <- function(q) as.vector(mu(q))*rho*q

# protein concentrations
prot <- function(q) tau(t,rho*b(q))*v(q)

# internal concentrations "i" of metabolites "m" and proteome "p"
ci <- function(q) rho*M%*%q

# biomass fractions
b <- function(q) M%*%q

# proteome fractions
phi <- function(q) (tau(t,rho*b(q))*v(q))/(ci(q)[p])

# Now special functions only if there is a "rRNA" in the model, for ribosome composition
if ("rRNA" %in% i_reactant) rprna <- function(q) {
  
  rRNA <- ci(q)[grep("rRNA",i_reactant)] 

  # rRNA bound to rProtein
  brRNA <- rRNA*(rRNA/(rRNA + KA[grep("rRNA",reactant),r]) )
  
  # ribosome protein concentration
  rp <- prot(q)[r]
  
  # ratio between ribosomal protein and ribosomal rna (rP/(rP + brRNA))
  return( rp/(rp + rRNA) ) 

}

if ("rRNase" %in% reaction) {
  # extra constraint on rRNA degradation
  j_rRNAse <- grep("rRNase", reaction) # column on this reaction, input
  i_rRNA   <- grep("rRNA", i_reactant) # row of this 
  Km       <- K[i_rRNA + n_a,j_rRNAse]  # Km of this reaction
  KI_rP    <- 100  # "inhibition constant" 
  v_rRNAse <- function(q) kcatf[j_rRNAse]*prot(q)[j_rRNAse]*( ci(q)[i_rRNA]/(Km + ci(q)[i_rRNA] ) )*( KI_rP/(KI_rP + phi(q)[r]*ci(q)[p])  )
  
  extra_ineq <- function(q) v_rRNAse(q) - v(q)[j_rRNAse]
  
}

# Optimization #################################################################

# equality constraints (density constraint)
g <- function(q) s%*%q - 1

# equality constraints if also constraint on ribosome composition

if (ribcomp > 0 & "rRNA" %in% i_reactant) g <- function(q) c(s%*%q - 1, rprna(q) - ribcomp )
  
# inequality constraints (min c, max c, min phi, max phi)

h <- function(q) c( ci(q) - min_c, max_c - ci(q), phi(q) - min_phi, max_phi - phi(q))

# Indirect elasticities ########################################################

E <- function(q) rho*dtau(t,rho*b(q))%*%M

# Gradient of negative_mu ######################################################

negative_dmu <- function(q) -((mu(q)^2)/b(q)[p])*(M[p,]/mu(q) - 
                                      t(q)%*%(rho*dtau(t,rho*b(q))%*%M) - tau(t,rho*b(q)))

# starts loop for the optimization at each environmental condition #############
q_opt  <- matrix(rep(0,r*n_conditions),ncol=r)
mu_opt <- rep(0,n_conditions)
otime  <- rep(0,n_conditions)
conv   <- rep(0,n_conditions)
iter   <- rep(0,n_conditions)
solver_cond <- rep(0,n_conditions)
dmu_opt  <- matrix(rep(0,r*n_conditions),ncol=r)

t <- 0

suppressMessages(
for (cond in 1:n_conditions) {
  
  t <- cond - 1
  
  a  <- at(t)
  
  # Upper bounds 
  upper_q <- max_q
  
  # lower bounds 
  lower_q <- min_q
  
  # Optimization
  solver      <- "SLSQP" # Local solvers: "COBYLA", "LBFGS", "MMA", or "SLSQP"
  
  # measuring the total optimization time
  st <- system.time({
    
    # optimization without gradients 
    res <- auglag(q0, fn = negative_mu, lower = lower_q, 
                  upper = upper_q, heq = g, hin = h, 
                  localsolver = c(solver), localtol = 1e-8,
                  control = list(maxeval = 100000))
  })
  
  # solution
  q0 <- res$par
  
  # solution
  q_opt[cond,] <- q0
  
  # optimal mu
  mu_opt[cond] <- mu(q0)
  
  # convergence (codes: "-1" means optimization problem, "5" means optimization 
  # stop because maxeval, "4" means optimization stopped because xtol_rel or 
  # xtol_abs (above) was reached)
  conv[cond] <- res$convergence
  
  # optimization time
  otime[cond] <- signif(st[[1]], digits= 4)
  
  # number of iterations
  iter[cond] <- res$iter  
  
  print(paste("optimization: ",cond, "/", n_conditions,", optimization time: ",
              otime[cond]," s, convergence: ", conv[cond],", growth rate: ",
              signif(mu_opt[cond], digits= 3), sep=""))
  }
)

print("----------------- Finished --------------------")

# produces c_opt
c_opt <- matrix(rep(0,p*n_conditions),ncol=p)
for (cond in 1:n_conditions)  c_opt[cond,] <- ci(q_opt[cond,])

# mean optimization time
mean_time <- signif(mean(otime), digits= 3)
