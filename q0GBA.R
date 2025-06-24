# Growth optimization

#source("GBA_Kinetics.R")

# Functions #################################################################### 

# mu ############## (growth rate)
mu <- function(q) as.numeric(M[p,r]*q[r]/(tau(0,rho*b(q))%*%q)) 

# negative mu (for minimization)
negative_mu <- function(q) -as.numeric(M[p,r]*q[r]/(tau(0,rho*b(q))%*%q)) 

# biomass fractions
b <- function(q) M%*%q

# Optimization #################################################################

# equality constraints (density constraint)

g <- function(q) s%*%q - 1

# inequality constraints (non negativity of concentrations). Here we assume ####
# mu ~ 1 for a faster calculation of p, instead of using the function p(f) #####

h <- function(q) c( rho*b(q), rho*tau(0,rho*b(q))*q )

# Gradient of negative_mu ######################################################

negative_dmu <- function(q) -((mu(q)^2)/b(q)[p])*(M[p,]/mu(q) - t(q)%*%(rho*dtau(0,rho*b(q))%*%M) - tau(0,rho*b(q)))

# starts loop for the optimization at each environmental condition #############

maxbound <- 10*max(q0_sgm)

x  <- at(0) 

# Upper bounds 
upper_q <- rep(maxbound,r)

# lower bounds 
lower_q <- rep(-maxbound,r)

# Optimization
solver      <- "SLSQP" # Local solvers: "COBYLA", "LBFGS", "MMA", or "SLSQP"

# optimization (by default the function auglag minimizes, so we use 
# negative_mu instead of mu as the objective function)

res <- auglag(q0_sgm, fn = negative_mu, gr = negative_dmu, lower = lower_q, 
                                upper = upper_q, heq = g, hin = h, 
                                localsolver = c(solver), localtol = 1e-8,
                                control = list(maxeval = 100000)) 
q0 <- res$par
