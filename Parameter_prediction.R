# Kinetic parameter prediction

mu_data  <- exp_data[1]

phi_data <- exp_data[-1]

cm0 <-  as.numeric(rho*mu_data*M%*%q0_sgm)[-p]

esat <- predict.parameters

# Estimates Km values assuming typical ratio Km = cm/esat for substrates and 
# Km = cm*esat for products

# forces some reactions to be irreversible if kcatb = 0 in the file
rev <- as.vector(1*(kcatb > 0))

K[(n_a+1):(p+n_a-1),] <- ceiling( diag(cm0/esat)%*%(1*(M[-p,] < 0)) + is.reversible*diag(cm0)%*%(1*(M[-p,] > 0))%*%diag(rev) )

# Estimates total saturation of each reaction assuming irrev. with cm = esat*Km

sat <- (esat/(esat + 1))^(colSums(1*(M[-p,] < 0)))

# assumes transporters are completely saturated

sat[1:n_tr] <- rep(1,n_tr)

b0 <- M%*%q0_sgm

kcatf <- t(ceiling(mu_data*q0_sgm/(b0[p]*phi_data*sat)))

# if some kcat was rounded to zero
kcatf[kcatf == 0] <- 1

# estimating kcatb
kcatb <- round(kcatf/5)

kcatb[kcatb == 0] <- 1

kcatb <- is.reversible*rev*kcatb

# ribosome is irreversible
kcatb[r] <- 0

# first separate Km matrices for substrates and for products
KS <- K*(Mtotal<0)
KP <- K*(Mtotal>0)

KS[KS == 0] <- Inf

KP[KP == 0] <- Inf

KI[KI == 0] <- Inf

