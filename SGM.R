# Singular Growth Modes ########################################################

rM <- rankMatrix(M)[1]

# Singular Value Decomposition of M
U   <- svd(M,nu=p,nv=r)$u 

V <- svd(M,nu=p,nv=r)$v 

U[abs(U) < 1e-10] <- 0

V[abs(V) < 1e-10] <- 0

sigma <- svd(M)$d[1:rM]

Sigma <- matrix(rep(0,p*r),ncol=r)

Sigma[1:rM,1:rM] <- diag(sigma)

correction <- sign(colSums(U%*%Sigma)) #sign(V[r,1:rM])

for (i in 1:rM) U[,i] <- U[,i]*correction[i]

for (i in 1:rM) V[,i] <- V[,i]*correction[i]

################################################################################

# definition of mode values ####################################################

alpha <- 1:r
beta  <- 1:p

# test whether there is SGM with all non-negative c

signU <- 1*(U > 0)

signV <- 1*(V > 0)

is.BGS <- as.numeric(colSums(signU) == p)

pBGS <- grep(p,colSums(signU))

# estimate f0 based on SGM
# if the is a BGM (SGM with non-negative concentrations)
if (length(pBGS) == 1) {

  q0_sgm <- V[,pBGS]/(sum(U[,pBGS])*sigma[pBGS])

  b0_sgm <- M%*%q0_sgm
  
  z <- V%*%q0_sgm
  
  #pie(q0_sgm,reaction)
  
  #pie(b0_sgm,i_reactant)

} 

q0 <- q0_sgm
