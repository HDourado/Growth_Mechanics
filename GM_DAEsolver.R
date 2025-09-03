# Growth optimization

# Functions ####################################################################

# mu ############## (growth rate)
mu <- function(q,dq,t) as.numeric( (M[p,r]*q[r] - tau(t,rho*b(q))%*%dq )/(tau(t,rho*b(q))%*%q)) 

# fluxes mu*rho*f
v <- function(q,dq,t) as.vector(rho*(mu(q,dq,t)*q + dq))

# "distance" R to q0_OGS #######################################################
R <- function(q,q0_OGS) sqrt(sum( (q[-1] - q0_OGS[-1] )^2))

# protein concentrations
prot <- function(q,dq,t) tau(t,rho*b(q))*v(q,dq,t)

# proteome fractions
phi <- function(q,dq,t) (tau(t,rho*b(q))*v(q,dq,t))/(rho*b(q)[p])

# biomass fractions
b <- function(q) M%*%q

# indirect elasticity matrix E 
dtaudq <- function(q,t) rho*dtau(t,rho*b(q))%*%M

# dtau/dt_j
dtaudtj <- function(q,dq,t,j) as.numeric(dq%*%(t(dtaudq(q,t))[,j]) + t(dadt(t))%*%(t( dtauda(t,rho*b(q)) )[,j]) )

# dtau/dt
dtaudt <- function(q,dq,t) {
  
  tj <- 1:r
  
  for (j in 1:r) tj[j] <- dtaudtj(q,dq,t,j) 
    
return(tj)
  
}

# "electric" field = grad varphi (valid only for 1 transporter, to be generalized)
E_i <- function(q,t) {
  
  grad <- rep(0,r-1)
  
  for (j in 1:r) {
  
  grad[j] <- M[p,j] - mu(q,0*q,t)*tau(t,rho*b(q))[j] - t(mu(q,0*q,t)*q)%*%(  
    
    c(dtaudq(q,t)[,j] + dtaudq(q,t)%*%q*s[j]) )
  
  }
  
  grad <- - grad/as.numeric(tau(t,rho*b(q))%*%q)
  
  return(grad[-1] )
  
}

# Energy ############## 
varphi <- function(q,t) as.numeric( (M[p,r]*q[r] )/(tau(t,rho*b(q))%*%q)) 

# A_i ############## ("magnetic" potential) 
A_i <- function(q,t) tau(t,rho*b(q))[-1]/as.numeric(tau(t,rho*b(q))[1]/s[1] +  t(tau(t,rho*b(q))[-1])%*%q[-1] )  # (varphi(q,t)/M[p,r])*tau(t,rho*b(q))[-1]/q[r] 

# starts loop for the optimization at each environmental condition #############

# number of time points
nt <- (totalT/delta_t) + 1

qt     <- matrix(rep(0,r*nt),ncol=r)
dqt    <- matrix(rep(0,r*nt),ncol=r)
kappat     <- matrix(rep(0,r*nt),ncol=r)
mu_av  <- 0
otime  <- rep(0,nt)
conv   <- rep(0,nt)
iter   <- rep(0,nt)

# EOM for reaction j (= 0 at optimal state)
motionj <- function(q,dq,t,j) M[p,j] - mu(q,dq,t)*tau(t,rho*b(q))[j] - t(dq + mu(q,dq,t)*q)%*%(  
  
           c(dtaudq(q,t)[,j] + dtaudq(q,t)%*%q*s[j]) ) + dtaudtj(q,dq,t,j) - 
  
           (tau(t,rho*b(q))[j]/(t(q)%*%tau(t,rho*b(q))))*( t(dq)%*%tau(t,rho*b(q)) + 
                                                       
                                                       t(q)%*%dtaudt(q,dq,t) )  

# system of algebraic differential equations (equations of motion + density constraint)
Res_DAE <- function(t,q,dq, pars) 
{
  with (as.list(c(q, dq, pars)), {
    
    a <- at(t)
    
    da <- dadt(t)
    
    res <- numeric(r)
    
    # factor (1000) multiplying density constraint has to be the same as 1/rtol!
    
    for (j in 1:r) {
      
      res[j] <- M[p,j] - mu(q,dq,t)*tau(t,rho*b(q))[j] - t(dq + mu(q,dq,t)*q)%*%(  
      
      c(dtaudq(q,t)[,j] + dtaudq(q,t)%*%q*s[j]) ) + dtaudtj(q,dq,t,j) - 
      
      (tau(t,rho*b(q))[j]/(t(q)%*%tau(t,rho*b(q))))*( t(dq)%*%tau(t,rho*b(q)) + 
                                                  
                        t(q)%*%dtaudt(q,dq,t) )  
    }
    
    res[1:n_tr] <- res[1:n_tr] - rep((s%*%q - 1)/rtol,n_tr) # adds continuity constraint scaled by rtol
     
    return(list(res))
  })
}

t       <- 0

dqt[1,] <- dq0

names(q0) <- paste("q0_",1:r,sep="")

names(dq0) <- paste("dq0_",1:r,sep="")

time <- seq(0,totalT,delta_t)

# refines q0 from GBA, finds OGS for the the initial x, assuming balanced growth dq0 = 0
DAE_OGS <- suppressMessages(daspk(y = q0, dy = dq0*0, times = time[1:5], res = Res_DAE, parms = c(1), 
                               rtol = rtol,estini=2))

q0_OGS <- DAE_OGS[2,-1]

# calculates now solution for first 5 time points, in order to find first a consistent initial condition 
# (then we take q at the second time point)
DAE0 <- suppressMessages(daspk(y = q0_OGS, dy = dq0, times = time[1:5], res = Res_DAE, parms = c(1), 
             rtol = rtol,estini=2))

q0_consistent <- DAE0[2,-1]

# estimate derivatives dq. Takes second point 
dqt0 <- (DAE0[-1,-1] -  DAE0[-5,-1])/delta_t

dq0_consistent  <- dqt0[2,]

dq0_consistent[abs(dq0_consistent) < 1e-7] <- 0 # removes numerical artifacts

# "cleans" numerical noise in q0_1 if there is only 1 transporter
if (n_tr == 1) q0[1] <- 1/s[1]

# measuring the total optimization time
st <- system.time({
  
  # solves system of DAE
  
  DAE <- daspk(y = q0_consistent, dy = dq0, times = time, res = Res_DAE, parms = c(1), 
               rtol = rtol,estini=2)  # estini = 2 : estimates q0 from dq0 
  
})

nt <- dim(DAE)[1] 

qt <- DAE[,-1]

# pick second time point for q0, to avoid numerical issues with first point
q0 <- qt[2,]

mu0 <- mu(q0,dq0,0)

tau0 <- tau(0,rho*b(q0))

# "cleans" numerical noise in dqt_1 if there is only 1 transporter
if (n_tr == 1) qt[,1] <- 1/s[1]

nt <- dim(qt)[1] 

time <- DAE[,1]

# estimate derivatives dq. For first point in time, dq0 = 0
dqt <- rbind(dq0,(qt[2:nt,] -  qt[1:(nt-1),])/delta_t)

# smooths second point
dqt[2,] <- (dqt[1,] + dqt[3,])/2

# "cleans" numerical noise in dqt_0 if there is only 1 transporter
if (n_tr == 1) dqt[,1] <- 0

ddqt <- dqt*0
# estimate second derivatives ddq. For first point in time, ddq0 = 0
ddqt <- rbind(dq0*0,(dqt[2:nt,] -  dqt[1:(nt-1),])/delta_t)

ddqt[1,] <- ddqt[2,] 

mu_opt <- 0
for (nc in 1:nt) {
  
  t <- delta_t*(nc - 1)
  
  a <- at(t)
  
  da <- dadt(t)
  
  mu_opt[nc] <- mu(qt[nc,],dqt[nc,],t)
  
}

# SGM
kappat <- t(V)%*%t(qt)

a1 <- 0
for (t in 1:nt) a1[t] <-  at(time[t])[1]

phit   <- matrix(rep(0,r*nt),ncol=r)
for (t in 1:nt) phit[t,] <- phi(qt[t,],dqt[t,],time[t])

# produces cmt
cit <- matrix(rep(0,p*nt),ncol=p)
for (i in 1:nt) cit[i,] <- rho*b(qt[i,])

cmt <- cit[,-p]

# biomass of metabolites
bmt <- cmt/rho

metabolite <- reactant[n_a+1:(n_a+p-2)]

# proteins
pt   <- matrix(rep(0,r*nt),ncol=r)
for (t in 1:nt) pt[t,] <- prot(qt[t,],dqt[t,],time[t])

# estimate dpdt
dpt <- (pt - rbind(pt[1,],pt[1:(nt-1),]))/delta_t

# smooths second point
dpt[2,] <- (dpt[1,] + dpt[3,])/2

# fluxes
vt   <- matrix(rep(0,r*nt),ncol=r)
for (t in 1:nt) vt[t,] <- v(qt[t,],dqt[t,],time[t])

# chit
chit   <- matrix(rep(0,r*nt),ncol=r)
for (t in 1:nt) chit[t,] <- (dpt[t,] + mu_opt[t]*pt[t,])/(M[p,r]*vt[t,r])

# repeats first 3 chit due to numerical artifacts
for (t in 1:3) chit[t,] <- chit[4,]

# "Energy"
varphit <- 0
for (t in 1:nt) varphit[t] <- varphi(qt[t,],t)

varphi0 <- varphi(q0,0)

# grad varphi t
E_it <- matrix(rep(0,(r-1)*nt),ncol=r-1)
for (t in 1:nt) E_it[t,] <- E_i(qt[t,],t)

# velocity t
normdqt <- 0
for (t in 1:nt) normdqt[t] <- sqrt(sum(dqt[t,-1]^2))

# momentum t
normA_it <- 0
for (t in 1:nt) normA_it[t] <- sqrt(sum(A_i(qt[t,],t)^2))

normE_it <- 0
for (t in 1:nt) normE_it[t] <- sqrt(sum(E_i(qt[t,],t)^2))

# distance to OGS
Rt <- 0
for (t in 1:nt) Rt[t] <- R(qt[t,],q0_OGS)

muOGS <- mu(q0_OGS, 0*dq0,0)

# momentum
A_it <- qt[,-1]*0
if (r > 2 ) for (t in 1:nt) A_it[t,] <- A_i(qt[t,],t)
if (r == 2) for (t in 1:nt) A_it[t]  <- A_i(qt[t,],t)

# estimate derivatives dA_i. For first point in time, dA_i0 = 0
if (r > 2)  dA_it <- rbind(A_it[1,]*0,(A_it[2:nt,] -  A_it[1:(nt-1),])/delta_t)
if (r == 2) dA_it <- rbind(A_it[1]*0,(A_it[2:nt] -  A_it[1:(nt-1)])/delta_t)

# smooths second point
if (r > 2 ) dA_it[2,] <- (dA_it[1,] + dA_it[3,])/2
if (r == 2) dA_it[2] <- (dA_it[1] + dA_it[3])/2

# estimate derivatives dmu. For first point in time, dmu0 = 0
dmut <- rbind(mu_opt[1]*0,(mu_opt[2:nt] -  mu_opt[1:(nt-1)])/delta_t)

# smooths second point
dmut[2] <- (dmut[1] + dmut[3])/2

# estimate second derivatives ddmu. For first point in time, ddmu0 = 0
ddmut <- rbind(dmut[1]*0,(dmut[2:nt] -  dmut[1:(nt-1)])/delta_t)

# smooths second point
ddmut[2] <- (ddmut[1] + ddmut[3])/2

# plots ########################################################################
par(mfcol = c(6, 1), mar = numeric(4), oma = c(4, 4, .5, .5), 
    mgp = c(2, .6, 0))

# averaging function
average <- function(input) {
  
  av <- 0
  
  for (t in 1:nt) av[t] <- sum(input[1:t])/t
  
  return(av)
}

# average mu
amu <- mu_opt
for (t in 1:nt) amu[t] <- mean(mu_opt[1:t])

cols <- colorRampPalette(c( "green","darkgreen"))

cole <- colorRampPalette(c( "deepskyblue","navyblue"))

# number of enzymatic reactions
ne <- r - 1 - n_tr

colorj <- c(cols(n_tr),cole(ne),2)

colm <- colorRampPalette(c( "lightgrey","grey50"))

colorm <- colm(p-1)

if (p > 4) colorm <- rainbow(p-1)

coloraxis <- "darkgrey"
colorlab <- "black"

# function to calculate y-axis limits in the plot automatically
ylimits <- function(something) {
  
  yrange <- max(something) - min(something)
  
  if (yrange < 0.0001) yrange <- 0.5
  
  ymax <- max(something) + 0.1*yrange
  
  ymin <- max( min(something) - 0.1*yrange, 0)
  
  return(c(ymin,ymax))
}

# function to calculate y-axis limits in the plot automaticaly
ylimitslab <- function(something) {
  
  yrange <- max(something) - min(something)
  
  if (yrange < 0.0001) yrange <- 0.5
  
  ymax <- max(something) + 0.5*yrange
  
  ymin <- max( min(something) - 0.1*yrange, 0)
  
  return(c(ymin,ymax))
}

# is environment constant ######################################################
const_env <<- 0

# picks many time points
testx <- 0
for (i in 1:nt) testx[i] <- sum(dadt(time[i]))

if (sum(testx == 0) == nt) const_env <<- 1

# plots
plot(time,a1,xlim=c(0,totalT), ylim=ylimits(a1),
     ylab=bquote("Media composition " ~ a^1 ~~ (a.u.)), yaxs="i", xaxs="i",
     type = "l", axes = FALSE)
axis(2L)
box()
# growth rate
plot(time,mu_opt, xlim=c(0,totalT), ylim=ylimits(mu_opt),
     yaxs="i", xaxs="i", ylab=bquote("Specific growth rate " ~ mu ~~ (h^-1)),type = "l", axes = FALSE)
lines(time,amu,lty=2)
if (const_env==1) abline(a=muOGS,b=0,col="grey")
axis(2L)
box()
matplot(time,phit[,r],xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorj[r],ylim=ylimits(phit[,r]), 
        ylab=bquote("Proteome fractions " ~  phi ), yaxs="i", xaxs="i", axes = FALSE)
axis(2L)
box()
matplot(time,phit[,-r],xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorj[-r],ylim=ylimitslab(phit[,-r]), 
        ylab=bquote("Proteome fractions " ~  phi ), yaxs="i", xaxs="i", axes = FALSE)
legend("topright",reaction[-r],col=colorj[-r],pch=16,cex=1,horiz=T)
axis(2L)
box()
matplot(time,bmt,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, 
        frame.plot=F, col= colorm,ylim=ylimitslab(bmt), 
        xlab=bquote("Time" ~ (h)), yaxs="i", xaxs="i", axes = FALSE)
axis(2L)
box()
legend("topright",metabolite,col=colorm,pch=16,cex=1,horiz=T)
matplot(time,chit,xlim=c(0,totalT),pch=16, mgp=c(2, 0.5, 0), col.lab = colorlab, frame.plot=F, 
        col= colorj,ylim= ylimitslab(chit), xlab=bquote("Time" ~ (h)), 
        ylab=bquote("Ribosome allocation " ~  chi ), yaxs="i", xaxs="i")
box()
legend("topright",reaction,col=colorj,pch=16, horiz = T)
mtext(bquote("Time" ~ (h)), side = 1, outer = TRUE, line = 2.2)
mtext(bquote(a[C] ~~ (g/L)), side = 2, outer = TRUE, line = 2.2, adj = 9.5/10)
mtext(bquote(mu ~~ (h^-1)), side = 2, outer = TRUE, line = 2.2, adj = 7.8/10)
mtext(bquote(phi[r])       , side = 2, outer = TRUE, line = 2.2, adj = 5.9/10)
mtext(bquote(phi)         , side = 2, outer = TRUE, line = 2.2, adj = 4.25/10)
mtext(bquote(b           ), side = 2, outer = TRUE, line = 2.2, adj = 2.4/10)
mtext(bquote(chi         ), side = 2, outer = TRUE, line = 2.2, adj = 0.85/10)

################################################################################

# mean optimization time
mean_time <- signif(mean(otime), digits= 3)
