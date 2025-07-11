# Loads model saved as .ods file ###############################################

setwd(paste(directory,"/Models",sep=""))

odsfile <- paste(modelname,".ods", sep = "")

sheets <- list_ods_sheets(odsfile)

nsheets    <- length(sheets)

# Position of parameters
posM          <- (1:nsheets)[sheets == "M"]
posK          <- (1:nsheets)[sheets == "K"]
posKI         <- (1:nsheets)[sheets == "KI"]
posKA         <- (1:nsheets)[sheets == "KA"]
poskcat       <- (1:nsheets)[sheets == "kcat"]
posrho        <- (1:nsheets)[sheets == "rho"]
posat         <- (1:nsheets)[sheets == "at"]
posq          <- (1:nsheets)[sheets == "q"]
posphi        <- (1:nsheets)[sheets == "phi"]
posrmw        <- (1:nsheets)[sheets == "rmw"]
posminphi     <- (1:nsheets)[sheets == "min_phi"]
posmaxphi     <- (1:nsheets)[sheets == "max_phi"]
posminf       <- (1:nsheets)[sheets == "min_q"]
posmaxf       <- (1:nsheets)[sheets == "max_q"]
posminc       <- (1:nsheets)[sheets == "min_c"]
posmaxc       <- (1:nsheets)[sheets == "max_c"]
posribcomp    <- (1:nsheets)[sheets == "ribcomp"]

# Reads q0, qT, T ##############################################################

#q0  <- as.numeric(read_ods(odsfile, sheet= posq)[1,-1])
dq0 <- as.numeric(read_ods(odsfile, sheet= posq)[1,-1]) 

# Reads a(t) function ##########################################################

eval(parse(text=noquote(read_ods(odsfile, sheet= posat)[1,2]))) 

# Reads dadt(t) function #######################################################

eval(parse(text=noquote(read_ods(odsfile, sheet= posat)[2,2])))

# Getting data from sheets

# Mass fraction matrix Mtotal including external reactants #####################

Mtotal <- as.matrix(read_ods(odsfile, sheet= posM)[,-1])

# reaction and reactant names ##################################################

reaction <- colnames(Mtotal)
reactant <- unlist(read_ods(odsfile, sheet= posM)[,1])

rownames(Mtotal) <- reactant

# Density rho ##################################################################

rho <- as.numeric(read_ods(odsfile, sheet= posrho)[1,2])

# Michaelis constant matrix K ##################################################
if (length(posK) > 0) {
  
  K <- as.matrix(read_ods(odsfile, sheet= posK)[,-1])
  
} else K <- 0.1*(Mtotal<0)

rownames(K) <- reactant

# inhibition constant matrix KI ################################################
if (length(posKI) > 0) {
  
  KI <- as.matrix(read_ods(odsfile, sheet= posKI)[,-1])
  
} else KI <- 0*K

rownames(KI) <- reactant

# activation constant matrix KA ################################################
if (length(posKA) > 0) {
  
  KA <- as.matrix(read_ods(odsfile, sheet= posKA)[,-1])
  
} else KA <- 0*K

rownames(KA) <- reactant

# kcat
kcatf <- as.numeric(read_ods(odsfile, sheet= poskcat)[1,-1])
kcatb <- as.numeric(read_ods(odsfile, sheet= poskcat)[2,-1])

# Experimental data ############################################################

exp_data <- as.matrix(read_ods(odsfile, sheet= posphi))

# reactant molecular weights ###################################################

if (length(posrmw) > 0) rmw <- as.matrix(read_ods(odsfile, sheet= posrmw)[,-1])

# reads ribosome composition constraint ########################################

if (length(posribcomp) > 0) {
  
  ribcomp <- as.numeric(read_ods(odsfile, sheet= posribcomp)[1,1])
  
} else  ribcomp <- 0

# Definitions ##################################################################

# number of external reactants
n_a <- length(at(0))

# internal matrix M
M <- Mtotal[-c(1:n_a),]

# names of internal reactants
i_reactant <- reactant[-c(1:n_a)]

# number of internal reactants
p <- dim(M)[1]

# number of reactions
r <- dim(M)[2]

# the sum of each M column 
s <- colSums(M)

# delete numerical artifacts
s[abs(s) < 1e-10] <- 0

# number of transport reactions
n_tr <- length(c(1:(r-1))[s[1:(r-1)] != 0] )

# bounds for minimal and maximal c #############################################

if (length(posminc) > 0) {
  
  min_c <- as.numeric(as.matrix(read_ods(odsfile, sheet= posminc)[,2]))
  
} else min_c <- rep(0,p)

if (length(posmaxc) > 0) {
  
  max_c <- as.numeric(as.matrix(read_ods(odsfile, sheet= posmaxc)[,2]))
  
} else max_c <- rep(rho,p)

# bounds for minimal and maximal f #############################################
if (length(posminf) > 0) {
  
  min_q <- as.numeric(read_ods(odsfile, sheet= posminf)[1,])
  
} else min_q <- rep(0,r)

if (length(posmaxf) > 0) {
  
  max_q <- as.numeric(read_ods(odsfile, sheet= posmaxf)[1,])
  
} else max_q <- rep(100,r)

# bounds for minimal phi #######################################################
if (length(posminphi) > 0) {
  
  min_phi <- as.numeric(read_ods(odsfile, sheet= posminphi)[1,])
  
} else min_phi <- rep(0,r)

# bounds for maximal phi
if (length(posmaxphi) > 0) {
  
  max_phi <- as.numeric(read_ods(odsfile, sheet= posmaxphi)[1,])
  
} else max_phi <- rep(1,r)


setwd(directory)
