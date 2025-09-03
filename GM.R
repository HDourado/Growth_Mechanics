# GM solver ####################################################################

# Singular Growth Modes ########################################################

source("SGM.R")

# Predicts kinetic parameters based on mu and phi data #########################

if (predict.parameters > 0 & r > 2) source("Parameter_prediction.R")

# kinetics #####################################################################

source("Kinetics.R")

# Refines initial q0 using GBA #################################################

suppressMessages(source('q0GBA.R'))

# Solve ########################################################################

if (r-rM == 0) suppressMessages( source("GM_DAEsolver.R") )

if (r-rM > 0) stop("The matrix M is not full column rank.")

# Exporting results ############################################################

# Exporting csv file with results #######

source("GM_Exportcsv.R")

# Plots #################################

source("GM_Frequency.R")

#source("GM_Plots.R")

source("GM_main.R")

source("GM_averages.R")

#if (r == 3) source("GM_fields.R")

# Print results #########################

print("--------------------------------------------")

if (min(cit) < 0) print("negative c!")

if (min(phit) < 0) print("negative phi!") 

if (min(chit[-c(1:5),]) < 0) print("negative chi!") 

################################################################################

print(paste("nu      =", nu) )

print(paste("mu*     =", muOGS))

print(paste("Lambda  =", Lambda))

print(paste("varphi0 =", varphi0))

# ratio between cell cycle frequency/biomass allocation frequency
print(paste("nu_c/nu =", Lambda/log(2,base = exp(1))/nu) )

print("--------------------------------------------")
