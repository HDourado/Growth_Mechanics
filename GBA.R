# GBA ##########################################################################

# Reads model saved as .ods file ###############################################

suppressMessages(source("Readmodelods.R"))

if (is.reversible == 1) modelname <- paste(modelname,"_rev",sep="")

# Medium concentrations ########################################################

# Sigmoid function for decreasing concentration of first medium concentration a_1, other a = 10
if (n_a == 1) at <-function(t) c(-9.999*t^2/(0.1 + t^2) + 10)
if (n_a >  1) at <-function(t) c(-9.999*t^2/(0.1 + t^2) + 10 , rep(10,n_tr-1))

n_conditions <<- 20

# kinetics #####################################################################

source("Kinetics.R")

# Singular Growth Modes ########################################################

source('SGM.R')

# Predicts kinetic parameters based on mu and phi data #########################

if (predict.parameters > 0) source("Parameter_prediction.R")

# Optimization on f ############################################################

source("GBA_solver.R") 

# Exporting results ############################################################

# Exporting csv file with results #######

source("GBA_Exportcsv.R")

# Plots #################################

source("GBA_Plots.R")

# dev.off()
