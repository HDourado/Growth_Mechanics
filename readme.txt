# Functions for GBA and GM

# 1) Run the "Interface.R" script

# 2) For GBA, run GBA(modelname,predict.parameters,is.reversible) , where

# modelname: name of a document in the Models folder in quotes, e.g. "L3"

# predict.parameters: 0 = does not predict, take from file, 1,2,3,... = predicts using this parameter as K_ratio,
# as described in "methods" of https://www.biorxiv.org/content/10.1101/2025.06.24.661369v1.

# is.reversible (only relevant if what to predict kinetic parameters): 0 = not , 1 = yes.

# Then check the folder "Results GBA" for result plots and data

# Example: 
# GBA("L3",3,0)

# 3) For GM, run GM(modelname,predict.parameters,is.reversible,delta_t,totalT) , where

# modelname: name of a document in the Models folder in quotes, e.g. "L3"

# predict.parameters: 0 = does not predict, take from file, 1,2,3,... = predicts using this parameter as K_ratio,
# as described in "methods" of https://www.biorxiv.org/content/10.1101/2025.06.24.661369v1.

# is.reversible (only relevant if what to predict kinetic parameters): 0 = not , 1 = yes.

# delta = interval of time in hours

# totalT: total time in hours

# Then check the folder "Results GM" for result plots and data

# Example:
# GM("L3",3,0,0.001,3)

# for accurate Lambda calculations, use smaller dt, e.g. dt = 0.00001
# GM("L3",3,0,0.00001,3)


