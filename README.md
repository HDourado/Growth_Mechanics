# Growth_Mechanics
R codes for numerical simulations in Growth Mechanics and Growth Balance Analysis

1) Run "Interface.R" 

2) Use the function (modelname,predict.parameters,is.reversible,delta_t,totalT), where the inputs are

# modelname: name of a document in the Models file in quotes, e.g. "L3"

# predict.parameters: 0 = does not predict, take from file, 1,2,3,... = predicts using this parameter as K_ratio, as described in "methods".

# is.reversible (only relevant if what to predict kinetic parameters): 0 = not , 1 = yes.

# delta = interval of time in hours

# totalT: total time in hours

# Example:
# GM("L3",3,0,0.001,3)
