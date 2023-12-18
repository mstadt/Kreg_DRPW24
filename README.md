# K regulation W24 WiM DRP
This is the K regulation model used for the W24 WiM DRP.

# Driver files
### driverSS.m
Run this script to compute the steady state result for the baseline
K model.

### driverMeal.m
Run this script to run a single meal simulation. 
Option to save the simulation results to /MealSim/.

### plot_3MealTypes.m
This script can be used to plot 3 simulations from driverMeal.m.

# Functions
### kreg_eqns 
This file contains the model equations.

### set_params
This file sets the parameter values. 
Use **par2vector** to change to a vector.
