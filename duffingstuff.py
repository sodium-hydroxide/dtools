# Set of functions for duffing system written by Noa3h J. Blair
# Imports packages used
import numpy as np
import sympy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import dtools as dt

# Produces poincare section for the duffing equation when the initial conditions are (0.5,0.5)
# the parameter g gives the strength of the driving force 
def duffing_run(g):
    q, p = dt.poincare(dt.duffing,0.1,0.1,g,50)
    return q, p

# Gets the values for a given value of g
q, p  = duffing_run(20.3)

# generates a plot of the the separatrix for the unforced problem. Points inside the separatrix are negative energy,
# those outside are positive energy
Q, P = np.meshgrid(np.arange(-1.5,1.5,0.025),np.arange(-1.5,1.5,0.025))
H = 0.5*(P**2) - 0.5*(Q**2) + 0.25*(Q**4)
plt.contour(Q,P,H,0)

# Plots a scatter plot of the points
plt.scatter(q,p, s = 10)
# shows the plot
plt.show()