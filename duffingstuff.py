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
    q, p = dt.poincare(dt.duffing,0.5,0.0,g,100)
    return q, p

q, p  = duffing_run(0.0)

plt.scatter(q,p)

plt.show()