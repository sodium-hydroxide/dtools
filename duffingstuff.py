# Set of functions for duffing system written by Noa3h J. Blair
# Imports packages used
import numpy as np
import sympy as sp
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

import dynamicstools as dt

#dt.duffing_poincare_plot(0.5,0.3,1.0,500)


y = np.zeros(20)
x = np.zeros(20)

for j in range(20):
    x[j] = j + 1
   # y[j] = 4*j + 1 + 4
    A = 1
    for k in range(j):
        A = A*(4*k + 1)
    y[j] = A

print(x)
print(y)