import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import Hamiltonian as H
import K_sum as K

ta = np.linspace(0,1,400)
tau = np.linspace(0,1,400)
k = np.full(100,1)
b = 1

def plotting(r,g):
    A =[]
    for i in range(len(ta)):
        plot = H.Chi_sp(1,r,k,g,8,ta[i])
        A.append(plot)

    return(A)

df = np.column_stack((tau,plotting(1,0.1)))
np.savetxt('NCA_py.txt',df)