import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import Hamiltonian as H
import K_sum as K

ta = np.linspace(0,1,400)
tau = np.linspace(0,1,400)
k = np.full(100,1)
b = 1

garray = [0 for i in range(25)]
gsqaure = [0 for i in range(len(garray))]
for i in range(len(garray)):
    if (i<20):
        garray[i+1] = round(garray[i]+0.05,2)
        gsqaure[i+1] = round(garray[i+1]**2,4)
    elif (i<24):
        gsqaure[i+1] = gsqaure[i] + 1


def plotting(r,g):
    A =[]
    for i in range(len(ta)):
        plot = H.Chi_sp(1,r,k,g,1,ta[i])
        A.append(plot)

    return(A)

for i in range(len(garray)):
    a = tau
    b = plotting(1,gsqaure[i])

    df = np.column_stack((a,b))

    np.savetxt('Exacttest_beta0_001_g_{}.txt'.format(gsqaure[i]),df)