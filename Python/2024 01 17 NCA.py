import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import Hamiltonian as H
import K_sum as K

betamaxarray = [0 for i in range(11)]

for j in range(11):
    if j == 0:
        betamaxarray[j] = 0.001
    if j > 0:
        betamaxarray[j] = round(j * 0.1,1)

for b in betamaxarray:

    ta = np.linspace(0,b,400)
    tau = np.linspace(0,b,400)
    k = np.full(100,1)

    garray = [0 for i in range(25)]
    gsqaure = [0 for i in range(len(garray))]
    for i in range(len(garray)):
        if (i<20):
            garray[i+1] = round(garray[i]+0.05,2)
            gsqaure[i+1] = round(garray[i+1]**2,4)
        elif (i<24):
            gsqaure[i+1] = gsqaure[i] + 1


    def plotting(r,g,beta):
        A =[]
        for i in range(len(ta)):
            plot = H.Chi_sp(beta,r,k,g,1,ta[i])
            A.append(plot)

        return(A)

    for i in range(len(garray)):
        a = tau
        b = plotting(1,gsqaure[i],b)

        df = np.column_stack((a,b))

        np.savetxt('Exacttest_beta0_{}_g_{}.txt'.format(b,gsqaure[i]),df)
