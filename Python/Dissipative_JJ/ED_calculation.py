import numpy as np
import scipy as sp
import Hamiltonian as H

################# setting beta value #####################

betamaxarray = [0 for i in range(5)]

for j in range(1,len(betamaxarray),1):
    betamaxarray[j] = round(0.5 + betamaxarray[j-1],2)
    betamaxarray[j-1] = betamaxarray[j]

'''
for j in range(11):
    if j == 0:
        betamaxarray[j] = 0.001
    if j > 0:
        betamaxarray[j] = round(j * 0.1,1)
'''

################ calculate coupling value for each beta #########################

for b in betamaxarray:

    ############### for b value ########################

    ta = np.linspace(0,b,201)
    tau = np.linspace(0,b,201)
    k = np.full(100,1)

    ############## coupling in one bath ####################
    garray = [0 for i in range(25)]
    gsqaure = [0 for i in range(len(garray))]
    for i in range(len(garray)):
        if (i<20):
            garray[i+1] = round(garray[i]+0.05,2)
            gsqaure[i+1] = round(garray[i+1]**2,4)
        elif (i<24):
            gsqaure[i+1] = gsqaure[i] + 1
    
    newgarray = np.linspace(1,2,11)

    ############## plotting function #######################

    def plotting(gamma,g,beta):
        A =[]
        for i in range(len(ta)):
            plot = H.Chi_sp(beta,gamma,k,g,7,ta[i])
            A.append(plot)

        return(A)

    for i in gsqaure:
        domain = tau
        value = np.array(plotting(1,i,b))

        df = np.column_stack((domain.real,value.real,value.imag))

        np.savetxt('Exacttest_grid200_size21_beta{}_g_{}.txt'.format(b,i),df)