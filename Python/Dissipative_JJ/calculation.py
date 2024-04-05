import numpy as np
import scipy as sp
import cal_parameter as param
import Hamiltonian as H
import K_sum as K

size = [21, 30, 42, 51] 

g_arr = [0 for i in range(21)]

for s in size:
    for i in range(21):
        if i==0:
            g_arr[i] = i
        if i>0:
            g_arr[i] = round(g_arr[i-1] + 0.01,2)



    garray = [0 for i in range(25)]
    gsquare = [0 for i in range(len(garray))]
    
    for i in range(len(garray)):
        if (i<20):
            garray[i+1] = round(garray[i]+0.05,2)
            gsquare[i+1] = round(garray[i+1]**2,4)
        elif (i<24):
            gsquare[i+1] = round(gsquare[i] + 0.2,1)


    for i in g_arr:
        HM = H.Hamiltonian_Matrix(param.gamma,s,i,1)

        domain = param.tau_grid
        value = H.Chi_sp(HM,param.tau_grid)

        df = np.column_stack((domain.real,value.real,value.imag))

        np.savetxt('Exacttest_size{}_beta{}_omega{}_g_{}.txt'.format(s,param.beta,param.omega,i),df)


    for i in gsquare:
        HM = H.Hamiltonian_Matrix(param.gamma,s,i,1)

        domain = param.tau_grid
        value = H.Chi_sp(HM,param.tau_grid)

        df = np.column_stack((domain.real,value.real,value.imag))

        np.savetxt('Exacttest_size{}_beta{}_omega{}_g_{}.txt'.format(s,param.beta,param.omega,i),df)