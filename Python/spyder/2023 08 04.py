import numpy as np
import matplotlib.pyplot as plt

def standard_even(r,z):

    A = [[0 for j in range(z)] for k in range(z)]

    for i in range(z):
        for j in range(z):
            try:
                if i==j:
                    A[i][j] = (i)**2
                if abs(i-j) == 1:
                    A[i][j] = -r/2
                else:
                    None
            except:
                None

    test_A = np.array(A)
    test_A[0][1] = -r/(2*np.sqrt(2))
    test_A[1][0] = -r/(2*np.sqrt(2))
    A_eig=np.linalg.eig(test_A)

    #print(A_eig[0])

    return A_eig[0]

def standard_odd(r,z):

    A = [[0 for j in range(z)] for k in range(z)]

    for i in range(z):
        for j in range(z):
            try:
                if i==j:
                    A[i][j] = (i+1)**2
                if abs(i-j) == 1:
                    A[i][j] = -r/2
                else:
                    None
            except:
                None

    test_A = np.array(A)
    A_eig=np.linalg.eig(test_A)

    #print(A_eig[0])

    return A_eig[0]

def eigenplotter(r,z):
    ar = np.linspace(-r,r,r*100)
    arr = []
    arr_e = []

    for i in ar:
        arr.append(standard_odd(i,z))
        arr_e.append(standard_even(i,z))

    plt.figure(figsize=(10,8))
    plt.plot(r,arr)
    plt.plot(r,arr_e)
    #plt.grid(True,linestyle='--')
    plt.show()
    
eigenplotter(5,20)
