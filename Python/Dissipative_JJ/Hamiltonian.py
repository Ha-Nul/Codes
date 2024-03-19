import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import cmath
import random

np.set_printoptions(threshold=784,linewidth=np.inf)

## Full source code available : 2023 09 06.ipynb

#####################Local Hamiltonian##########################

def Local_odd(gamma, n):
    '''Independent Function, Creates H_0 Hamiltonian in sine basis,
        n is dimensinon of matrix'''

    A = [[0 for j in range(n)] for k in range(n)]

    for i in range(n):
        for j in range(n):
            try:
                if i==j:
                    A[i][j] = ((i+1)**2)
                if abs(i-j) == 1:
                    A[i][j] = -gamma/2
                else:
                    None
            except:
                None

    return np.array(A)

def Local_odd_Eigenvec(gamma: float, n: int ,m: int):
    ''' Requires Local_odd, m is mth eigenvector of H_0 Matrix, corresponds with mth eigenvalue'''

    A = Local_odd(gamma,n)
    A_eigvec = np.linalg.eig(A)
    A_trans = np.transpose(A_eigvec[1])

    return A_trans[m]

def Local_odd_Eigenval(gamma: float,n: int,m: int):
    ''' Requires Local_odd, m it mth eigenvalue of H_0 Matrix, corresponds with mth eigenvector'''

    A = Local_odd(gamma,n)
    A_eigval = np.linalg.eig(A)

    return A_eigval[0][m]

def Local_even(gamma:float, n:int):
    '''Independent Function, Creates H_0 Hamiltonian in cosine basis
    n = dimension of matrix'''

    A = [[0 for j in range(n)] for k in range(n)]

    for i in range(n):
        for j in range(n):
            try:
                if i==j:
                    A[i][j] = (i)**2
                if abs(i-j) == 1:
                    A[i][j] = -gamma/2
                else:
                    None
            except:
                None

    np_A = np.array(A)
    np_A[0][1] = -gamma/(np.sqrt(2))
    np_A[1][0] = -gamma/(np.sqrt(2))

    return np_A

def Local_even_Eigenvec(gamma: float,n: int ,m: int):
    ''' Requires Local_even, m is mth eigenvector of H_0 Matrix, corresponds with mth eigenvalue'''

    A = Local_even(gamma,n)
    A_eig = np.linalg.eig(A)
    A_trans = np.transpose(A_eig[1])

    return A_trans[m]

def Local_even_Eigenval(gamma: float,n: int,m: int):
    ''' Requires Local_even, m it mth eigenvalue of H_0 Matrix, corresponds with mth eigenvector'''

    A = Local_even(gamma,n)
    A_eig = np.linalg.eig(A)

    return A_eig[0][m]

def Elements(d,x,y):

    x_trans = np.transpose(x)

    arr = np.matmul(d,y)
    arr2 = np.matmul(x_trans,arr)

    return arr2

#################################################################
#################### Prerequisites for Interacting Hamiltonian #######################

def INT_odd_Eigenvec(gamma, n, m):
    LOC_VEC_O = Local_odd_Eigenvec(gamma, n, m)
    INT_ARR_O = []

    for i in range(n):
        INT_ARR_O.append(-1j * (i+1) * LOC_VEC_O[i])
    
    return np.array(INT_ARR_O)

def INT_even_Eigenvec(gamma, n, m):
    LOC_VEC_E = Local_even_Eigenvec(gamma, n, m)
    INT_ARR_E = []

    for i in range(n):
        INT_ARR_E.append(1j * i * LOC_VEC_E[i])

    return np.array(INT_ARR_E)

#################### Main Hamiltonian ##########################

def Hamiltonian_Matrix(gamma: float,n: int, g: float,omega: float):
    '''Return Hamiltonian H in matrix form, requires prerequisites,
    r : Value of Gamma, z : dimension of Matrix, g : coupling strength, omega : frequency '''

    ######### Elements of Hamiltonian matrix ###############
    
    # Eigenvalues for local and bath part
    LOC_EV_ODD = Local_odd_Eigenval(gamma,n,0)
    LOC_EV_EVE_g = Local_even_Eigenval(gamma,n,0)
    LOC_EV_EVE_s = Local_even_Eigenval(gamma,n,1)

    # Eigenvectors for N matrix
    LOC_ODD = Local_odd_Eigenvec(gamma,n,0)
    LOC_EVE_g = Local_even_Eigenvec(gamma,n,0)
    LOC_EVE_s = Local_even_Eigenvec(gamma,n,0)

    N_ODD = INT_odd_Eigenvec(gamma, n, 0)
    N_EVE_g = INT_even_Eigenvec(gamma, n, 0)
    N_EVE_s = INT_even_Eigenvec(gamma, n, 1)

    N_ELE_og = np.dot(LOC_ODD,N_EVE_g)
    N_ELE_go = -N_ELE_og #first state dot product to ground state
    N_ELE_os = np.dot(LOC_ODD,N_EVE_s)
    N_ELE_so = -N_ELE_os
    

    # N_Matrix (interaction Hamiltonian)
    A = [[0 for i in range(n)] for j in range(n)]

    for i in range (n):
        for j in range (n):
            try:
                if i==j and i%3 == 0:
                   A[i][j+4] = (np.sqrt((j+4)//3)) * g * N_ELE_go
                if i==j and i%3 == 1:
                    A[i][j+2] = (np.sqrt((j+2)//3)) * g * N_ELE_og
                    A[i][j+4] = (np.sqrt((j+4)//3)) * g * N_ELE_os
                if i==j and i%3 == 2:
                    A[i][j+2] = (np.sqrt((j+2)//3)) * g * N_ELE_so
            except:
                None

    np_A = np.array(A)
    np_B = np.transpose(np_A.copy())
    
    # Hamiltonian_local , Hamiltonian_bath
    for i in range(n):
        for j in range(n):
            try:
                if i==j and i%3==0:
                    A[i][j] = LOC_EV_EVE_g + (j//3)
                elif i==j and i%3==1:
                    A[i][j] = LOC_EV_ODD + (j//3)
                elif i==j and i%3==2:
                    A[i][j] = LOC_EV_EVE_s + (j//3)
                else:
                    A[i][j] = 0
            except:
                None
    
    np_C = np.array(A)
    Array = np_A + np_B + np_C

    return Array

def Hamiltonian_Matrix_Eigenval(gamma,n,g,omega,i):
    A_eigval = np.linalg.eig(Hamiltonian_Matrix(gamma,n,g,omega))[0][i]
    return A_eigval

def Hamiltonian_Matrix_Eigenvec(gamma,n,g,omega,i):
    A = np.linalg.eig(Hamiltonian_Matrix(gamma,n,g,omega))[1].T
    A_eigvec = A[i]
    return A_eigvec

## Lie group ####################################
def Lie_group(x):
    sigma = [[0 for i in range(3)] for j in range(3)]

    if x == 1:
        sigma[0][1] = 1
        sigma[1][0] = 1
        return np.array(sigma)
    if x == 2:
        sigma[0][1] = 1j
        sigma[1][0] = -1j
        return np.array(sigma)
    if x == 3:
        sigma[0][0] = 1
        sigma[1][1] = 1
        return np.array(sigma)
    if x == 4:
        sigma[0][2] = 1
        sigma[2][0] = 1
        return np.array(sigma)
    if x == 5:
        sigma[0][2] = -1j
        sigma[2][0] = 1j
        return np.array(sigma)
    if x == 6:
        sigma[1][2] = 1
        sigma[2][1] = 1
        return np.array(sigma)
    if x == 7:
        sigma[1][2] = -1j
        sigma[2][1] = 1j
        return np.array(sigma)
    if x == 8:
        sigma[0][0] = 1/(3**0.5)
        sigma[2][2] = -2/(3**0.5)
        return np.array(sigma)

def Lie_tensorproduct(x,y):
    A = np.identity(y)
    Lie_Tens = np.kron(A,Lie_group(x))
    return Lie_Tens
## Lie group ####################################
## Spectral Function ############################
def Spectral_Function(beta: float, gamma :float ,n: int,g: float,matsufreq: float,eta:float,omega: float):
    '''Spectral Density of Correlation_Function. r = Value of gamma, z = dimension/2 of Matrix,
    g = Coupling strength, matsufreq = Matsubara frequency of given Function, eta = Value for analytic continuation'''
    
    # Tr Z
    eig_val = np.array([[Local_even_Eigenval(gamma,3,0),0,0],
         [0,Local_odd_Eigenval(gamma,3,0),0],
         [0,0,Local_even_Eigenval(gamma,3,1)]])
    exp_val = sp.linalg.expm(-beta*eig_val)
    Z = np.sum(np.array(exp_val))

    # Spectral
    A = []
    for i in range(3*n):
        for j in range(3*n):
            expec_nm = Elements(Lie_tensorproduct(1,n),Hamiltonian_Matrix_Eigenvec(gamma,3*n,g,omega,i),Hamiltonian_Matrix_Eigenvec(gamma,3*n,g,omega,j))
            conju = np.conjugate(expec_nm)

            n = Hamiltonian_Matrix_Eigenval(gamma,3*n,g,omega,i)
            m = Hamiltonian_Matrix_Eigenval(gamma,3*n,g,omega,j)

            denom = (matsufreq + n - m)**2 + eta**2
            #sign of numerator e^E_m may can change
            numer = np.exp(-beta*n)-np.exp(-beta*m)
            value = (conju*expec_nm*numer*2*eta)/denom

            #print(denom)

            A.append(value)
                    
    #np.imag() was not used.
    np_A = np.array(A)
    sum_A = np.sum(np_A)

    return sum_A/Z

## Spectral function ##########################

## Chi_function ###############################

def Chi_sp(beta : float, gamma : float, omega: float, g: float, n: int , tau: float):
    '''Chi_function, r : value of gamma, z : dimension/3 of hilbert space, g : coupling strength, omega : frequency'''

    # Tr Z
    Exp_val = sp.linalg.expm(-beta*Hamiltonian_Matrix(gamma,3*n,g,omega))
    Z = np.trace(Exp_val)

    # Chi_calculation
    A = []
    for i in range(3*n):
        for j in range(3*n):
            E_N = Hamiltonian_Matrix_Eigenval(gamma,3*n,g,omega,i)
            E_M = Hamiltonian_Matrix_Eigenval(gamma,3*n,g,omega,j)
            Exp_val_chi = np.exp(-(beta-tau) * E_N - tau * E_M)

            Expec_NM = Elements(Lie_tensorproduct(1,n),Hamiltonian_Matrix_Eigenvec(gamma,3*n,g,omega,i),Hamiltonian_Matrix_Eigenvec(gamma,3*n,g,omega,j))
            #Expec_NM = Hamiltonian_Matrix_Eigenvec(gamma,3*z,g,omega,j).T @ Lie_group(1) @ Hamiltonian_Matrix_Eigenvec(gamma,3*z,g,omega,i)

            A.append(Exp_val_chi * (Expec_NM**2))
    
    #print(A)
    # Chi_total
    np_A = np.array(A)
    sum_A = np.sum(np_A)

    return sum_A/Z
