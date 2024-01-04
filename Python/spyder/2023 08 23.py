import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import cmath

np.set_printoptions(threshold=784,linewidth=np.inf)

## prerequisites ##########################

def Standard_even(r: float ,z: int):
    '''Independent Function, Creates H_0 Hamiltonian in cosine basis
    r = Value of Gamma, z = Dimension of Matrix'''

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
    test_A[0][1] = -r/(np.sqrt(2))
    test_A[1][0] = -r/(np.sqrt(2))

    return test_A

def Standard_even_Eigenvec(r,z,m):

    Matrix = Standard_even(r,z)
    A_eig = np.linalg.eig(Matrix)
    A_trans = np.transpose(A_eig[1])

    return A_trans[m]

def eigenval_even(r,z,m):

    Matrix = Standard_even(r,z)
    A_eig = np.linalg.eig(Matrix)

def eigenval_even(r,z,m):

    Matrix = Standard_even(r,z)
    A_eig = np.linalg.eig(Matrix)

    return A_eig[0][m]

## Prerequisites #############################

def Elements(d,x,y):

    #x_trans = np.transpose(x)

    arr = np.matmul(d,y)
    arr2 = np.matmul(x,arr)

    return arr2

def Hamiltonian_Matrix(r,z,g,omega):

    A = [[0 for j in range(z)] for k in range(z)]

    a_up = g*Elements(Standard_even(0,z),Standard_even_Eigenvec(r,z,0),Standard_even_Eigenvec(r,z,0))
    a_down = g*Elements(Standard_even(0,z),Standard_even_Eigenvec(r,z,1),Standard_even_Eigenvec(r,z,1))
    a_diagonal = g*Elements(Standard_even(0,z),Standard_even_Eigenvec(r,z,0),Standard_even_Eigenvec(r,z,1))


    for i in range(z):
        for j in range(z):
            try:

                if i==j and j%2 == 0:
                    A[i][j] = i*omega/2 + eigenval_even(r,z,0)
                    A[i+1][j+1] = i*omega/2 + eigenval_even(r,z,1)

                if (j-i) == 2 and j%2 == 0 or j%2 == 2: #파이썬은 행렬 index가 0부터 시작함.
                    A[i][j] = a_up*np.sqrt(j/2)
                    A[i+1][j+1] = a_down*np.sqrt(j/2)
                    A[i][j+1] = a_diagonal*np.sqrt(j/2)
                    A[i+1][j] = a_diagonal*np.sqrt(j/2)

                if z%2 != 0 and i == z-2 and j == z-1 :
                    A[i][j] = a_diagonal*np.sqrt(j/2)

                else:
                    None
            except:
                None

    np_A = np.array(A)
    np_B = np.array(A)
    sum_A = np.transpose(np_B)

    for i in range(z):
        for j in range(z):
            if i == j:
                sum_A[i][j] = 0

    test_A = sum_A + np_A


    return test_A

def Hamiltonian_Matrix_Eigenval(r,z,g,omega):
    
    A = [[0 for i in range(z)]for j in range(z)]

    for i in range(z):
        for j in range(z):
            if i == j:
                A[i][j] = np.linalg.eig(Hamiltonian_Matrix(r,z,g,omega))[0][i]

    return np.array(A)

def Hamiltonian_Matrix_Eigenvec(r,z,g,omega):
    return np.linalg.eig(Hamiltonian_Matrix(r,z,g,omega))[1].T

def pauli(x):
    sigma = [[0 for i in range(2)] for j in range(2)]

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
    
def Pauli_tensorproduct(x,y):
    A = np.identity(y)

    Pauli_Tens = np.kron(A,pauli(x))

    return Pauli_Tens

def Spin_correlation(beta,r,tau,z,g):

    Correlation = []
    Z = []

    Inner_e_1 = sp.linalg.expm(-beta*Hamiltonian_Matrix_Eigenval(r,2*z,g,1))
    Inner_e_2 = sp.linalg.expm(tau*Hamiltonian_Matrix_Eigenval(r,2*z,g,1))
    Inner_e_3 = sp.linalg.expm(-tau*Hamiltonian_Matrix_Eigenval(r,2*z,g,1))

    matmul1 = np.matmul(Inner_e_3,Pauli_tensorproduct(1,z))
    matmul2 = np.matmul(Pauli_tensorproduct(1,z),matmul1)
    matmul3 = np.matmul(Inner_e_2,matmul2)
    matmul4 = np.matmul(Inner_e_1,matmul3)

    for i in range(2*z):
        Correlation.append(Elements(matmul4,Hamiltonian_Matrix_Eigenvec(r,2*z,g,1)[i],Hamiltonian_Matrix_Eigenvec(r,2*z,g,1)[i]))
        Z.append(sp.linalg.expm(-beta*Hamiltonian_Matrix_Eigenval(r,2*z,g,1)[i][i]))

    np_Cor = np.array(Correlation)
    np_Z = np.array(Z)
    
    sum_Cor = np.sum(np_Cor)
    sum_Z = np.sum(np_Z)

    return sum_Cor/sum_Z

def Time(x):
    Time = []
    
    tau = np.linspace(0,1,200)

    for i in range(len(tau)):
        a = Spin_correlation(0,1,tau[i],2,x)
        Time.append(a)

    return Time

def Spectral_Function(beta,r,z,g,omega,eta):
    Z = []
    for i in range(2*z):
        Z.append(sp.linalg.expm(-beta*Hamiltonian_Matrix(r,2*z,g,1)[i][i]))
    
    Z = np.array(Z)
    Z_sum = np.sum(Z)

    A = [[0 for i in range(2*z)] for j in range(2*z)]
    
    for i in range(2*z):
        for j in range(2*z):
            expec = Elements(Pauli_tensorproduct(1,z),Hamiltonian_Matrix_Eigenvec(r,2*z,g,1)[i],Hamiltonian_Matrix_Eigenvec(r,2*z,g,1)[j])          
            conju = np.conjugate(expec)

            n = Hamiltonian_Matrix(r,2*z,g,1)[i][i]
            m = Hamiltonian_Matrix(r,2*z,g,1)[j][j]

            denom =  (omega + n - m)**2 + eta**2
            #sign of numerator e^E_m may can change 
            numer = np.exp(-beta*n)-np.exp(-beta*m)

            value = expec*conju*numer*2*eta/(denom)

            A[i][j] = value
            
    np_A = np.array(A)
    sum_A = np.sum(A)

    return sum_A/Z_sum


