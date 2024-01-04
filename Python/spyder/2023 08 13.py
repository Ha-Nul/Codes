import numpy as np
import matplotlib.pyplot as plt
import math

r = np.linspace(0,30,200)


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
    test_A[0][1] = -r/(2)
    test_A[1][0] = -r/(np.sqrt(2))

    return test_A

def eigenvec_even(r,z,m):

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
    test_A[0][1] = -r/(2)
    test_A[1][0] = -r/(np.sqrt(2))
    A_eig=np.linalg.eig(test_A)

    return A_eig[1][m]

def Elements(d,x,y):

    N = standard_even(0,d)

    arr = np.matmul(N,y)
    arr2 = np.matmul(x,arr)

    return arr2

def Mat(r,d):
    Mat = []

    for i in r:
        one = eigenvec_even(i,d,0)
        two = eigenvec_even(i,d,1)

        rix = [[Elements(d,one,one),Elements(d,one,two)],[Elements(d,two,one),Elements(d,two,two)]]
        Mat.append(rix)

        #print(rix)

    return Mats

oone =[]
ttwo =[]
thr =[]
four =[]

for j in range(len(r)):
    oone.append(Mat(r,4)[j][0][0])
    ttwo.append(Mat(r,4)[j][0][1])
    thr.append(Mat(r,4)[j][1][0])
    four.append(Mat(r,4)[j][1][1])


fig, ax = plt.subplots(figsize=(10,8))

plt.plot(r,oone)
plt.plot(r,ttwo)
plt.plot(r,thr)
plt.plot(r,four)
plt.show()