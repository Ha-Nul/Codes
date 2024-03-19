import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import cmath
import random

np.set_printoptions(threshold=784,linewidth=np.inf)

## Full source code available : 2023 09 06.ipynb

## K sum prerequisites ####################

## k linspace and tau grid requried #################

#unit = 1e-21 

def bose_dist(x,b):

    T = 273
    #boltz = ct.k*Ts

    return 1/(np.exp(x*b)-1)

def green(tau,k,b):

    #for i in range(len(tau)):
        #if tau[i] > 0:
            #return (bose_dist(k)+1)*np.exp(-k*tau)
        #if tau[i] < 0:
            #return (bose_dist(k))*np.exp(-k*tau)
        return ((bose_dist(k,b)+1)*np.exp(-k*tau)) + (bose_dist(k,b))*np.exp(k*tau)

def omega(v,k):
    return v*np.abs(k)

def coupling(v,g,W,k):
    w = omega(v,k)
    cut_off = W
    return g*np.sqrt(w/(1+(w/cut_off)**2))

def interact(tau,v,g,W,k,b):
    g_k = np.abs(coupling(v,g,W,k))**2

    n = len(k)

    k_sum = np.zeros(n)
    t_array = np.zeros(n)

    for j in range(n):
        t = tau[j]
        for i in range(n):
            k_sum[i] = g_k[i] * green(t,omega(v,k),b)[i]
        t_array[j] = k_sum[j]#np.sum(k_sum)
        k_sum = np.zeros(len(k))
    
    #print(t_array)
    return t_array

