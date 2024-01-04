import numpy as np
import scipy.constants as ct
import random
import matplotlib.pyplot as plt

import json

k = np.linspace(1,10,500)
tau = np.linspace(0,1,500)


unit = 1e-21 

def bose_dist(x):

    T = 273
    boltz = ct.k*T

    return 1/(np.exp(x*unit/boltz)-1)

def green(tau):
    
    if tau[i] > 0:
        return (bose_dist(k)+1)*np.exp(-k*tau)
    if tau[i] < 0:
        return (bose_dist(k))*np.exp(-k*tau)
    
def omega(v):
    return v*k

def coupling(v,g,W):
    w = omega(v)
    cut_off = W
    return g*np.sqrt(w/(1+(w/cut_off)**2))

def interact(tau):

    file_path = "./data.json"

    with open(file_path, 'r') as file:
        data = json.load(file)
    
    v = data['k']
    g = data['g']
    W = data['W']

    g_k = np.abs(coupling(v,g,W)**2)

    n = len(k)

    k_sum = np.zeros(n)
    t_array = np.zeros(n)

    for j in range(n):
        for i in range(n):
            k_sum[i] = g_k[i] * green(tau)[i]
        t_array[j] = np.sum(k_sum)
        k_sum = np.zeros(len(k))
    
    #print(t_array)
    return t_array
