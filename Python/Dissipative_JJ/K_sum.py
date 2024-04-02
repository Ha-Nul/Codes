import numpy as np
import scipy as sp
import cmath
import random

np.set_printoptions(threshold=784,linewidth=np.inf)

## Hamiltonian source code available : 2023 09 06.ipynb
## pre sourcecode available : 2023 05 18 function.ipynb

## K sum prerequisites ####################

## k linspace and tau grid requried #################

#unit = plank const

mode_grid = np.linspace(1,100,100)

OMG_Arr = np.zeros(len(mode_grid)) # not shibainu
g_Arr = np.zeros(len(mode_grid))

def Tilde_g_cal_function(alpha: float, k_cutoff: float, mode_grid):
    '''Calculates the g value for coupling.'''
    nu = np.pi * k_cutoff / alpha

    for i in range(len(mode_grid)):
        OMG_Arr[i]  = (k_cutoff * mode_grid[i]/mode_grid[len(mode_grid)-1])
        g_Arr[i] = (np.sqrt((2 * k_cutoff / (alpha * mode_grid[len(mode_grid)-1])) * (OMG_Arr[i] / (1 + (nu * OMG_Arr[i]/k_cutoff)**2))))

def Interact_V(mode_grid, tau_grid):
    INT_Arr = np.zeros(len(tau_grid))

    for i in range(len(tau_grid)):
        for j in range(len(mode_grid)):
            INT_Arr[i] += g_Arr[j]**2 * np.cosh((tau_grid[i] - tau_grid[len(tau_grid)-1]/2) * OMG_Arr[j]) / np.sinh(tau_grid[len(tau_grid)-1] * OMG_Arr[j] / 2)
    
    return INT_Arr

def Tilde_g(alpha,k_cutoff,mode_grid):
    Tilde_g_cal_function(alpha,k_cutoff,mode_grid)
    return np.sum(g_Arr)

