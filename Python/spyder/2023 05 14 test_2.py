import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as ct
from scipy import signal

mu0 = ct.mu_0
eps0 = ct.epsilon_0
c0 = ct.speed_of_light
imp0 = np.sqrt(mu0/eps0)

jmax = 500
jsource = 10
nmax = 3000

Ex = np.zeros(jmax)
Hz = np.zeros(jmax)

Ex_prev = np.zeros(jmax)
Hz_prev = np.zeros(jmax)

lambda_min = 350e-9

dx = lambda_min/20
dt = dx/c0

eps=np.ones(jmax)

unit = ct.kgf
conduct = 1.55*unit

for j in range(jmax):
        div = j//50
        a = div%2
        
        if a==1:
            eps[j] = np.sqrt(np.sqrt(eps0**2 + 16*np.pi**2*conduct**2))*eps0
        else:
            eps[j]=eps0

material_prof = eps > eps0

def Source_Function(t):
    lambda_0 = 550e-9
    w0 = 2*np.pi*c0/lambda_0
    t0 = w0*3
    
    return np.exp(-(t-t0)**2/w0**2)*np.sin(w0*t*dt)

def Function(t):
    n = 10

    return (1/t)*np.sin(n*(t))

for n in range(nmax):
    #boundary condition?
    Hz[jmax-1] = Hz_prev[jmax-2]

    for j in range(jmax-1):
        Hz[j] = Hz_prev[j] + dt/(dx*mu0)*(Ex[j+1]-Ex[j])
        Hz_prev[j] = Hz[j]
    fft_Hz=np.fft.fft(Hz)
    #boundary condition?
    Ex[0] = Ex_prev[1]

    for j in range(1,jmax):
        Ex[j] = Ex_prev[j] + dt/(dx*eps[j]) * (Hz[j] - Hz[j-1])
        Ex_prev[j] = Ex[j]
    fft_Ex=np.fft.fft(Ex)
    print(fft_Ex)

    Ex[jsource] += Function(n+1)
    Ex_prev[jsource] = Ex[jsource]

    if n%20 == 0:
        plt.subplot(2,1,1)
        plt.plot(Ex)
        plt.plot(material_prof)
        plt.ylim([-1,1])

        plt.subplot(2,1,2)
        plt.plot(fft_Ex)
        plt.ylim([-1,1])

        plt.show()
        plt.close()
    