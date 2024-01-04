import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as ct

eps0 = ct.epsilon_0
mu0 = ct.mu_0
c0 = ct.speed_of_light
imp0 = np.sqrt(mu0/eps0)

jmax = 500
Hmax = 500
jsource = 10
nmax = 3000

#t=np.linspace(0,1,500)

Ex = np.zeros(jmax)
Hz = np.zeros(jmax)
Ex_prev = np.zeros(jmax)
Hz_prev = np.zeros(jmax)

#block to analyze response
Ex_sum = np.zeros(nmax)
response = []
response_x = np.arange(nmax)
frequency = np.ones(3000)

#frequency domain
for i in range(1,3000):
    frequency[i] = 2*np.pi/i


lambda_min = 350e-9
dx = lambda_min/20
dt = dx / c0

eps=np.full(500,eps0)

unit = ct.kgf
conduct = 1.55*unit

eps[35:60]= np.sqrt(np.sqrt(eps0**2 + 16*np.pi**2*conduct**2))*eps0

eps_m = np.tile(eps,5)

material_prof = eps_m > eps0

def Source_Function(t):
    lambda_0 = 550e-9
    w0 = 2*np.pi*c0/lambda_0
    tau = 30
    t0 = tau*3
    
    return np.exp(-(t-t0)**2/tau**2)*np.sin(w0*t*dt)

def Function(t):
    n = 10

    return (1/(t-0.9))*np.sin(n*(t))


for n in range(nmax):
    #update magnetic field boundaries
    Hz[jmax-1] = Hz_prev[jmax-2]
    Hz[0] = Hz[jmax-1]
    
    for j in range(jmax-1):
        Hz[j] = Hz_prev[j] + dt/(dx*mu0) * (Ex[j+1] - Ex[j])
        Hz_prev[j] = Hz[j]
    #Magenetic field source
    Hz[jsource-1] -= Source_Function(n)/imp0
    Hz_prev[jsource - 1] = Hz[jsource-1]

    #update magnetic field boundaries
    Ex[0] = Ex_prev[499] 
    #Brilloin zone boundary condition
    #Ex[499]=Ex[0]
    #Update electric field source
    for j in range(1,jmax):
        Ex[j] = Ex_prev[j] + dt/(dx*eps_m[j]) * (Hz[j] - Hz[j-1])
        Ex_prev[j] = Ex[j]
        #for calculate response
        Ex_sum[j] = Ex[j]
        
    response.append(Ex_sum)
    #Electric field source
    Ex[jsource] += Source_Function(n+1)
    Ex_prev[jsource] = Ex[jsource]
    fft_Ex=np.fft.fft(Ex)

    if n%10 == 0:
        plt.subplot(2,1,1)
        plt.plot(Ex)
        plt.plot(material_prof)
        plt.ylim([-1,1])

        plt.subplot(2,1,2)
        plt.plot(fft_Ex,color='b')
        plt.ylim([-1,1])

        plt.show()
        plt.close()

#response transform block

#for i in range(1,500):
    #response[:][i]=response[:][i-1]+response[:][i]
#n_response = np.array(response[:][498])

#plt.plot(n_response)
#plt.show()

#analyze frequency domain
#plt.subplot(2,1,1)
#plt.plot(response_x,response)
#plt.show()

#plt.subplot(2,1,2)
#plt.plot(frequency,np.fft.fft(response))
#plt.xlim(0.04,0.12)
#plt.show()