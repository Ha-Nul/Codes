import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as ct

eps0 = ct.epsilon_0
mu0 = ct.mu_0
c0 = ct.speed_of_light
imp0 = np.sqrt(mu0/eps0)

jmax = 500
jsource = 2
nmax = 4000

Ex = np.zeros(jmax)
Hz = np.zeros(jmax)
Ex_prev = np.zeros(jmax)
Hz_prev = np.zeros(jmax)

lambda_min = 500e-9
dx = lambda_min/20
dt = dx / c0

eps=np.full(500,eps0)

eps[243:274]=4.243*eps0
#eps[149]=1
eps[275:315]=4.4873*eps0
eps[316:341]=36*eps0
eps[342:384]=84*eps0
eps[385:410]=36*eps0
eps[411:452]=84*eps0


#eps_m = np.tile(eps,5)
eps_m = eps
material_prof = eps_m > eps0

response = np.zeros(nmax)

def Source_Function(t):
    lambda_0 = 550e-9
    w0 = 2*np.pi*c0/lambda_0
    tau = 30
    t0 = tau*3
    
    return np.exp(-(t-t0)**2/tau**2)*np.sin(w0*t*dt)


for n in range(nmax):
    #update magnetic field boundaries
    Hz[jmax-1] = Hz_prev[jmax-2]
    #Hz[0] = Hz[jmax-1]
    
    for j in range(jmax-1):
        Hz[j] = Hz_prev[j] + dt/(dx*mu0) * (Ex[j+1] - Ex[j])
        Hz_prev[j] = Hz[j]
    #Magenetic field source
    Hz[jsource-1] -= Source_Function(n)/imp0
    Hz_prev[jsource - 1] = Hz[jsource -1]

    #update magnetic field boundaries
    Ex[0] = Ex_prev[1]
    #Ex[0] = Ex[499]
    #Update electric field source
    for j in range(1,jmax):
        Ex[j] = Ex_prev[j] + dt/(dx*eps_m[j]) * (Hz[j] - Hz[j-1])
        Ex_prev[j] = Ex[j]
    #Electric field source
    Ex[jsource] += Source_Function(n+1)/imp0
    Ex_prev[jsource] = Ex[jsource]
    
    fft_Ex=np.fft.fft(Ex)
    response[n]=Ex[250]

    if n%4000 == 0:
        plt.subplot(2,1,1)
        plt.plot(Ex,color='limegreen')
        plt.plot(material_prof,color='orange')
        plt.ylabel("Ez")
        plt.ylim([-1,1])
    
        plt.subplot(2,1,2)
        plt.plot(fft_Ex,color='b')
        plt.ylabel("Amplitute")
        plt.xlabel("Position")
        plt.plot()
        plt.ylim([-1,1])
    
        plt.show()
        plt.close()

resp_fft=np.fft.fft(response)
amp = abs(resp_fft)*(2/len(resp_fft))
freq = np.fft.fftfreq((len(resp_fft)),nmax)

#array for data
array_1 = []

index1 = np.where(amp >= 0.001)


for k in index1:
    array_1.append(amp[k])
    
#plt.subplot(2,1,1)
#plt.plot(response)
#plt.xlabel("time (a.u)")
#plt.ylabel("E (a.u)")


plt.plot(freq,amp)
plt.xlim(0.6e-5,1.6e-5)
plt.xlabel("frequency")
plt.ylabel("Amplitute")

plt.show()
plt.close()