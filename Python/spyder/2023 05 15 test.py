import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(-20,20,100)
y=np.zeros(100)

for i in range(100):
    y[i]=x[i]**2
    
    
    if i <= 100:
    
        plt.plot(x,y)
        plt.show()
        plt.close()