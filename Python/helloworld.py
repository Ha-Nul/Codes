import numpy as np
import math
import matplotlib.pyplot as plt

beta = 1
t = 0.1
const = 0
N = 90

w_values = np.linspace(-10, -1, num=50)

a, b, c = 1, t+beta, t/beta+1

def hype_boy(w):
    result = 1
    for n in range(N):
        nume = math.prod([a + i for i in range(n)])*math.prod([b + i for i in range(n)])
        denom = math.prod([c + i for i in range(n)]) * math.factorial(n+1)
        result += (np.exp(beta*w)**(n+1))*nume/denom
    
    return result

def eqn(w):
    ans = (-1)*np.exp(w*t)*hype_boy(w)/t+const
    return ans 

result = [eqn(w) for w in w_values]



plt.plot(w_values, result)
plt.xlabel('w')
plt.ylabel('Equation Result')
plt.title('Graph of the Equation')
plt.grid(True)
plt.show()

    
