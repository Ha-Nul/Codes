import numpy as np
import scipy as sp

import matplotlib.pyplot as plt

def coordvec(coord):
  qarr = np.array([coord[np.random.randint(0,len(coord))],coord[np.random.randint(0,len(coord))],coord[np.random.randint(0,len(coord))]])
  return qarr

def momenvec(momen):
  marr = np.array([momen[np.random.randint(0,len(momen))],momen[np.random.randint(0,len(momen))],momen[np.random.randint(0,len(momen))]])
  return marr

def Hamiltonian(x,y,mass,spconst):
    return x**2/(2*mass) + 0.5*spconst*y**2

def Hamconst(coord,momen,mass,spconst,Energy):
  qarr = coordvec(coord)
  marr = momenvec(momen)

  qnorm = np.linalg.norm(qarr)
  mnorm = np.linalg.norm(marr)
  #print(qnorm,mnorm)

  Hamcond = Hamiltonian(mnorm,qnorm,mass,spconst)
  #print(Hamcond)

  if Hamcond < (Energy + 0.1) and Hamcond > (Energy - 0.1):
    return (qnorm, mnorm, Hamcond)
    # Convert numpy arrays to tuples of floats
    #return (tuple(float(x) for x in qarr), tuple(float(x) for x in marr))
  else:
    None