#date:-28/09/2022
#Title:-Shooting method to Radial Equation of Hydrogen Atom
#Name:-Pawanpreet Kaur
#Roll no.:-2020PHY1092
import cmath
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.constants import hbar,m_e ,e#(Js),kg
k=9*10**9 #Nm**2/C**2
def RK4(func, X0,tmin,tmax,N,E):
    h=(tmax-tmin)/N
    t = np.linspace(tmin,tmax,N+1)
    X  = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        k1 = func(t[i],X[i],E)
        k2 = func( t[i] + h/2,X[i] + (h* k1)/2,E)
        k3 = func( t[i] + h/2,X[i] + h/2* k2,E)
        k4 = func(t[i] + h,X[i] + h   * k3,E)
        X[i+1] = X[i] + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X,t

def Sol1(x,Var,E):
    y1,y2=Var
    f1=y2
    f2=-((2*m_e*(E+(e*e*k/x)))/(hbar*hbar))
    return np.array([f1,f2])

IC=[0,1]
l=1*10**(-10)
r=5.29177210903e-11
E=-13.6*e

res1=RK4(Sol1,IC,0.1,l,100,E)
d1,d2=res1[0].T
d=np.sum(d1*d1)
r=np.arange(0.1*l,l,len(d1))
plt.scatter(res1[1],d1*d1/(d))
plt.xlabel("X (nm)")
plt.ylabel('Probability')
plt.title('Probability of finding the electron')
plt.grid()
plt.show()

