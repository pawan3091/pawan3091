'''
Pawanpreet kaur
2020PHY1092
'''

import numpy as np

def Eul(func, X0 ,tmin,tmax,N):
    h=(tmax-tmin)/N
    t = np.linspace(tmin,tmax,N+1)
    X  = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        X[i+1] = X[i] + func(t[i],X[i]) * h
    return X,t

def RK2(func, X0,tmin,tmax,N):
    h=(tmax-tmin)/(N)
    t = np.linspace(tmin,tmax,N+1)
    X  = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        k1 =h* func(t[i],X[i])
        k2 = h*func( t[i] + h,X[i] +  k1)
        X[i+1] = X[i] +  (k1 +k2 )/2
    return X,t
    
def RK4(func, X0,tmin,tmax,N):
    h=(tmax-tmin)/N
    t = np.linspace(tmin,tmax,N+1)
    X  = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        k1 = func(t[i],X[i])
        k2 = func( t[i] + h/2,X[i] + (h* k1)/2)
        k3 = func( t[i] + h/2,X[i] + h/2* k2)
        k4 = func(t[i] + h,X[i] + h   * k3)
        X[i+1] = X[i] + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X,t
    


