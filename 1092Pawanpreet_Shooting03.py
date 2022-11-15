'''
date:-28/09/2022
Title:-Shooting method to solve harmonic potential
NAME=PAWANPREET KAUR
COLLEGE ROLL NO. 2020PHY1092
UNIVERSITY ROLL NO. 20068567038
'''
import cmath
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from scipy.constants import hbar,m_e ,e#(Js),kg
const=9*10**9 #Nm**2/C**2
l0 = 50*10**-10
x = 200*10**-10
def RK4(func, X0,tmin,tmax,N,K):
    h=(tmax-tmin)/N
    t = np.linspace(tmin,tmax,N+1)
    X  = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        k1 = func(t[i],X[i],K[i])
        k2 = func( t[i] + h/2,X[i] + (h* k1)/2,K[i])
        k3 = func( t[i] + h/2,X[i] + h/2* k2,K[i])
        k4 = func(t[i] + h,X[i] + h   * k3,K[i])
        X[i+1] = X[i] + h / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X,t
    
def potential(n,V0):
    
    return V0*(e*e*const/n)
'''
    list=[]
    for i in n:

        if -x<=i<=-(l0):
             list.append(V0)
        elif -(l0)<i<(l0):
            list.append(0) 
        elif (l0)<=i<=x:
            list.append(V0)   
    return np.array(list)      
'''
def Sol1(x,Var,E,V):
    y1,y2=Var
    f1=y2
    f2=-((2*m_e*(E+V))/(hbar*hbar))
    return np.array([f1,f2])





def func(x,E,l0,m,N):
    f=x/l0
    epsilon_0 = hbar**2/(2*m*l0**2)
    print('epsilon',epsilon_0)
    eta1=np.linspace(-f,0, N)
    eta2=np.linspace(f, 0, N)

    V1=potential(eta1,V0)
    V2=potential(eta2,V0)
    small_e= E/epsilon_0
    v1 = V1/epsilon_0
    v2 = V2/epsilon_0
    
    kappa1=[]
    kappa2=[]
    omega1=[]
    omega2=[]
    k1=[]
    k2=[]

    for vi in v1:
    
        if vi>=small_e:
            kappa1.append(np.sqrt(vi-small_e))
            omega1.append(0)
            k1.append(np.sqrt(vi-small_e))

        else:
            kappa1.append(0)
            omega1.append(np.sqrt(small_e-vi))
            k1.append(np.sqrt(small_e-vi))

    for vi in v2:
     
        if vi>=small_e:
            kappa2.append(np.sqrt(vi-small_e))
            omega2.append(0)
            k2.append(np.sqrt(vi-small_e))

        else:
            kappa2.append(0)
            omega2.append(np.sqrt(small_e-vi))
            k2.append(np.sqrt(small_e-vi))
    

    data1={"eta1":eta1,"kappa1":kappa1,"omega1":omega1,"k1":k1}
    #print(pd.DataFrame(data1))

    data2={"eta2":eta1,"kappa2":kappa2,"omega2":omega2,"k2":k2}
    #print(pd.DataFrame(data2))
   
    
    #IC1=np.array([0,10])
    #IC2=np.array([0,10])
    IC1=np.real(np.array([0,1j*k1[0]*np.exp(1j*k1[0]*f)]))
    IC2=-np.real(np.array([0,-1j*k2[0]*np.exp(1j*k2[0]*(-f))]))

    res1=RK4(Sol1, IC2,0,-f,N,k1)
    res2=RK4(Sol1, IC1,0,f,N,k2)

    print('Initial Condition 1',IC1)
    print('Initial Condition 2',IC2)
    
    so1=res1[0].T
    #xs2=[-1]
    ls=so1[0][-1]
    ls_p=so1[1][-1]
    so2=res2[0].T
    rs=so2[0][-1]
    rs_p=so2[1][-1]
    #xs1=s1[-1]
    #print(np.array(k2)**2)
    return res1,res2,ls,rs,ls_p,rs_p,abs(ls-rs),ls_p-rs_p

def secant(e1,e2,tol):
    a=func(x,e1,l0,m,N)
    b=func(x,e2,l0,m,N)
    s1=(a[2]/a[4]-a[3]/a[5])
    s2=(b[2]/b[4]-b[3]/b[5])
    
    n=0
    while abs((e2-e1))>tol:
        e3=e2-(s1*(e2-e1)/(s2-s1))
        
        e1=e2
        e2=e3
        a=func(x,e1,l0,m,N)
        b=func(x,e2,l0,m,N)
        s1=abs(a[2]/a[4]-a[3]/a[5])
        s2=abs(b[2]/b[4]-b[3]/b[5])
        n+=1
        if n>1000:
            print('The itiration is too large')
            break
        else:
            continue
    print('The Root is:  ',e1)
    print('Number of Iterations:  ',n)
    return e1
'''
def secant(x1, x2, tol):
    n = 0; xm = 0; x0 = 0; c = 0;
    a=func(x,x1,l0,m,N)
    b=func(x,x2,l0,m,N)
    s1=abs(a[2]/a[4]-a[3]/a[5])
    s2=abs(b[2]/b[4]-b[3]/b[5])
    if (s1* s2< 0):
        while True:
             
            # calculate the intermediate value
            x0 = ((x1 * s2 - x2 * s1)/(s2- s1));
 
            # check ishoot x0 is root oshoot
            # equation or not
            c = s1 * s2;
 
            # update the value of interval
            x1 = x2;
            x2 = x0;
 
            # update number of iteration
            n += 1;
 
            # if x0 is the root of equation
            # then break the loop
            if (c == 0):
                break;
            xm = ((x1 * s2- x2 * s1)/(s2 - s1));
             
            if(abs(xm - x0) < tol):
                break;
         
        print("Root of the given equation =",round(x0, 6));
        print("No. of iterations = ", n);
    return x0
'''


V0=-5*e
m=m_e 

E1 =-1*e
E2 = 1*e

gg=x/l0
N=100
tol=0.000001
d=np.linspace(-gg,gg,100)
plt.plot(d,potential(d,V0)/e,label='Potential Function')

E=secant(E1,E2,tol)
print('E',E)
d1=func(x,E,l0,m,N)
s1,s2=d1[0],d1[1]

so1=s1[0].T
xs1=s1[-1]

so2=s2[0].T

xs2=s2[-1]


print('Phi from left at Origin',so1[0][-1])
print('Phi from right at Origin',so2[0][-1])

plt.plot(xs1,so1[0],label='Phi Left')
plt.plot(xs2,so2[0],label='Phi Right')
plt.xlabel("eta")
plt.legend()
plt.grid()
plt.show()

