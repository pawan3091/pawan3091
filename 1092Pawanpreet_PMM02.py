#date:-24/08/2022
#Title:-Propegation matrix to solve Step Potential
#Name:-Pawanpreet KAur
#Roll no.:-2020PHY1092
#Reference https://phys.libretexts.org/Bookshelves/Quantum_Mechanics/Quantum_Mechanics_III_(Chong)/06%3A_Appendices/6.02%3A_B-_The_Transfer_Matrix_Method

import numpy as np
from scipy.constants import hbar,m_e ,e,h#(Js),kg
import sys
import matplotlib.pyplot as plt
import pandas as pd

def k(E,V):
    if V>E:
        sys.exit("E should be greater than or equal to V")
    else:
        pass
    return np.sqrt(2*m*(E-V))/hbar

def p_free(k_j,L_j):
    s=np.zeros((2,2),dtype='complex')
    s[0][0]=np.exp(1j*k_j*L_j)
    s[1][1]=np.exp(-1j*k_j*L_j)
    return s

def p_step(k_j1,k_j):
    s=np.zeros((2,2),dtype='complex')
    s[0][0]=1+(k_j/k_j1)
    s[0][1]=1-(k_j/k_j1)
    s[1][0]=1-(k_j/k_j1)
    s[1][1]=1+(k_j/k_j1)
    d=np.dot(0.5,s)
    return d


m=m_e #mass of particle
V0=0      
V=5*1.6*10**(-19) #J (height of barrier)
L=10e-10 #(m)
E=1.01*V #J
V_array=[V0,V,V0]
E_array=np.linspace(E,2*V,100) #array of E values

def main(E,V_array,L):
    k_V=k(E,V_array[1])
    k_V0=k(E,V_array[0])
    p1=p_step(k_V0,k_V)
    p2=p_free(k_V,L)
    p3=p_step(k_V,k_V0)
    
    p=np.matmul(p3,p2)
    p=np.matmul(p,p1)
  
    r=p[1][0]/p[1][1]
    t=np.linalg.det(p)/p[1][1]
    
    R=abs(r)**2
    T=abs(t)**2

    return  R,T



s=main(E,V_array,L)
print("V = 5eV , E=5.05 eV , L= 10 Angstrom")
print("Reflection Coefficient =",s[0])
print("Transmission Coefficient =",s[1])

c1=[] #ref coeff
c2=[] #  trans coef
for E in E_array:
    d=main(E,V_array,L)
    c1.append(d[0])
    c2.append(d[1])
    
plt.plot(c2,E_array/(1.6*10**(-19)),label="Transmission Coefficient")
plt.plot(c1,E_array/(1.6*10**(-19)),label="Reflection Coefficient")
plt.plot(np.array(c1)+np.array(c2),E_array/(1.6*10**(-19)),label="Transmission Coefficient+Reflection Coefficient")
plt.xlabel(" Coefficient")
plt.ylabel("Energy (in eV)")
plt.legend()
plt.grid()
plt.title("V = 5eV , E0=5.05 eV")
plt.show()

energy=E_array/(1.6*10**(-19))
f={"Energy (in eV)":energy,"R ":c1,"T ":c2,"R+T":np.array(c1)+np.array(c2)}
print(pd.DataFrame(f))

d=[]
la_1=[]
la_p=[]
o=[]
k_=[]
k_p=[]
for i,E_ in zip(c2,E_array):
    if i>=0.9995:
        k1=np.sqrt(2*m*E_)/hbar
        k1_p=np.sqrt(2*m*(E_-V))/hbar
        k_.append(k1)
        k_p.append(k1_p)
        d.append(i)
        o.append(E_)
        la_1.append(2*np.pi/k1)
        la_p.append(2*np.pi/k1_p)

x=np.array(o)/(1.6*10**(-19))
r1=np.array(la_1)/L
r2=np.array(la_p)/L
f={"Energy (in eV)":x,"T":d,"k":k_,"lambda":la_1,"k'":k1_p,"lambda'":la_p,"l/w":r1,"l'/w":r2}
print(pd.DataFrame(f))
    





