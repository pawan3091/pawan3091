'''
date:-26/07/2022
Title:-Garrett Method to solve finite potential well
NAME=PAWANPREET KAUR
COLLEGE ROLL NO. 2020PHY1092
UNIVERSITY ROLL NO. 20068567038
'''


import numpy as np
from scipy.constants import hbar,m_e #(Js),kg
import sys
import matplotlib.pyplot as plt
import pandas as pd

def Energy(const,L):      #energy 
    return const/(L**2)

def delta(En): #penetration depth
    if V>En:     
         pass
    else:
        sys.exit("V should be greater than En")   #if V<En , we will get imaginary delta
    delta=hbar/np.sqrt(2*m*(V-En))
    return delta

def Ln(L,delta):  #length of infinite potential well
    return L+2*delta

def fun(En,V,m,L,n,tol):
    w=0         #setting up index to check whether tolerance has been achieved within given range of iterations
    i_array=[]      #array for storing number of iterations
    En_array=[En]    #array for storing value of energy
    delta_array=[]      #array for storing value of delta
    L_array=[]            #array for storing value of length
    ratio=[]      #array for storing value of ratio
    i=1
    while i<=10000:
          i_array.append(i)
          Delta=delta(En)
          delta_array.append(Delta)
          Lnew=Ln(L,Delta)          #calculating new length using L and delta
          L_array.append(Lnew)
          En=Energy(const,Lnew)     #calculating new energy using new length
          En_array.append(En)
          del_E=abs(En_array[-1]-En_array[-2])
          r=del_E/En_array[-2]      # calculating relative error
          ratio.append(r)
          i+=1
          if r<tol: #break if tolerance is achieved
            w=1
            break
          else:
            continue

    if w==1:
        message="tolerance achieved within given range of iterations"
    if w==0:
        message="tolerance not achieved within given range of iterations"

    return En_array,i_array,ratio,i,delta_array,L_array,message

def energy_plot(i_array,En_array,n):   #function for plotting En vs i for different n
    if n==1:
        c='r'
    elif n==2:
        c='green'
        
    elif n==3:
        c='magenta'
    elif n==4:
        c='purple'
    elif n==5:
        c='blue'
    elif n==6:
        c='brown'
    d=[0]+i_array
    plt.plot(d,np.array(En_array)/(1.6e-19),marker='o',c=c,label="$E_{0}$".format(n))

#input commands  
rest_mass_energy=float(input("Enter the rest mass energy of particle (in MeV) = "))
potential=float(input("Enter the potetial (in eV)= "))
len=float(input("Enter the length of potential barrier (in Angstrom)= "))

#m=rest_mass_energy*1.6e-19/((3e8)**2) #kg
#V=potential*10**6*1.6*10**(-19) #J
m=rest_mass_energy*10**6*1.6e-19/((3e8)**2) #kg
V=potential*1.6*10**(-19) #J
L=len*10**(-10) #(m)


#plotting En vs i for different values of n
for n in [1,2,3,4,5,6]:
    const=n*n*(np.pi**2)*hbar**2/(2*m)
    En=Energy(const,L)
    En_array,i_array,ratio,i,delta_array,L_array,message=fun(En,V,m,L,n,1e-9)
    energy_plot(i_array,En_array,n)

plt.ylabel("Energy $E_n$ (in eV)",c='r')
plt.xlabel("Iterations (i)",c='r')
plt.title("$E_n$ vs i plot",c='r')
plt.legend()
plt.grid()
plt.savefig("Figure1.png")
plt.show()

print()
print("#################################################################")
print()
print("FOR n=6   and i_max (maximum number of iterations)=10,000")

#table for printing iterations and ratio 
print("For n=6 " ,message)
s= {"iteraions":i_array,"Ratio":ratio}
print(pd.DataFrame(s))


#plotting  ratio vs iterations for n=6
plt.plot(i_array,ratio)
plt.xlabel("Number of iterations(i)",c='r')
plt.ylabel("Ratio $\dfrac{\Delta E}{E_i}$",c='r')
plt.title("$\dfrac{\Delta E}{E_i}$ vs i plot (n=6)",c='r')
plt.grid()
plt.savefig("Figure2.png")
plt.show()


r_array=[1e-5,1e-4,1e-3,1e-2,1e-1]
it=[]
msg=[]
for tol in r_array:
    En_array,i_array,ratio,i,delta_array,L_array,message=fun(En,V,m,L,n,tol)
    it.append(i)
    msg.append(message)

print()
print("#################################################################")
print()
print("For n=6 and i_max (Maximum number of iterations) = 10,000")
s= {"Ratio":r_array,"iteraions":it,"Message":msg}
print(pd.DataFrame(s))
print()
print("#################################################################")
print()

#plotting number of iterations vs tolerance 
plt.scatter(r_array,it)
plt.xlabel("tolerance",c='r')
plt.ylabel("Iterations (i)",c='r')
plt.title("Number of iterations vs tolerance plot (n=6)",c='r')
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.savefig("Figure3.png")
plt.show()

