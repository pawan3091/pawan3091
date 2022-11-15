'''
date:-17/08/2022
Title:-Garrett Method to solve nested finite potential well
NAME=PAWANPREET KAUR
COLLEGE ROLL NO. 2020PHY1092
UNIVERSITY ROLL NO. 20068567038
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
h=6.626e-34 #Js
h_bar=1.05466e-34 #Js

def E_n(n,l,m):               #for the energy calculation
    E_0=((h)**2)/(8*m*(l**2))
    E_n=(n**2)*(E_0)
    return E_n

def delta_n(m,V,E):
    if E>V:
        return ("V should be greater than E.")
    else:
        Del_n=(1.05466e-34)/(np.sqrt(2*m*(V-E)))
    return Del_n

def L_n(l,del_n):            #New Length of the well
    l_n=l+(2*del_n)
    return l_n

def fun(n,l,m,V,tol,imax):
    l1=[];l3=[];l2=[]

    y1=E_n(n,l,m)
    y2=delta_n(m,V,y1)
    y3=L_n(l,y2)
    y4=E_n(n,y3,m)

    err=(abs(y4-y1))/y1

    l1.append(y4)
    l3.append(err)
    
    for i in range(imax):
        
        if err<tol:
            print("Tolerance achieved in",i+1,"iteration")
            break
        else: 
            
            del_new=delta_n(m,V,l1[-1])
            l_new=L_n(l,del_new)
            E_n_new=E_n(n,l_new,m)
            

            err= abs(E_n_new - l1[-1])/l1[-1]
            l2.append(del_new)
            l1.append(E_n_new)
            l3.append(err)
    #print(l1)
    #print(l3)
    return l1,l3

def table(n,l,m,V,tol,nmax):
    y,y1=fun(n,l,m,V,tol,nmax)
    y=np.array(y)/(1.6*10**(-19))
    data1 = {"E_N": y, "error": y1}
    print(pd.DataFrame(data1))
    lis1=[]
    for i in range(len(y)):
        lis1.append(i)

    plt.plot(lis1,y,label="E_n",ls="dashdot",marker="o")
    plt.title("E_n v/s N")
    plt.xlabel("N(iteration)")
    plt.ylabel("E_n")
    plt.grid()
    plt.legend()
    plt.show()
l=10*1e-10
m=(0.5*((10**6)*1.6e-19))/(9*(10**16))
V1=4*1.6e-19
tol=0.1e-9
nmax=1000

def new_func(v1,v2,l1,l2,n,m,imax,tol,key):
    E=0
    if l1!=l2:
        if v1!=v2:
            if key==1:
                E=E_n(n,l1,m)
                if E<v1:
                    ans=table(n,l1,m,v1,tol,imax)
                elif E>=v1:
                    print("l1 changes to l2")
                    E=E_n(n,l2,m)
                    if E>=v2:
                        print('Energy is greater then v2 so the particle is free')
                    elif v1<E<v2:
                        ans=table(n,l2,m,v2,tol,imax)
                    elif  E<v1:
                        sys.exit('invalid values of potentials')
            elif key==2:                      
                E=E_n(n,l2,m)
                if E>=v2:
                        print('Energy is greater then v2 so the particle is free')
                if v1<=E<v2:
                    ans=table(n,l2,m,v2,tol,imax)
                elif E<v1:
                    print("l2 changes to l1")
                    E=E_n(n,l1,m)
                    if E<v1:
                        ans=table(n,l1,m,v1,tol,imax)
                    elif  E>=v1:
                        sys.exit('invalid values of potentials')

        elif v1==v2:
            if l1==l2:
                ans=table(n,l2,m,v2,tol,imax)
            elif l1!=l2:
                ans=table(n,l1,m,v1,tol,imax)
    elif l1==l2:
        if v1==v2:
           ans=table(n,l2,m,v2,tol,imax)               
        elif v1!=v2:
            v=max(v1,v2)
            ans=table(n,l2,m,v,tol,imax)  


l1=4.9*1e-10
l2=5*1e-10
m=(0.5*((10**6)*1.6e-19))/(9*(10**16))
V1=6*1.6e-19
V2=10*1.6e-19
tol=0.1e-9
nmax=1000
new_func(V1,V2,l1,l2,1,m,nmax,tol,2)

