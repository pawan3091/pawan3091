import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from scipy.integrate import simps
import pandas as pd

h_cut=1.0545*10e-34
m=3.32*10**(-27)


def potential(r,l):
    return (-2/(r) + (l*(l+1))/(r**2))

def v_effective(r,l):
    return (-2/(r) + (l*(l+1))/(r**2) + (h_cut**2)/(2*m)*((l*(l+1))/(r**2)))


def matrix(ri ,rf, N , l):
    h=(rf-ri)/(N + 1)
    X=np.linspace(ri,rf ,N + 2)
    x = X[1 : -1]
    a=np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i==j:
                a[i][j]=(-2 )* (-1/h**2) + potential(x[i] , l)
            if j==i+1:
                a[i][j]=-1/(h**2)
            if j==i-1:
                a[i][j]=-1/(h**2)
    return a ,x

eigen_function=[]

l=[1,2,3]
r=np.linspace(0.1,10,100)
for i in l:
    plt.plot(r,potential(r,i) ,label='Potential' )
    plt.plot(r,v_effective(r,i) , label='Effective Potential')
plt.grid()
plt.xlabel("r")
plt.ylabel("V(r)")
plt.legend()
plt.title("Potential Plot")
plt.show()


l=[0,1,2]
for j in l:
    analytical=[]
    eigen_value=[]
    for i in range (10):
        a,b=matrix(0,200,500,j)
        u, v = eigh(a)
        d=u[i]
        eigen_value.append(d)
        a=v[:,i]
        c=simps(np.power(a,2),b)
        u_n=(a)/np.sqrt((c))
        eigen_function.append(u_n)
        if j==1:
            analytical.append(-1/((i+2)**2))
        elif j==2:
            analytical.append(-1/((i+3)**2))
        else:
            analytical.append(-1/((i+1)**2))
    data={"Eigen Values " :eigen_value , "Analytical Values":analytical }
    print("Eigen Values for l =" + str(j))
    print(pd.DataFrame(data))
# print("Eigen Function" ,  eigen_function)
a,b=matrix(0,150,500,0)
u, v = eigh(a)
zero=[]
for i in range (4):
    d=u[i]
    eigen_value.append(d)
    a=v[:,i]
    c=simps(np.power(a,2),b)
    u_n=(a)/np.sqrt((c))
    zero.append(u_n[0])
    eigen_function.append(u_n)
    plt.plot(b, u_n,label='Radial Wavefunction for n = '+str(i))
    plt.xlabel("\u03BE")
    plt.legend()
    plt.ylabel("u(\u03BE)")
plt.title("Radial Waveform")
plt.grid()
plt.show()


a=v[:,0]
c=simps(np.power(a,2),b)
u_n=(a)/np.sqrt((c))
zero.append(u_n[0])
eigen_function.append(u_n)
plt.plot(b, u_n**2,label='Radial Wavefunction for n = '+str(0))
plt.xlabel("\u03BE")
plt.legend()
plt.ylabel("u(\u03BE)")
plt.title("Probability Density for different l = 0")
plt.grid()
plt.show()

# plt.plot(b,zero)
# plt.show()


a,b=matrix(0,150,500,1)
u, v = eigh(a)
zero=[]
for i in range (2):
    d=u[i]
    eigen_value.append(d)
    a=v[:,i]
    c=simps(np.power(a,2),b)
    u_n=(a)/np.sqrt((c))
    zero.append(u_n[0])
    eigen_function.append(u_n)
    plt.plot(b, u_n**2,label='Radial Wavefunction for n = '+str(i))
    plt.xlabel("\u03BE")
    plt.legend()
    plt.ylabel("u(\u03BE)")
plt.title("Probability Density for different l = 1")
plt.grid()
plt.show()



a,b=matrix(0,150,500,2)
u, v = eigh(a)
zero=[]
for i in range (3):
    d=u[i]
    eigen_value.append(d)
    a=v[:,i]
    c=simps(np.power(a,2),b)
    u_n=(a)/np.sqrt((c))
    zero.append(u_n[0])
    eigen_function.append(u_n)
    plt.plot(b, u_n**2,label='Radial Wavefunction for n = '+str(i))
    plt.xlabel("\u03BE")
    plt.legend()
    plt.ylabel("u(\u03BE)")
plt.title("Probability Density for different l = 2")
plt.grid()
plt.show()