import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from scipy.integrate import simps
import pandas as pd

def potential(r,l,ratio):
    return ((-2/(r)*np.exp(-r/ratio)) + (l*(l+1))/(2*r**2))
def potential_1(r):
    return -1/r
def matrix(ri ,rf, N , l,ratio):
    h=(rf-ri)/(N + 1)
    X=np.linspace(ri,rf ,N + 2)
    x = X[1 : -1]
    a=np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i==j:
                a[i][j]=(-2 )* (-1/h**2) + potential(x[i] , l ,ratio)
            if j==i+1:
                a[i][j]=-1/(h**2)
            if j==i-1:
                a[i][j]=-1/(h**2)
    return a ,x
eigen_function=[]

ratio=[2,5,10,20,100]
eigen_value0 = []
eigen_vector=[]
x=[]
for j in ratio:
    a,b=matrix(0,10,500,0,j)
    u, v = eigh(a)
    eigen_value0.append(u[0])
    eigen_vector.append(v[:,0])
    x.append(b)
    plt.scatter(b,v[:,0],label='for ratio ='+str(j))
    plt.legend()
plt.grid()
plt.xlabel("\u03BE")
plt.ylabel("u(\u03BE)")
plt.title("Radial Wavefunction ")
plt.show()

plt.scatter(ratio, eigen_value0 , label='Eigen value')
plt.grid()
plt.legend()
plt.xlabel("Ratio")
plt.ylabel("Eigen value of ground state")
plt.show()
data={"Ratio" :ratio , "Eigen Value" : eigen_value0}
print(pd.DataFrame(data))

for j in ratio:
    a,b=matrix(0,10,500,0,j)
    u, v = eigh(a)
    eigen_value0.append(u[0])
    eigen_vector.append(v[:,0])
    x.append(b)
    plt.scatter(b,np.power(v[:,0],2),label='Probability for ratio ='+str(j))
    plt.legend()
plt.grid()
plt.xlabel("\u03BE")
plt.ylabel("$u^2(\u03BE)$")
plt.title("Probability Radial Wavefunction ")
plt.show()

r=np.linspace(0.1,20,500)
for i in ratio:
    plt.plot(r,potential(r,0,i) ,label='Potential at ratio ='+str(i) )
    plt.plot(r,potential_1(r) ,label='Coulomb Potential at ratio ='+str(i))
plt.grid()
plt.xlabel("r")
plt.ylabel("V(r)")
plt.legend()
plt.title("Potential Plot")
plt.show()


