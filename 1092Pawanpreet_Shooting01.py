'''
#date:-7/09/2022
#Title:-Shooting Method to solve a constant Potential
NAME=PAWANPREET KAUR
COLLEGE ROLL NO. 2020PHY1092
UNIVERSITY ROLL NO. 20068567038
'''

import matplotlib.pyplot as plt 
import numpy as np 
import scipy.constants as const 
import pandas as pd 
from tabulate import tabulate 
import math

h_bar = const.hbar # J/sec
c = const.c # m/s
h = const.h # Js
l0 = 5e-10 # m
e = const.e # C
m = const.m_e # kg
a = 0
b = 20*1e-10
eeta_left = np.linspace(-b/l0,a/l0,500)
eeta_right = np.linspace(a/l0,b/l0,500)
eeta = np.concatenate((eeta_left,eeta_right))
V0 = 9*e 
E01 = -12*e
E02 = 17*e 
ε0 = h_bar*h_bar/(2*m*l0*l0)
v = V0/ε0
e1 = E01/ε0
e2 = E02/ε0




def numerov(phi,eeta,K):
		h = abs(eeta[0]-eeta[1])
		for i in range(1,len(eeta)-1):
			phi.append((2-(h*K)**2)*phi[i] - phi[i-1])
		return phi

def phi(K,eeta,case):
	if case=='l':
		return np.exp(complex(0,1)*K*eeta)
	elif case=='r':
		return np.exp(complex(0,-1)*K*eeta)

def shoot(e,case):
	if e>v:
		omega = np.sqrt(e-v)
		K=omega
	else :
		kappa = np.sqrt(v-e)
		K=complex(0,1)*kappa

	phi_initial_left = [np.exp(complex(0,1)*K*eeta_left[0]), np.exp(complex(0,1)*K*eeta_left[1])]
	phi_initial_right = [np.exp(complex(0,-1)*K*eeta_right[-1]), np.exp(complex(0,-1)*K*eeta_right[-2])]
	phi_estimate_left = numerov(phi_initial_left,eeta_left,K)
	phi_estimate_right = numerov(phi_initial_right,eeta_right,K)

	phi_estimate_left = np.array(phi_estimate_left)
	phi_estimate_right = np.array(phi_estimate_right)
	phi_estimate = np.concatenate((phi_estimate_left,np.flipud(phi_estimate_right)))
	if case==0:
		return eeta,phi_estimate
	elif case==1:
		return phi_estimate_right[0]-phi_estimate_left[-1]


def secant(x1, x2, tol):
    n = 0; xm = 0; x0 = 0; c = 0;
    if (shoot(x1,1) * shoot(x2,1) < 0):
        while True:
             
            # calculate the intermediate value
            x0 = ((x1 * shoot(x2,1) - x2 * shoot(x1,1))/(shoot(x2,1) - shoot(x1,1)));
 
            # check ishoot x0 is root oshoot
            # equation or not
            c = shoot(x1,1) * shoot(x0,1);
 
            # update the value of interval
            x1 = x2;
            x2 = x0;
 
            # update number of iteration
            n += 1;
 
            # if x0 is the root of equation
            # then break the loop
            if (c == 0):
                break;
            xm = ((x1 * shoot(x2,1) - x2 * shoot(x1,1))/(shoot(x2,1) - shoot(x1,1)));
             
            if(abs(xm - x0) < tol):
                break;
         
        print("Root of the given equation =",round(x0, 6));
        print("No. of iterations = ", n);
        x,y = shoot(x0,0)
        fig, ax = plt.subplots()
        ax.plot(x,y)
        ax.set_title("Wavefunction v/s Eeta")
        ax.set_ylabel("wavefunction")
        ax.set_xlabel('eeta')
        ax.grid() 
    else:
        print("Can not find a root in ","the given interval");

#eeta, phi_estimate = shoot(e)
secant(e1,e2,0.000001)
plt.show()
