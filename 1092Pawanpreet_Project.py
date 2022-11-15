'''
PAWANPREET KAUR
College Roll No. =2020PHY1092
University Roll No. = 20068567038
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as integrate
import imageio.v2 as imageio


a=1 #arbitrarily taken 1
A=(2*a/np.pi)**(1/4) #constant A calculated using normalisation 

#function to calculate psi(x,0)
def psi_0(x1,a):
        psi=[]
        for x in x1:        
           psi.append(A*np.exp(-a*x*x))
        return x1,np.array(psi)                        

#function to calculate phi(k)
def phi_k(x,psi_,k):
    value=[]
    for i in k:
        value_=np.exp(-1j*i*x)
        m_v=psi_*value_/np.sqrt(2*np.pi)
        d=integrate.simps(m_v,x)
        value.append(d)
    a=integrate.simps(np.power(value,2),x)
    value1=(value)/(a)
    return k,np.array(value1)
    
#function to calculate phi(x,t)    
def time_evolution(phi_,x,t,k):
    psi_x_t=[]
    for i in range(len(x)):
        factor=np.exp((-1j*(k*x[i] -(k**2)*t)))
        value=np.array(phi_)*factor/np.sqrt(2*np.pi)
        d=integrate.simps(value,k)
        psi_x_t.append(d)
    norm=integrate.simps(np.power(psi_x_t,2),k)   
    psi_norm=psi_x_t/np.sqrt(norm)
    return psi_norm
    
N=1000
x1=np.linspace(-4,4,N)   
x,psi_=psi_0(x1,1)
k=np.linspace(-10,10,N)
k,phi_=phi_k(x,psi_,k)

#probability density at t=0
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Probability Density plot at t=0')
ax1.plot(x,np.power(psi_,2))   #position basis
ax1.grid()
ax1.set(xlabel="x",ylabel="$U(x)^2$",title="Probability Density in position basis")

ax2.plot(k, abs(phi_**2))  #momentum basis
ax2.grid()
ax2.set(xlabel="k",ylabel="$U(k)^2$",title="Probability Density in momentum basis")
plt.show()
#plt.savefig('Prob_density.png', dpi=72)
#plt.close()


t1=np.linspace(0,5,50)   #time frames
dpi = 72
s=np.arange(0,50,1)

#plotting psi(x,t) in different time frames
for (i,f) in zip(t1,s):
    m=time_evolution(phi_,x,i,k)
    fig = plt.figure(figsize=(600/dpi, 450/dpi), dpi=dpi)
    ax = fig.add_subplot(111)
    
    ax.plot(k, abs(np.power(m,2)), c='b', lw=3, alpha=0.8)
    ax.set_ylim(0, 0.4)
    ax.grid()
    # Fill under the line with reduced opacity
    ax.fill_between(k,abs(np.power(m,2)), facecolor='b', alpha=0.5)
    ax.yaxis.set_tick_params(width=2, length=5)
    ax.spines['left'].set_position('center')
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_linewidth(2)
    ax.annotate(text='$\ tau=${:.2f} '.format(i), xy=(0.8, 0.8),xycoords='axes fraction', ha='center', va='center')
    ax.set_xlabel('x')
    ax.set_ylabel("psi(x,t)")
    ax.set_label("NUMERICAL")
    plt.savefig('psi2-{0}.png'.format(f), dpi=dpi)
    plt.close()

# Creating the animation from png images
with imageio.get_writer('psi2_num.gif', mode='i') as writer:
    for i in s:
        image = imageio.imread('psi2-{0}.png'.format(i))
        writer.append_data(image)


#ANALYTIC SOLUTION
'''
The code below generates the frames in the above animation of |Ψ(x,t)|^2 for an electron with a=1bohr^(−1/2). We will work in atomic units so me=1 and hbar=1, but convert the time to attoseconds (as) for the annotation.
'''
# Grid of times t1 is in atomic units
hbar, Eh = 1.054571726e-34, 4.35974417e-18 # hbar and hartree in SI units for the time conversion
def plot_psi2(ax, i, t, psi2):
    # Plot |psi|^2 on Axes ax for frame i, time t
    ax.plot(x, psi2, c='r', lw=3, alpha=0.8)
    # Fill under the line with reduced opacity
    ax.fill_between(x, psi2, facecolor='r', alpha=0.5)
    ax.set_ylim(0, 0.8)
    ax.grid()
    ax.yaxis.set_tick_params(width=2, length=5)
    ax.spines['left'].set_position('center')
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['bottom'].set_linewidth(2)
    ax.xaxis.set_tick_params(width=2, length=5, direction='out')
    ax.yaxis.set_ticklabels([])

    # Add x-axis label and annotate with time in attoseconds
    t_in_as = t * hbar/Eh * 1.e18
    ax.annotate(text='{:.2f} as'.format(t_in_as), xy=(0.8, 0.8),xycoords='axes fraction', ha='center', va='center')
    ax.set_xlabel('$x$ / bohr')
    ax.set_label("ANALYTIC")

# Creating the animation frames at 72 dpi, 600x450 pixels as PNG images

for i, t in enumerate(t1):
    w = np.sqrt(a/(1+(2*a*t)**2))
    psi2 = np.sqrt(2/np.pi) * w * np.exp(-2*w*x**2)
    fig = plt.figure(figsize=(600/dpi, 450/dpi), dpi=dpi)
    ax = fig.add_subplot(111)
    plot_psi2(ax, i, t, psi2)
    plt.savefig('psi2_ana-{0}.png'.format(i), dpi=dpi)
    plt.close()
  
   
with imageio.get_writer('psi2_ana.gif', mode='i') as writer:
    for i in s:
        image = imageio.imread('psi2_ana-{0}.png'.format(i))
        writer.append_data(image)
