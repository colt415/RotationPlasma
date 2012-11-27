#postprocessing for diffusion problem

from boutdata import collect
from boutdata import *
from boututils import *

import os
import subprocess

import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np
from numpy import *
from scipy import fft

Ni0=collect('Ni0', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
Ni=collect('Ni', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
phi0=collect('phi0', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
phi=collect('phi', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
rho=collect('rho', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
#print T.shape

#showdata(rho[:,:,1,:])

#subprocess.Popen('rm *.png',shell=True)
#print T[50,33,2,:]

r=linspace(0.95,1.045,20)
delta=2*pi/64
theta=np.linspace(0,2*pi,64)
[R,Theta]=meshgrid(r,theta)

t=linspace(0,201*0.02,201)

n=200
Ts=0.05
Fs=1/Ts
deltaf=Fs/n
fre=np.arange(0,(n/2)*deltaf,deltaf)
print fre.shape
FFT=fft(Ni[:,10,2,0])
print FFT.shape
#plt.plot(r,Ni0[:,2])
plt.plot(t,Ni[:,10,2,0])
plt.grid(True)
#plt.subplot(1,2,2)
#plt.plot(fre,FFT[n/2:n])
#plt.grid(True)
plt.show()


r1=np.arange(0.05,0.95,0.05)
theta1=np.arange(0,2*pi+pi/32,pi/32.0)
[R1,Theta1]=meshgrid(r1,theta1)
A=zeros([len(theta1),len(r1)])

fig,ax=plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
CS=ax.contourf(Theta,R,rho[0,:,2,:].T,20)
CS1=ax.contourf(Theta1,R1,A[:,:])
CB=fig.colorbar(CS,shrink=0.8,extend='both')
plt.title(r'Diffusion with Gaussian initial condition, $A=1.0, \sigma=0.32$',fontsize=15)
plt.show()

"""
pp=PdfPages('.pdf')

for i in range(101):
  number=str(i+1).zfill(3)
  fig,ax=plt.subplots(subplot_kw=dict(projection='polar'))

  ax.set_theta_zero_location("N")
  ax.set_theta_direction(-1)

  CS=ax.contourf(Theta,R,T[i,:,2,:].T,20)
  CB=fig.colorbar(CS,shrink=0.8,extend='both')
  plt.title(r'Diffusion with Gaussian initial condition, $A=1.0, \sigma=0.32  $'+ number,fontsize=15)
  plt.savefig(pp,format='pdf')

pp.close()
"""
