#postprocessing for diffusion problem

from boutdata import collect
#from boutdata import *
#from boututils import *

import os
import subprocess

import matplotlib.pyplot as plt
from matplotlib import *
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np
from numpy import *

Ni=collect('Ni', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
phi=collect('phi', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
rho=collect('rho', path="/home/colt/Research/BOUT/examples/RotationPlasma/data")
#print T.shape

#showdata(T[:,2,:,:])

#subprocess.Popen('rm *.png',shell=True)
#print T[50,33,2,:]
r=np.arange(0.98,1.02,0.01)
delta=2*pi/64
theta=np.linspace(0,2*pi,64)

[R,Theta]=meshgrid(r,theta)
"""
plt.figure()
plt.plot(r,Ni[0,:,2,10])
plt.grid(True)
plt.show()
"""

fig,ax=plt.subplots(subplot_kw=dict(projection='polar'))

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

CS=ax.contourf(Theta,R,rho[10,:,2,:].T,20)
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
