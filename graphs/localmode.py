#plot local modes with rotation frequency \Omega=-4

import matplotlib.pyplot as plt
from matplotlib import *
import numpy as np
from numpy import *

m=np.arange(2,11,1.0)

Omega=-4.0

plt.subplot(1,2,1)
plt.plot(m,m*Omega-2*Omega/m,'b',linewidth=1.5)
#plt.legend((r'$\omega$',r'$\gamma_-$',r'$\gamma_+$'),'best',shadow=True)
plt.title(r'Local modes with $\Omega=-4$',fontsize=18)
plt.xlabel(r'm',fontsize=18)
plt.ylabel(r'$\omega\left[\frac{c_s \rho_s}{a^2}\right]$',rotation='horizontal',fontsize=20)
plt.xlim([1,10])
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(m,(2*Omega/m)*sqrt(m**2/2-1),'r--',m,-(2*Omega/m)*sqrt(m**2/2-1),'r',linewidth=1.5)
plt.legend((r'$\gamma_-$',r'$\gamma_+$'),'best',shadow=True)
plt.title(r'Local modes with $\Omega=-4$',fontsize=18)
plt.xlabel(r'm',fontsize=15)
plt.ylabel(r'$\gamma\left[\frac{c_s \rho_s}{a^2}\right]$',rotation='horizontal',fontsize=20)
plt.grid(True)
plt.xlim([1,10])

plt.show()
