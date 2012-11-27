#generate 2-d cylindrical grid 5*16*65

from netCDF4 import Dataset

import numpy as np
from numpy import *

grid=Dataset("rotation_local.nc",'w')

x=grid.createDimension('x',20)
y=grid.createDimension('y',5)

grid.createVariable('nx','i4',())
grid.createVariable('ny','i4',())
grid.variables['nx'][:]=20
grid.variables['ny'][:]=5
nx=20
ny=5

dx=grid.createVariable('dx','f4',('x','y',))
dy=grid.createVariable('dy','f4',('x','y',))
dx[:,:]=0.005
dy[:,:]=0.01

#contravariant metric tensor

g11=grid.createVariable('g11','f4',('x','y',))
g22=grid.createVariable('g22','f4',('x','y',))
g33=grid.createVariable('g33','f4',('x','y',))
g12=grid.createVariable('g12','f4',('x','y',))
g13=grid.createVariable('g13','f4',('x','y',))
g23=grid.createVariable('g23','f4',('x','y',))

g12[:,:]=0.0
g13[:,:]=0.0
g23[:,:]=0.0

g11[:,:]=1.0
g22[:,:]=1.0

r=linspace(0.95,1.045,20)

for i in range(nx):
  g33[i,:]=r[i]**(-2)

#background rotation frequency
Omega=grid.createVariable('Omega','f4',())
grid.variables['Omega'][:]=-4.0
Omega=-4.0

#background potential
phi0=grid.createVariable('phi0','f4',('x','y',))

for i in range(nx):
  phi0[i,:]=2.0+Omega*r[i]**2/2

#background density
Ni0=grid.createVariable('Ni0','f4',('x','y',))

for i in range(nx):
  Ni0[i,:]=exp(-r[i]**2)

#background vorticity
rho0=grid.createVariable('rho0','f4',("x",'y',))
rho0[:,:]=Omega

grid.close()
