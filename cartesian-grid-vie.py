# FFT-accelerated VIE solver using a Cartesian grid. 
# Currently using "DDA" evaluation of all the integrals.

import numpy as np
from scipy.special import hankel1

ko = 20  # Wavenumber
rad = 1. # radius of circle
lam = 2*np.pi/ko
refInd = 1.3
n_per_lam = 20  # Pixels per wavelength

h_temp = lam / n_per_lam # temp pixel dimension

wx = 2 * rad
wy = 2 * rad

# How many points in x and y directions
M = np.int(np.ceil(wx / h_temp))
N = np.int(np.ceil(wy / h_temp))

dx = wx/M
dy = wy/N

A = dx * dy      # pixel area
a = np.sqrt(A / np.pi) # radius of equivalent-area circle

# Get coordinates of points on grid    
# FIX ME: I thought this complex number stuff was elegant at first, but it's
# actually just annoying. Switch to using meshgrid
x = np.zeros((M*N, 1), dtype=np.complex128)
counter= 0
for j in range(N):
    for i in range(M):
        x[counter] = -wx/2 + dx/2+dx*i \
            - 1j*wy/2 + 1j * (dy/2+dy*j)
        counter=counter+1

x_coord = (np.arange(M)+1) * dx - dx/2 - wx/2
y_coord = (np.arange(N)+1) * dy - dy/2 - wy/2

perm = np.ones(M*N)  # permittivities

# FIX ME: include more geoemetries. Just circle at the moment.
idx = np.where(np.abs(x)<=rad) # locate indices of points inside circle