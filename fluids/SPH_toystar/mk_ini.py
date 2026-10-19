'''
#!/usr/bin/env python
'''

# create initial particle distribution for toy star

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

R = 0.75
# number of particles, use N=100 for debugging, N=1000 for production run
N = 1000
M = 2.0

# EOS params (use the same in your SPH program for EOS)
n = 1.0
k = 0.1

# this is the lambda for the artificial gravitational force \vec{a} = -\lambda \vec{x}.
lambd = 2*k * (1+n) * np.pi**(-3/(2*n)) * (M * gamma(2.5+n)/(R**3*gamma(1+n)))**(1./n) * R**(-2)

# create sphere filled with arbitrary particles
# particle mass
mass = M/N

x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)
vx = np.zeros(N)
vy = np.zeros(N)
vz = np.zeros(N)
m = np.zeros(N)
i = 0
while i < N:
    phi = np.random.uniform(high=2.*np.pi)
    costheta = np.random.uniform(low=-1.0, high=1.0)
    u = np.random.uniform()
    theta = np.arccos(costheta)
    r = R * np.cbrt(u)

    x[i] = r * np.sin(theta) * np.cos(phi)
    y[i] = r * np.sin(theta) * np.sin(phi)
    z[i] = r * np.cos(theta)

    m[i] = mass
    i += 1


fig, ax = plt.subplots(1,3)
for _ax in ax:
    _ax.set_aspect('equal')
    _ax.set_xlim(-1,1)
    _ax.set_ylim(-1,1)

ax[0].scatter(x, y, s=1)
ax[1].scatter(x, z, s=1)
ax[2].scatter(y, z, s=1)
fig.show()


np.savetxt("random_distribution.dat", np.transpose([x,y,z,vx,vy,vz,m]))
