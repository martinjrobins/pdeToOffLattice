# -*- coding: utf-8 -*-
"""
Created on Wed Dec 24 12:07:07 2014

@author: robinsonm
"""

from fipy import *
from scipy.special import erf, erfc
from math import sqrt
import pyTyche as tyche
import numpy as np
import matplotlib.pyplot as plt

nx = 51
D = 1.
L = 1.
dx = L/nx
conversion_rate = 500.

lam = 10.0**6
k = 0.1
beta = sqrt(k/D)
N = lam/(2*beta*D)
threshold = 5*10**4

mol_dt = 10.0**(-4)
timeStepDuration = mol_dt*5.
steps = 200

#############
# PDE stuff #
#############

mesh = Grid2D(nx=nx,ny=nx)
phi = CellVariable(name="solution variable",  mesh=mesh, value=0.)
total = CellVariable(name="total variable",  mesh=mesh, value=0.)
mask = CellVariable(name="mask",  mesh=mesh, value=0.)

baseEq = TransientTerm() == DiffusionTerm(coeff=D) - ImplicitSourceTerm(coeff=k)

#####################
# Off-lattice stuff #
#####################

A = tyche.new_species([D,D,0])
dummy = tyche.new_species([0,0,0])
not_dummy = tyche.new_species([0,0,0])

grid = tyche.new_structured_grid([0,0,0],[L,L,1],[dx,dx,1])
A.set_grid(grid)
dummy.set_grid(grid)
not_dummy.set_grid(grid)

sink = tyche.new_uni_reaction(conversion_rate,[[A,dummy.pde()],[A.pde()]])
source = tyche.new_zero_reaction_lattice(conversion_rate,[[A.pde(),not_dummy.pde()],[A]])
uni = tyche.new_uni_reaction(k,[[A],[]])
flux = tyche.new_zero_reaction(2*pi*lam/(dx*dx),[L/2-dx/2,L/2-dx/2,0],[L/2+dx/2,L/2+dx/2,1])

diffusion = tyche.new_diffusion()
algorithm = tyche.group([diffusion,flux,uni,sink,source])
algorithm.add_species(A)


############
# Plotting #
############

def concentration_gradient(x,y,t):
    r = np.sqrt((x-L/2.)**2+(y-L/2.)**2)
    exact = (lam/(D*beta)) * (
               np.exp(-beta*(r))
                - 0.5*np.exp(-beta*(r))*erfc((2.0*beta*D*t-(r))/np.sqrt(4.0*D*t))
                - 0.5*np.exp(beta*(r))*erfc((2.0*beta*D*t+(r))/np.sqrt(4.0*D*t))
        )
    return exact

plt.figure()
cv_indicies = range((nx*nx-1)/2-(nx-1)/2,(nx*nx-1)/2+(nx-1)/2+1)
x = np.arange(0,L,dx)+dx/2
print x
xv, yv = np.meshgrid(x, x)
t = 0
analytical_2d = concentration_gradient(xv,yv,t)

plt.subplot(3, 1, 1)

plot_pde, = plt.plot(x,phi.value[cv_indicies],linewidth=2,label='PDE')
off_lattice_concentration = A.get_concentration([0,0,0],[L,L,1],[nx,nx,1])
total.setValue(phi.value+numerix.reshape(off_lattice_concentration[:,:,0],[nx*nx]))
plot_total, = plt.plot(x,total.value[cv_indicies],linewidth=2,label='Total')
plot_off_lattice = plt.bar(x-dx/2,off_lattice_concentration[(nx-1)/2,:,0],width=dx)
plot_analytical, = plt.plot(x,analytical_2d[(nx-1)/2,:],linewidth=2,linestyle='--',label='Analytical')
plt.legend()
plt.ylim(0,N*1.5)

plt.subplot(3, 1, 2)
plot_analytical_2d = plt.imshow(analytical_2d, interpolation='nearest', 
                            origin='bottom', 
                            vmin=0,
                            vmax=1.5*N, 
                            cmap='jet')

plt.colorbar(plot_analytical_2d)


plt.subplot(3, 1, 3)
plot_total_2d = plt.imshow(numerix.reshape(total.value,[nx,nx]), interpolation='nearest', 
                            origin='bottom', 
                            vmin=0,
                            vmax=1.5*N,  
                            cmap='jet')

plt.colorbar(plot_total_2d)


#############
# Time Loop #
#############
for step in range(steps):
    print "doing step ",step
    
    #plot
    plt.savefig("twoDimUniReactionMoving/twoDimUniReactionMoving%04d.png"%step)  
    
    mask = (total > threshold)
    mask = total > 0
    not_mask = mask==False
    
    #set off-lattice generators 
    phiOld = (phi.value*not_mask).value + 0.0
    A.set_pde(numerix.reshape(phiOld,[nx,nx,1]))
    dummy.set_pde(numerix.reshape((mask*1.0).value,[nx,nx,1]))
    not_dummy.set_pde(numerix.reshape((not_mask*1.0).value,[nx,nx,1]))

    #integrate pde model
    eq = baseEq + conversion_rate*ImplicitSourceTerm(coeff=mask==False)
    eq.solve(var=phi, dt=timeStepDuration)    
    
    #integrate off-lattice model
    t = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    #transfer sink particles to pde 
    phiNew = phi.value + np.reshape(A.get_pde()[:,:,0],[nx*nx]) - phiOld
    phi.setValue(phiNew)
    
    #update plotting
    analytical_2d = concentration_gradient(xv,yv,t)
    plot_analytical_2d.set_data(analytical_2d)
    plot_analytical.set_ydata(analytical_2d[(nx-1)/2,:])
    plot_pde.set_ydata(phiNew[cv_indicies])
    off_lattice_concentration = A.get_concentration([0,0,0],[L,L,1],[nx,nx,1])
    total.setValue(phi.value+np.reshape(off_lattice_concentration[:,:,0],[nx*nx]))
    
    plot_total_2d.set_data(numerix.reshape(total.value,[nx,nx]))
    plot_total.set_ydata(total.value[cv_indicies])
    for rect, height in zip(plot_off_lattice, off_lattice_concentration[(nx-1)/2,:,0]):
        rect.set_height(height)
    
    

