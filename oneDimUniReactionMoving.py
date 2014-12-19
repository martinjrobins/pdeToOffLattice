# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:05:32 2014

@author: robinsonm
"""

from fipy import *
from scipy.special import erf, erfc
from math import sqrt
import pyTyche as tyche
import numpy as np
import matplotlib.pyplot as plt

nx = 50
D = 1.
L = 1.
dx = L/nx
conversion_rate = 500.

interface = L/2.0
lam = 10.0**6
k = 10.
beta = sqrt(k/D)
N = lam/(2*beta*D)
threshold = 5*10**4

mol_dt = 10.0**(-4)
timeStepDuration = mol_dt*10.
steps = 200


#############
# PDE stuff #
#############

mesh = Grid1D(nx=nx, dx=dx)
phi = CellVariable(name="solution variable",  mesh=mesh, value=0.)
baseEq = TransientTerm() == DiffusionTerm(coeff=D) - ImplicitSourceTerm(coeff=k)

#####################
# Off-lattice stuff #
#####################

A = tyche.new_species([D,0,0])
grid = tyche.new_structured_grid([0,0,0],[L,1,1],[dx,1,1])
A.set_grid(grid)

xlow = tyche.new_xplane(0,1)
xhigh = tyche.new_xplane(L,-1)
xminboundary = tyche.new_reflective_boundary(xlow)
xmaxboundary = tyche.new_reflective_boundary(xhigh)

sink = tyche.new_uni_reaction(conversion_rate,[[A,A.lattice()],[A.pde()]])

source = tyche.new_zero_reaction_lattice(conversion_rate,[[A.pde()],[A]])


uni = tyche.new_uni_reaction(k,[[A],[]])

flux = tyche.new_zero_reaction(lam/dx,[0,0,0],[dx,1,1])

diffusion = tyche.new_diffusion()
algorithm = tyche.group([diffusion,flux,uni,sink,source,xminboundary])
algorithm.add_species(A)


############
# Plotting #
############

def concentration_gradient(x,t):
    exact = (lam/(D*beta)) * (
               np.exp(-beta*(x))
                - 0.5*np.exp(-beta*(x))*erfc((2.0*beta*D*t-(x))/np.sqrt(4.0*D*t))
                - 0.5*np.exp(beta*(x))*erfc((2.0*beta*D*t+(x))/np.sqrt(4.0*D*t))
        )
    return exact

plt.figure()
x = np.arange(dx/2.,L,dx)
t = 0
analytical = concentration_gradient(x,t)
plot_pde, = plt.plot(x,phi.value,linewidth=2,label='PDE')
off_lattice_concentration = A.get_concentration([0,0,0],[L,1,1],[nx,1,1])
plot_total, = plt.plot(x,phi.value+off_lattice_concentration[:,0,0],linewidth=2,label='Total')
plot_off_lattice = plt.bar(x-dx/2,off_lattice_concentration[:,0,0],width=dx)
plot_analytical, = plt.plot(x,analytical,linewidth=2,linestyle='--',label='Analytical')
plt.legend()
plt.ylim(0,N*1.5)


#############
# Time Loop #
#############
for step in range(steps):
    print "doing step ",step
    
    #plot
    plt.savefig("oneDimUniReactionMoving/oneDimUniReactionMoving%04d.png"%step)  
    
    mask = phi > threshold
    
    #set off-lattice generators 
    phiOld = phi.value*(!mask)
    A.set_pde(np.reshape(phiOld,[nx,1,1]))
    A.set_lattice(np.reshape(mask,[nx,1,1]))
    
    #integrate pde model
    eq = baseEq - conversion_rate*ImplicitSourceTerm(mask)
    eq.solve(var=phi, dt=timeStepDuration)    
    
    #integrate off-lattice model    
    t = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    #transfer sink particles to pde  
    phiNew = phi.value + A.get_pde()[:,0,0] - phiOld
    phi.setValue(phiNew)
    
    #update plotting
    analytical = concentration_gradient(x,t)
    plot_analytical.set_ydata(analytical)
    plot_pde.set_ydata(phiNew)
    off_lattice_concentration = A.get_concentration([0,0,0],[L,1,1],[nx,1,1])
    plot_total.set_ydata(phiNew+off_lattice_concentration[:,0,0])
    for rect, height in zip(plot_off_lattice, off_lattice_concentration[:,0,0]):
        rect.set_height(height)
    
    
