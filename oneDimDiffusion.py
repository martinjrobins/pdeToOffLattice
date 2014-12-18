# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:05:32 2014

@author: robinsonm
"""

from fipy import *
from scipy.special import erf
from math import sqrt
import pyTyche as tyche
import numpy as np
import matplotlib.pyplot as plt

nx = 50
D = 1.
L = 1.
dx = L/nx
conversion_rate = 100.

interface = L/2.0
N = 1000.
valueLeft = N
valueRight = 0

timeStepDuration = 9 * dx**2 / (2 * D)
steps = 200
mol_dt = timeStepDuration/100.

############
# PDE stuff#
############

mesh = Grid1D(nx=nx, dx=dx)

phi = CellVariable(name="solution variable",  mesh=mesh, value=N)
#phi.constrain(valueRight, mesh.facesRight)
#phi.constrain(valueLeft, mesh.facesLeft)

eq = TransientTerm() == DiffusionTerm(coeff=D) - conversion_rate*ImplicitSourceTerm((mesh.x > interface))

####################
# Off-lattice stuff#
####################

A = tyche.new_species([D,0,0])
#A.fill_uniform([0,0,0],[L,0,0],1000)
grid = tyche.new_structured_grid([0,0,0],[L,1,1],[dx,1,1])
A.set_grid(grid)

xlow = tyche.new_xplane(0,1)
xhigh = tyche.new_xplane(L,-1)
xminboundary = tyche.new_reflective_boundary(xlow)
xmaxboundary = tyche.new_reflective_boundary(xhigh)

#set off-lattice sink and source
sink = tyche.new_uni_reaction(conversion_rate,[[A],[A.pde()]])
pde_region = tyche.new_xplane(interface,1)
sink.set_geometry(0,pde_region)
source = tyche.new_zero_reaction_lattice(conversion_rate,[[A.pde()],[A]])
offlattice_region = tyche.new_xplane(interface,-1)
source.set_geometry(offlattice_region)

diffusion = tyche.new_diffusion()
algorithm = tyche.group([diffusion,sink,source,xminboundary,xmaxboundary])
algorithm.add_species(A)


############
# Plotting #
############

plt.figure()
x = np.arange(dx/2.,L,dx)
t = 0
analytical = N*(1 - erf(x / (2 * sqrt(D * t))))
plot_pde, = plt.plot(x,phi.value,linewidth=2,label='PDE')
plot_analytical, = plt.plot(x,analytical,linewidth=2,label='Analytical')
off_lattice_concentration = A.get_concentration([0,0,0],[L,1,1],[nx,1,1])
plot_total, = plt.plot(x,phi.value+off_lattice_concentration[:,0,0],linewidth=2,label='Total')
plot_off_lattice = plt.bar(x-dx/2,off_lattice_concentration[:,0,0],width=dx)
plt.legend()
plt.ylim(0,1200)

for step in range(steps):
    print "doing step ",step
    print A
    plt.savefig("oneDimDiffusion%04d.png"%step)    
    
    
    #set off-lattice generators 
    A.set_pde(np.reshape(phi.value,[nx,1,1]))
    phiOld = phi.copy()

    #integrate pde model    
    eq.solve(var=phi, dt=timeStepDuration)    
    
    #integrate off-lattice model    
    t = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    #transfer sink particles to pde  
    print A.get_pde()[:,0,0] - phiOld
    phi.setValue(phi.value + A.get_pde()[:,0,0] - phiOld)
    
    analytical = N*(1 - erf(x / (2 * sqrt(D * t))))
    plot_analytical.set_ydata(analytical)
    plot_pde.set_ydata(phi.value)
    off_lattice_concentration = A.get_concentration([0,0,0],[L,1,1],[nx,1,1])
    plot_total.set_ydata(phi.value+off_lattice_concentration[:,0,0])
    for rect, height in zip(plot_off_lattice, off_lattice_concentration[:,0,0]):
        rect.set_height(height)
    
    
