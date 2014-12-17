# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 00:18:14 2014

@author: mrobins
"""

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

valueLeft = 1
valueRight = 0

timeStepDuration = 9 * dx**2 / (2 * D)
steps = 200
mol_dt = timeStepDuration/100.

A = tyche.new_species([D,0,0])
A.fill_uniform([0,0,0],[L,0,0],1000)
grid = tyche.new_structured_grid([0,0,0],[L,0,0],[dx,0,0])
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
algorithm = tyche.group([diffusion, sink, xminboundary, xmaxboundary])
algorithm.add_species(A)

mesh = Grid1D(nx=nx, dx=dx)
phi = CellVariable(name="solution variable",  mesh=mesh, value=1000.)

# setup plotting
plt.figure()
x = np.arange(0,1,dx)
plot_pde, = plt.plot(x,phi.value,linewidth=2,label='PDE')
off_lattice_concentration = A.get_concentration([0,0,0],[L,1,1],[nx,1,1])
plot_off_lattice = plt.bar(x,off_lattice_concentration[:,0,0],width=dx)
plt.legend()
for step in range(steps):
    plt.savefig("test_plot%04d.pdf"%step)
    
    print "doing step ",step

    #set off-lattice generators 
    #(set B compartments equal to pde)
    A.set_pde(np.reshape(phi.value,[nx,1,1]))

    #integrate off-lattice model    
    time = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    plot_pde.set_ydata(phi.value)
    off_lattice_concentration = A.get_concentration([0,0,0],[L,0,0],[nx,1,1])
    for rect, height in zip(plot_off_lattice, off_lattice_concentration[:,0,0]):
        rect.set_height(height)
    
    

